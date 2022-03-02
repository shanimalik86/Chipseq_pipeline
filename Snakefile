import re

configfile: "config.yaml"

DATA_DIR = config["data"]

if config["type"]=="SE":
	SAMPLES, = glob_wildcards(DATA_DIR + "/{sample,[^(/|INPUT)]+}.fastq.gz")
elif config["type"]=="PE":
        SAMPLES, = glob_wildcards(DATA_DIR + "/{sample,[^/]+}_R1.fastq.gz")
else:
     	raise ValueError('please specify only "SE" or "PE" for the "type" parameter in the config file.')


def macs_input(wildcards):
	names=wildcards.sample.split("_")
	return "filtered/" + names[0] + "_INPUT_" + names[2] + "_rmdup.bam"

rule all:
	input:
              	bam = expand("macs2/{sample}_peaks.xls", sample=SAMPLES)

if config["type"]=="SE":
        rule trim:
                input:
                      	DATA_DIR + "/{sample}.fastq.gz"
                output:
                       	"trimmomatic_fastq/{sample}.trimmed.fastq.gz"
                message:
                        "I am trimming"
                params:
                       	length=config["trim"]["length"],
                        trail=config["trim"]["trail"],
                        lead=config["trim"]["lead"],
                        window=config["trim"]["window"],
                        type=config["type"]
                shell:
                      	"""
                        mkdir -p trimmomatic_fastq
                        module load Trimmomatic
                        trimmomatic {params.type} -threads 2 -phred33 {input} {output} ILLUMINACLIP:/home/apps/software/Trimmomatic/0.38-Java-1.8.0_152/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:{params.window} LEADING:{params.lead} TRAILING:{params.trail} MINLEN:{params.length}
                        """
else:
     	rule trim:
                input:
                      	r1=DATA_DIR + "/{sample}_R1.fastq.gz",
                        r2=DATA_DIR + "/{sample}_R2.fastq.gz"
                output:
                       	TR1="trimmomatic_fastq/{sample}_R1.trimmed.fastq.gz",
                        TR2="trimmomatic_fastq/{sample}_R2.trimmed.fastq.gz",
                        TR1un="trimmomatic_fastq/{sample}_R1un.trimmed.fastq.gz",
                        TR2un="trimmomatic_fastq/{sample}_R2un.trimmed.fastq.gz"
                message:
                        "I am paired trimming"
                params:
                       	length=config["trim"]["length"],
                        trail=config["trim"]["trail"],
                        lead=config["trim"]["lead"],
                        window=config["trim"]["window"],
                        type=config["type"]
                shell:
                      	"""
                        mkdir -p trimmomatic_fastq
                        module load Trimmomatic
                        trimmomatic {params.type} -threads 2 -phred33 {input.r1} {input.r2} {output.TR1} {output.TR1un} {output.TR2} {output.TR2un} ILLUMINACLIP:/home/apps/software/Trimmomatic/0.38-Java-1.8.0_152/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:{params.window} LEADING:{params.lead} TRAILING:{params.trail} MINLEN:{params.length}
                        """

rule bowtie2_index:
	input:
		fasta = config["fasta"]
	output:
		touch("index.done")
	params:
		prefix = "genome"
	shell:
		"""
		module load Bowtie2
		bowtie2-build {input.fasta} {params.prefix}
		"""

rule bowtie2_mapping:
	input:
		trimfile = "trimmomatic_fastq/{sample}.trimmed.fastq.gz",
		index = rules.bowtie2_index.output
	output:
		alignment = "bowtie/{sample}.sam"
	message:
		"I am Bowtie2 Alignment"
	params:
		index = "genome"
	shell:
		"""
		mkdir -p bowtie
		module load Bowtie2
		bowtie2 -p 12 -q -x {params.index} -U {input.trimfile} -S {output.alignment}
		"""

rule bamconversion:
	input:
		"bowtie/{sample}.sam"
	output:
		"filtered/{sample}_unsorted.bam"
	shell:
		"""
		mkdir -p filtered
		module load SAMtools
		module load BEDTools
		samtools view -h -S -b {input} -o {output}
		"""
rule bam_sort:
	input:
		"filtered/{sample}_unsorted.bam"
	output:
		"filtered/{sample}_sorted.bam"
	shell:
		"""
		module load sambamba
		sambamba sort -t 12 -o {output} {input}
		"""

rule blacklist:
	input:
		bam = "filtered/{sample}_sorted.bam",
		black = config["black"]
	output:
		"filtered/{sample}_blacklisted.bam"
	shell:
		"""
		module load BEDTools
		bedtools intersect -v -abam {input.bam} -b {input.black} > {output}
		"""

rule filter:
	input:
		"filtered/{sample}_blacklisted.bam"
	output:
		"filtered/{sample}_unique.bam"
	shell:
		"""
		module load sambamba
		sambamba view -h --nthreads 12 -f bam -F "[XS] == null and not unmapped" {input} > {output}
		"""

rule remove_dup:
	input:
		"filtered/{sample}_unique.bam"
	output:
		"filtered/{sample}_rmdup.bam"
	shell:
		"""
		module load SAMtools
		samtools rmdup -s {input} {output}
		"""
rule bam_index:
	input:
		"filtered/{sample}_rmdup.bam"
	output:
		"filtered/{sample}_rmdup.bam.bai"
	shell:
		"""
		module load SAMtools
		samtools index {input}
		"""

rule callpeaks:
	input:
		t="filtered/{sample}_rmdup.bam",
		c=macs_input
	output:
		"macs2/{sample}_peaks.xls"
	shell:
		"""
		mkdir -p macs2
		module load MACS2/2.2.5-IGB-gcc-4.9.4-Python-3.6.1
		macs2 callpeak -t {input.t} -c {input.c} -f BAM --gsize 3.0e9 --broad --broad-cutoff 0.1  --nomodel --extsize 125 -n $prefix --bdg --outdir macs2
		"""
