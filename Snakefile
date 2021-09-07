import os
import re
#import pandas as pd
#======================================================
# Config files
#======================================================
configfile: "config.yaml"

#======================================================
# Global variables
#======================================================

RAW_DATA_DIR =config['input_dir']
RESULTS_DIR=config['results_dir'].rstrip("/")
LONG_ASSEMBLER=config['long_assembler']

if RESULTS_DIR == "" and not RAW_DATA_DIR == "":
	RESULTS_DIR=os.path.abspath(os.path.join(RAW_DATA_DIR, os.pardir))


REPRESENTATIVE_CONTIGS=config['representative_contigs'].rstrip("/")
VIRAL_CONTIGS=REPRESENTATIVE_CONTIGS
if VIRAL_CONTIGS == "":
	VIRAL_CONTIGS_BASE="positive_contigs"
	VIRAL_CONTIGS_DIR=RESULTS_DIR + "/04_VIRAL_ID"
	REPRESENTATIVE_CONTIGS_BASE="95-85_positive_contigs"
	REPRESENTATIVE_CONTIGS_DIR=RESULTS_DIR + "/05_vOTUs"
else:
	REPRESENTATIVE_CONTIGS_BASE=os.path.basename(os.path.abspath(VIRAL_CONTIGS)).split(".")[0]
	REPRESENTATIVE_CONTIGS_DIR=os.path.dirname(os.path.abspath(VIRAL_CONTIGS)).rstrip("/")
	VIRAL_CONTIGS_BASE=""
	VIRAL_CONTIGS_DIR=""
	if RESULTS_DIR== "":
		RESULTS_DIR=REPRESENTATIVE_CONTIGS_DIR

print(VIRAL_CONTIGS_DIR)



RULES_DIR = 'rules'

CONFIDENCE_TYPES=["high", "low"]
SAMPLING_TYPE=config['sampling'].split()
SAMPLES=""

SRA_list=config['sra_list'].split()

CONTAMINANTS=config['contaminants_list'].split()
CONTAMINANTS.append("GCF_000819615.1")


NANOPORE=False
TOMBO=False
PAIRED=False
READ_TYPES=[config['forward_tag']]
POOLED=config['nanopore_pooled']
NANOPORE_SAMPLES=""
if not RAW_DATA_DIR == "":
	RAW_DATA_DIR=RAW_DATA_DIR.rstrip("/")
	SAMPLES,=glob_wildcards(RAW_DATA_DIR + "/{sample}_" + str(config['forward_tag']) + ".fastq")
	NANOPORE_SAMPLES,=glob_wildcards(RAW_DATA_DIR + "/{sample}_" + str(config['nanopore_tag']) + ".fastq")

	for fname in os.listdir(RAW_DATA_DIR):
		if fname.endswith(str(config['reverse_tag']) + '.fastq'):
			PAIRED=True
		elif fname.endswith(str(config['nanopore_tag']) + '.fastq'):
			NANOPORE=True
		elif fname.endswith(str(config['nanopore_tag']) + '_fast5_single'):
			TOMBO=True
else:
	RAW_DATA_DIR=RESULTS_DIR+"/00_RAW_DATA"


#NANOPORE_SAMPLES=SAMPLES

if PAIRED:
	READ_TYPES.append(config['reverse_tag'])
if POOLED:
	print("Nanopore reads are from a pooled sample")
	NANOPORE_SAMPLES=config['nanopore_pooled_name']
if len(SAMPLES)==1:
	SAMPLING_TYPE=["tot"]
SAMPLING_TYPE_TOT=["tot"]

dir_list = ["RULES_DIR","ENVS_DIR", "ADAPTERS_DIR", "CONTAMINANTS_DIR","RAW_DATA_DIR", "QC_DIR", "CLEAN_DATA_DIR", "ASSEMBLY_DIR", "VIRAL_DIR", "vOUT_DIR", "MAPPING_DIR", "MMSEQS", "ANNOTATION", "ASSEMBLY_TEST", "BENCHMARKS"]
dir_names = ["rules", "../envs", "db/adapters",  RESULTS_DIR + "/db/contaminants" ,RAW_DATA_DIR, RESULTS_DIR + "/01_QC", RESULTS_DIR + "/02_CLEAN_DATA", RESULTS_DIR + "/03_CONTIGS", VIRAL_CONTIGS_DIR , REPRESENTATIVE_CONTIGS_DIR ,RESULTS_DIR + "/06_MAPPING", RESULTS_DIR + "/08_MMSEQS", RESULTS_DIR + "/07_ANNOTATION", RESULTS_DIR + "/08_ASSEMBLY_TEST", RESULTS_DIR + "/BENCHMARK"]
dirs_dict = dict(zip(dir_list, dir_names))

print("Read Types = " )
print(*READ_TYPES, sep = ", ")

print("Sample Names = ")
print(*SAMPLES, sep = ", ")

print("Contaminants = ")
print(*CONTAMINANTS, sep = ", ")

print("Reference contigs = ")
print(REPRESENTATIVE_CONTIGS_BASE)
print(REPRESENTATIVE_CONTIGS_DIR)

print("Results Dir = ")
print(RESULTS_DIR)
print("Nanopore = ")
print(NANOPORE)

print("FAST5 = ")
print(TOMBO)

print("Nanopore samples= ")
print(NANOPORE_SAMPLES)
#======================================================
# Rules
#======================================================


rule all:
	input:
		expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample_nanopore}_contigs_flye.fasta",sample_nanopore=SAMPLES )
rule qualityCheckNanopore:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"]+"/{sample_nanopore}_nanopore.fastq",
	output:
		nanoqc_dir=temp(directory(dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanoplot")),
		nanoqc=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_preQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/qualityCheckNanopore/{sample_nanopore}.tsv"
#	threads: 1
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.raw_fastq}
		mv {output.nanoqc_dir}/nanoQC.html {output.nanoqc}
		"""

rule remove_adapters_quality_nanopore:
	input:
		raw_data=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore.fastq",
	output:
		trimmed_data=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq",
		porechopped=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_porechopped.fastq"),
	params:
		headcrop=50,
		tailcrop=50,
		quality=10,
	message:
		"Trimming Nanopore Adapters with Porechop"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/remove_adapters_quality_nanopore/{sample_nanopore}.tsv"
	threads: 2
	shell:
		"""
		porechop -i {input.raw_data} -o {output.porechopped} --threads {threads}
		NanoFilt -q {params.quality} -l 1000 --headcrop {params.headcrop} --tailcrop {params.tailcrop} {output.porechopped} > {output.trimmed_data}
		"""

rule remove_contaminants_nanopore:
	input:
		trimmed_data=dirs_dict["RAW_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq",
		contaminants_fasta=expand(dirs_dict["CONTAMINANTS_DIR"] +"/{contaminants}.fasta",contaminants=CONTAMINANTS),
	output:
		fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.fastq"),
		size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.txt",
		phix_contaminants_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{sample_nanopore}_nanopore_contaminants.fasta",
	message:
		"Remove contamination with Porechop"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/remove_contaminants_nanopore/{sample_nanopore}.tsv"
	threads: 2
	shell:
		"""
		cat {input.contaminants_fasta} > {output.phix_contaminants_fasta}
		minimap2 -ax map-ont {output.phix_contaminants_fasta} {input.trimmed_data} | samtools fastq -n -f 4 - > {output.fastq}
		grep -c "^@" {output.fastq} > {output.size}
		"""


rule postQualityCheckNanopore:
	input:
		fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.fastq"),
	output:
		nanoqc_dir=temp(directory(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanoqc_post")),
		nanoqc=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_postQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/postQualityCheckNanopore/{sample_nanopore}.tsv"
#	threads: 1
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.fastq}
		mv {output.nanoqc_dir}/nanoQC.html {output.nanoqc}
		"""

rule qualityStatsNanopore:
	input:
		fastq=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}.fastq"),
	output:
		nanostats=dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanostats_postQC.html",
	message:
		"Performing nanoQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/env3.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/qualityStatsNanopore/{sample_nanopore}.tsv"
#	threads: 1
	shell:
		"""
		NanoStat --fastq {input.fastq} > {output.nanostats}
		"""


rule asemblyFlye:
	input:
		nanopore=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.fastq",
	output:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample}/assembly.fasta",
		scaffolds_final=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_flye.fasta"
	message:
		"Assembling Nanopore reads with Flye"
	params:
		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample}",
		genome_size="20m"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/asemblyFlye/{sample}.tsv"
	threads: 4
	shell:
		"""
		flye --nano-raw {input.nanopore} --out-dir {params.assembly_dir} --genome-size {params.genome_size} --meta --threads {threads}
		cp {output.scaffolds} {output.scaffolds_final}
		""rule
