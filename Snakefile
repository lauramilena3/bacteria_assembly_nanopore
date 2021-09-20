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

if RESULTS_DIR == "" and not RAW_DATA_DIR == "":
	RESULTS_DIR=os.path.abspath(os.path.join(RAW_DATA_DIR, os.pardir))



RULES_DIR = 'rules'

SAMPLING_TYPE=config['sampling'].split()
SAMPLES=""


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



dir_list = ["RULES_DIR","ENVS_DIR", "ADAPTERS_DIR", "CONTAMINANTS_DIR","RAW_DATA_DIR", "QC_DIR", "CLEAN_DATA_DIR", "ASSEMBLY_DIR",  "MAPPING_DIR", "MMSEQS", "ANNOTATION", "ASSEMBLY_TEST", "BENCHMARKS"]
dir_names = ["rules", "envs", "db/adapters",  RESULTS_DIR + "/db/contaminants" ,RAW_DATA_DIR, RESULTS_DIR + "/01_QC", RESULTS_DIR + "/02_CLEAN_DATA", RESULTS_DIR + "/03_CONTIGS" ,RESULTS_DIR + "/06_MAPPING", RESULTS_DIR + "/08_MMSEQS", RESULTS_DIR + "/07_ANNOTATION", RESULTS_DIR + "/08_ASSEMBLY_TEST", RESULTS_DIR + "/BENCHMARK"]
dirs_dict = dict(zip(dir_list, dir_names))

print("Read Types = " )
print(*READ_TYPES, sep = ", ")

print("Sample Names = ")
print(*SAMPLES, sep = ", ")

print("Contaminants = ")
print(*CONTAMINANTS, sep = ", ")

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


def inputAll(wildcards):
	inputs=[]
	inputs.append(dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html")
	inputs.append(dirs_dict["QC_DIR"]+ "/postQC_illumina_report.html")
	inputs.append(dirs_dict["ASSEMBLY_DIR"] + "/assembly_quast_report.txt")

	if NANOPORE:
		inputs.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_preQC.html",sample_nanopore=NANOPORE_SAMPLES))
		inputs.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_postQC.html", sample_nanopore=NANOPORE_SAMPLES))
		inputs.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_postQC.html", sample_nanopore=NANOPORE_SAMPLES))

	return inputs


rule all:
	input:
		inputAll,



include: os.path.join(RULES_DIR, '01_long_read_qc.smk')
include: os.path.join(RULES_DIR, '01_short_read_qc.smk')
include: os.path.join(RULES_DIR, '02_assembly.smk')
