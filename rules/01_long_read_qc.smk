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
		trimmed_data=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq",
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


rule postQualityCheckNanopore:
	input:
		fastq=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq",
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
		fastq=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq",
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
