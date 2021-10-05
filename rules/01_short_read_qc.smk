rule downloadContaminants:
	output:
		contaminant_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{contaminant}.fasta",
		contaminant_dir=temp(directory(dirs_dict["CONTAMINANTS_DIR"] +"/temp_{contaminant}")),
	message:
		"Downloading contaminant genomes"
	params:
		contaminants_dir=dirs_dict["CONTAMINANTS_DIR"],
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml",
	threads:
		16
	shell:
		"""
		mkdir {output.contaminant_dir}
		cd {output.contaminant_dir}
		wget $(esearch -db "assembly" -query {wildcards.contaminant} | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{{print $0"/"$NF"_genomic.fna.gz"}}')
		gunzip -f *gz
		cat *fna >> {output.contaminant_fasta}
		"""

rule qualityCheckIllumina:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"]+"/{sample}_{reads}.fastq"
	output:
		html=(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{reads}_fastqc.html"),
		zipped=temp(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{reads}_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/qualityCheckIllumina/{sample}_{reads}.tsv"
#	threads: 1
	shell:
		"""
		fastqc {input}
		"""

rule multiQC:
	input:
		html=expand(dirs_dict["RAW_DATA_DIR"]+"/{sample}_{reads}_fastqc.html", sample=SAMPLES, reads=READ_TYPES),
		zipped=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{reads}_fastqc.zip", sample=SAMPLES, reads=READ_TYPES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
	params:
		fastqc_dir=dirs_dict["RAW_DATA_DIR"],
		html_name="preQC_illumina_report.html",
		multiqc_dir=dirs_dict["QC_DIR"],
	message:
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/multiQC/multiqc.tsv"
	shell:
		"""
		multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule trim_adapters_quality_illumina_PE:
	input:
		forward=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq",
		reverse=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq",
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
		forward_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq"),
		reverse_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired.fastq"),
		trimmomatic_values=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_trimmomatic_values.txt"),
	params:
		adapters=dirs_dict["ADAPTERS_DIR"] + "/" + config['adapters_file']
	message:
		"Trimming Illumina Adapters with Trimmomatic"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/trim_adapters_quality_illumina_PE/{sample}.tsv"
	threads: 8
	shell:
		"""
		echo leading {config[trimmomatic_leading]} trailing {config[trimmomatic_trailing]} winsize {config[trimmomatic_window_size]} winqual {config[trimmomatic_window_quality]} minlnth {config[trimmomatic_minlen]} > {output.trimmomatic_values}
		trimmomatic PE -threads {threads} -phred33 {input.forward} {input.reverse} \
		{output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
		ILLUMINACLIP:{params.adapters}:2:30:10:2:true LEADING:{config[trimmomatic_leading]} TRAILING:{config[trimmomatic_trailing]} \
		SLIDINGWINDOW:{config[trimmomatic_window_size]}:{config[trimmomatic_window_quality]} MINLEN:{config[trimmomatic_minlen]}
		cat {output.forward_unpaired} {output.reverse_unpaired} > {output.unpaired}
		"""

# rule remove_contaminants_PE:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
# 		forward_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq",
# 		reverse_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq",
# 		contaminants_fasta=expand(dirs_dict["CONTAMINANTS_DIR"] +"/{contaminants}.fasta",contaminants=CONTAMINANTS),
# 	output:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_bbduk.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_bbduk.fastq"),
# 		merged_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.fastq"),
# 		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired_bbduk.fastq"),
# 		phix_contaminants_fasta=dirs_dict["CONTAMINANTS_DIR"] +"/{sample}_contaminants.fasta"
# 	message:
# 		"Removing phiX174 and user given contaminants with BBtools"
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/remove_contaminants_PE/{sample}.tsv"
# 	threads: 4
# 	resources:
# 		mem_gb=40
# 	shell:
# 		"""
# 		cat {input.contaminants_fasta} > {output.phix_contaminants_fasta}
# 		#PE
# 		#PAIRED
# 		bbduk.sh -Xmx{resources.mem_gb}g in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} ref={output.phix_contaminants_fasta} k=31 hdist=1 threads={threads}
# 		#UNPAIRED
# 		cat {input.forward_unpaired} {input.reverse_unpaired} > {output.merged_unpaired}
# 		bbduk.sh -Xmx{resources.mem_gb}g in={output.merged_unpaired} out={output.unpaired} ref={output.phix_contaminants_fasta} k=31 hdist=1 threads={threads}
# 		"""

# rule remove_human_PE:
# 	input:
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_bbduk.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_bbduk.fastq"),
# 		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired_bbduk.fastq"),
# 		kraken_db_human=(config['kraken_db_human']),
# 	output:
# 		forward_paired_temp=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R_1.fastq"),
# 		reverse_paired_temp=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R_2.fastq"),
# 		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.fastq"),
# 		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.fastq"),
# 		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.fastq",
# 		paired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_paired_clean.txt",
# 		unpaired_size=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.txt",
# 		kraken_output_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-kraken2-out_paired.txt",
# 		kraken_report_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-kraken2-report_paired.txt",
# 		kraken_output_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-kraken2-out_unpaired.txt",
# 		kraken_report_unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}-kraken2-report_unpaired.txt",
# 	message:
# 		"Removing human reads with Kraken"
# 	params:
# 		unclassified_name_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R#.fastq",
# 	conda:
# 		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
# 	threads: 4
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/remove_human_PE/{sample}.tsv"
# 	resources:
# 		mem_gb=40
# 	shell:
# 		"""
# 		#PAIRED
# 		kraken2 --paired --db {input.kraken_db_human} --threads {threads} --output {output.kraken_output_paired} \
# 				--report {output.kraken_report_paired} --unclassified-out {params.unclassified_name_paired} \
# 				{input.forward_paired} {input.reverse_paired}
# 		cp {output.forward_paired_temp} {output.forward_paired}
# 		cp {output.reverse_paired_temp} {output.reverse_paired}
# 		grep -c "^@" {output.forward_paired} > {output.paired_size}
# 		#UNPAIRED
# 		kraken2 --db {input.kraken_db_human} --threads {threads} --output {output.kraken_output_unpaired} \
# 				--report {output.kraken_report_unpaired} --unclassified-out {output.unpaired} {input.unpaired}
# 		grep -c "^@" {output.unpaired} > {output.unpaired_size} ||  echo "0" > {output.unpaired_size}
# 		"""


rule postQualityCheckIlluminaPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired.fastq"),
	output:
		html_forward=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_forward_paired_fastqc.html"),
		zipped_forward=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_forward_paired_fastqc.zip"),
		html_reverse=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_fastqc.html"),
		zipped_reverse=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_fastqc.zip"),
		html_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_unpaired_fastqc.html"),
		zipped_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_unpaired_fastqc.zip"),
	params:
		postQC_dir=dirs_dict["CLEAN_DATA_DIR"] +"/postQC",
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/postQualityCheckIlluminaPE/{sample}.tsv"
#	threads: 1
	shell:
		"""
		fastqc {input.forward_paired} -o {params.postQC_dir}
		fastqc {input.reverse_paired} -o {params.postQC_dir}
		fastqc {input.unpaired} -o {params.postQC_dir}
		"""

rule postMultiQC:
	input:
		html_forward=expand(dirs_dict["CLEAN_DATA_DIR"]  + "/postQC" + "/{sample}_forward_paired_fastqc.html", sample=SAMPLES),
		zipped_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_forward_paired_fastqc.zip", sample=SAMPLES),
		html_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_fastqc.html", sample=SAMPLES),
		zipped_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_reverse_paired_fastqc.zip", sample=SAMPLES),
		html_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/postQC" + "/{sample}_unpaired_fastqc.html", sample=SAMPLES),
		zipped_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"]  + "/postQC" + "/{sample}_unpaired_fastqc.zip", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/postQC_illumina_report.html"
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"] +  "/postQC",
		html_name="postQC_illumina_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/postMultiQC/multiqc.tsv"
	shell:
		"""
		multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule normalizeReads_PE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
		unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired.fastq"),
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.fastq",
	message:
		"Normalizing reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/normalizeReads_PE/{sample}.tsv"
	params:
		min_depth=config['min_norm'],
		max_depth=config['max_norm']
	threads: 4
	resources:
		mem_mb=6000
	shell:
		"""
		#PE
		#paired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in1={input.forward_paired} in2={input.reverse_paired} out1={output.forward_paired} out2={output.reverse_paired} \
		target={params.max_depth} mindepth={params.min_depth} t={threads}
		#unpaired
		bbnorm.sh -Xmx{resources.mem_mb}m ecc in={input.unpaired} out={output.unpaired} target={params.max_depth} mindepth={params.min_depth}
		"""

rule kmer_rarefraction:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
	output:
		histogram=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram.csv"),
	message:
		"Counting unique reads with BBtools"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kmer_rarefraction/{sample}.tsv"
	threads: 1
	shell:
		"""
		bbcountunique.sh in1={input.forward_paired} in2={input.reverse_paired} out={output.histogram} interval={config[kmer_window]}
		"""

rule plot_kmer:
	input:
		histograms=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kmer_histogram...csv", sample=SAMPLES),
	output:
		plot=(dirs_dict["CLEAN_DATA_DIR"] + "/kmer_rarefraction_plot.png"),
		svg=(dirs_dict["CLEAN_DATA_DIR"] + "/kmer_rarefraction_plot.svg"),

	message:
		"Plot unique reads with BBtools"
	threads: 1
	run:
		import pandas as pd
		import seaborn as sns; sns.set()
		import matplotlib.pyplot as plt

		plt.figure(figsize=(12,12))
		sns.set(font_scale=2)
		sns.set_style("whitegrid")

		read_max=0

		for h in input.histograms:
			df=pd.read_csv(h, sep="\t")
			df.columns=["count", "percent", "c", "d", "e", "f", "g", "h", "i", "j"]
			df=df[["count", "percent"]]
			ax = sns.lineplot(x="count", y="percent", data=df,err_style='band', label=h.split("/")[-1].split("_kmer")[0])
			read_max=max(read_max,df["count"].max())

		ax.set(ylim=(0, 100))
		ax.set(xlim=(0, read_max*1.2))

		ax.set_xlabel("Read count",fontsize=20)
		ax.set_ylabel("New k-mers (%)",fontsize=20)
		ax.figure.savefig(output.plot)
		ax.figure.savefig(output.svg, format="svg")

rule contaminants_KRAKEN:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq"),
		kraken_db=(config['kraken_db']),
	output:
		kraken_output=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output.csv"),
		kraken_report=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report.csv"),
		kraken_domain=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_domains.csv"),
	params:
		kraken_db=config['kraken_db']+ "/minikraken2_v2_8GB_201904_UPDATE",
	message:
		"Kraken contamination"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kmer_rarefraction/{sample}.tsv"
	threads: 4
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output} --report {output.kraken_report}

		grep -P 'D\t' {output.kraken_report} | sort -r > {output.kraken_domain}
		"""
