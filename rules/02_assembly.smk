rule asemblyFlye:
	input:
		fastq=dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_nanofilt.fastq",
	output:
		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample_nanopore}/assembly.fasta",
		scaffolds_final=dirs_dict["ASSEMBLY_DIR"] + "/{sample_nanopore}_contigs_flye.fasta"
	message:
		"Assembling Nanopore reads with Flye"
	params:
		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample_nanopore}",
		genome_size="20m"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/asemblyFlye/{sample_nanopore}.tsv"
	threads: 4
	shell:
		"""
		flye --nano-raw {input.fastq} --out-dir {params.assembly_dir} --genome-size {params.genome_size} --meta --threads {threads}
		cp {output.scaffolds} {output.scaffolds_final}
		"""

rule shortReadAsemblySpadesPE:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_norm.fastq"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_norm.fastq"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.fastq",
	output:
		scaffolds=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.fasta"),
		filtered_list=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades.filtered_list.txt"),
		assembly_graph=dirs_dict["ASSEMBLY_DIR"] +"/{sample}_assembly_graph_spades.fastg",
	params:
		raw_scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades.scaffolds.fasta",
		assembly_graph=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades.assembly_graph.fastg",
		assembly_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades"),
	message:
		"Assembling PE reads with metaSpades"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/shortReadAsemblySpadesPE/{sample}.tsv"
	threads: 8
	shell:
		"""
		spades.py  --pe1-1 {input.forward_paired} --pe1-2 {input.reverse_paired}  --pe1-s {input.unpaired} -o {params.assembly_dir} \
		-t {threads} --only-assembler --memory 350
		grep "^>" {params.raw_scaffolds} | sed s"/_/ /"g | awk '{{ if ($4 >= {config[min_len]} && $6 >= {config[min_cov]}) print $0 }}' \
		| sort -k 4 -n | sed s"/ /_/"g | sed 's/>//' > {output.filtered_list}
		seqtk subseq {params.raw_scaffolds} {output.filtered_list} > {output.scaffolds}
		cp {params.assembly_graph} {output.assembly_graph}
		sed "s/>/>{wildcards.sample}_/g" -i {output.scaffolds}
		"""

rule assemblyStatsILLUMINA:
	input:
		quast_dir=(config["quast_dir"]),
		scaffolds_spades=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.fasta", sample=SAMPLES)
	output:
		quast_report_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast"),
		quast_txt=dirs_dict["ASSEMBLY_DIR"] + "/assembly_quast_report.txt"
	message:
		"Creating assembly stats with quast"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/assemblyStatsILLUMINA.tsv"
	threads: 4
	shell:
		"""
		{input.quast_dir}/quast.py {input.scaffolds_spades} -o {output.quast_report_dir} --threads {threads}
		cp {output.quast_report_dir}/report.txt {output.quast_txt}
		"""
#
# rule busco_assesment:
# 	input:
# 		scaffolds_spades=expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.fasta", sample=SAMPLES)
# 	output:
# 		quast_report_dir=directory(dirs_dict["ASSEMBLY_DIR"] + "/statistics_quast"),
# 		quast_txt=dirs_dict["ASSEMBLY_DIR"] + "/assembly_quast_report.txt"
# 	message:
# 		"Creating assembly stats with quast"
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/assemblyStatsILLUMINA.tsv"
# 	threads: 4
# 	shell:
# 		"""
# 		{input.quast_dir}/quast.py {input.scaffolds_spades} -o {output.quast_report_dir} --threads {threads}
# 		cp {output.quast_report_dir}/report.txt {output.quast_txt}
# 		"""
