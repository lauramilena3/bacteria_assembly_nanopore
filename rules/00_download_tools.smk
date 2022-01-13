rule get_SRAToolkit:
	output:
		SRAToolkit_dir=directory("tools/sratoolkit.2.10.0-ubuntu64"),
	message:
		"Downloading VIBRANT"
	params:
		tools="tools",
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml",
	threads:
		16
	shell:
		"""
		cd {params.tools}
		wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz
		tar -xzf sratoolkit.2.10.0-ubuntu64.tar.gz
		"""
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

rule get_VIBRANT:
	output:
		VIBRANT_dir=directory(os.path.join(workflow.basedir, config['vibrant_dir'])),
	message:
		"Downloading VIBRANT"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/AnantharamanLab/VIBRANT
		chmod -R 744 VIBRANT
		cd VIBRANT/databases
		./VIBRANT_setup.py
		"""
rule getQUAST:
	output:
		quast_dir=directory(config["quast_dir"])
	message:
		"Downloading QUAST"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 4
	shell:
		"""
		curl -OL https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
		tar -xzf quast-5.0.2.tar.gz -C tools
		cd {config[quast_dir]}
		./setup.py install
		"""
rule get_mmseqs:
	output:
		mmseqs_dir=directory(os.path.join(workflow.basedir, config['mmseqs_dir'])),
		refseq=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna")),
		refseq_taxid=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna.taxidmapping")),
	message:
		"Downloading MMseqs2"
	params:
		taxdump=(os.path.join(workflow.basedir,"db/ncbi-taxdump/")),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 8
	shell:
		"""
		MM_dir={output.mmseqs_dir}
		echo $MM_dir
		if [ ! -d $MM_dir ]
		then
			mkdir -p tools
			cd tools
			git clone https://github.com/soedinglab/MMseqs2.git
			cd MMseqs2
			mkdir build
			cd build
			cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
			make -j {threads}
			make install
		fi
		#download taxdump
		cd ../../../db
		mkdir -p ncbi-taxdump
		cd ncbi-taxdump
		wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		tar -xzvf taxdump.tar.gz
		#download RefSeqViral
		wget ftp://ftp.ncbi.nlm.nih.gov/blast//db/ref_viruses_rep_genomes.tar.gz
		tar xvzf ref_viruses_rep_genomes.tar.gz
		blastdbcmd -db ref_viruses_rep_genomes -entry all > {output.refseq}
		blastdbcmd -db ref_viruses_rep_genomes -entry all -outfmt "%a %T" > {output.refseq_taxid}
		{output.mmseqs_dir}/build/bin/mmseqs createdb {output.refseq} RefSeqViral.fnaDB
		{output.mmseqs_dir}/build/bin/mmseqs createtaxdb RefSeqViral.fnaDB tmp --ncbi-tax-dump {params.taxdump} --tax-mapping-file {output.refseq_taxid}
		"""
rule get_ALE:
	output:
		ALE_dir=directory(config['ALE_dir']),
	message:
		"Downloading ALE"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/sc932/ALE.git
		cd ALE/src
		make
		"""
rule get_VIGA:
	output:
		VIGA_dir=directory(os.path.join(workflow.basedir, config['viga_dir'])),
		piler_dir=directory(os.path.join(workflow.basedir, config['piler_dir'])),
		trf_dir=directory(os.path.join(workflow.basedir, config['trf_dir'])),
	message:
		"Downloading MMseqs2"
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone --depth 1 https://github.com/lauramilena3/viga.git
		chmod 744 viga/create_dbs.sh viga/VIGA.py
		cd viga
		./create_dbs.sh
		cd ..
		wget https://www.drive5.com/pilercr/pilercr1.06.tar.gz --no-check-certificate
		tar -xzvf pilercr1.06.tar.gz
		cd pilercr1.06
		make
		cd ..
		mkdir TRF
		cd TRF
		wget http://tandem.bu.edu/trf/downloads/trf409.linux64
		wget http://tandem.bu.edu/irf/downloads/irf307.linux.exe
		mv trf409.linux64 trf
		mv irf307.linux.exe irf
		chmod 744 trf irf
		"""

# rule downloadViralTools:
# 	output:
# 		virSorter_dir=directory(config['virSorter_dir']),
# 		virFinder_dir=directory(config['virFinder_dir']),
# 	message:
# 		"Downloading required VirSorter and VirFinder"
# 	threads: 1
# 	shell:
# 		"""
# 		#VIRSORTER
# 		VS_dir="{config[virSorter_dir]}"
# 		echo $VS_dir
# 		if [ ! -d $VS_dir ]
# 		then
# 			mkdir -p tools
# 			cd tools
# 			git clone https://github.com/simroux/VirSorter.git
# 			cd VirSorter/Scripts
# 			make clean
# 			make
# 			cd ../../../
# 		fi
# 		#VIRFNDER
# 		VF_dir="{config[virFinder_dir]}"
# 		echo $VF_dir
#			if [ ! -d $VF_dir ]
# 		then
# 			if [ ! {config[operating_system]} == "linux" ]
# 			then
# 				curl -OL https://raw.github.com/jessieren/VirFinder/blob/master/mac/VirFinder_1.1.tar.gz?raw=true
# 			else
# 				curl -OL https://github.com/jessieren/VirFinder/blob/master/linux/VirFinder_1.1.tar.gz?raw=true
# 			fi
# 			mkdir -p {output.virFinder_dir}
# 			mv VirFinder*tar.gz* {output.virFinder_dir}/VirFinder_1.1.tar.gz
# 		fi
# 		"""
rule downloadVirSorterDB:
	output:
		virSorter_dir=directory(config['virSorter_dir']),
	message:
		"Downloading VirSorter database"
	threads: 8
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	params:
		virSorter_db="db/VirSorter"
	shell:
		"""
		#git clone https://github.com/jiarong/VirSorter2.git {output.virSorter_dir}
		#cd {output.virSorter_dir}
		#pip install .
		virsorter setup -d db -j {threads}
		mkdir {output.virSorter_dir}
		"""
rule downloadCheckvDB:
	output:
		checkv_db=directory(config['checkv_db']),
	message:
		"Downloading CheckV database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	shell:
		"""
		checkv download_database ./db
		"""

rule getKrakenTools:
	output:
		kraken_tools=directory(config['kraken_tools']),
	message:
		"Downloading KrakenTools"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/jenniferlu717/KrakenTools
		chmod 777 KrakenTools/*
		"""

rule downloadKrakenDB:
	output:
		kraken_db=directory(config['kraken_db']),
	message:
		"Downloading miniKraken database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	shell:
		"""
		wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
		mkdir {output.kraken_db}
		tar -xvf minikraken_8GB_202003.tgz -C {output.kraken_db}
		"""

rule downloadKrakenDB_human:
	output:
		kraken_db_human=directory(config['kraken_db_human']),
	message:
		"Downloading human Kraken database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	shell:
		"""
		kraken2-build --download-library human --db {output.kraken_db_human} --threads {threads} --use-ftp
		kraken2-build --download-taxonomy --db {output.kraken_db_human}
		kraken2-build --build --db {output.kraken_db_human} --threads {threads}
		kraken2-build --clean --db {output.kraken_db_human}
		"""

# rule downloadVirSorterDB:
# 	output:
# 		virSorter_db=directory(config['virSorter_db']),
# 	message:
# 		"Downloading VirSorter database"
# 	threads: 1
# 	params:
# 		virSorter_db="db/VirSorter"
# 	shell:
# 		"""
# 		VS_db="{config[virSorter_db]}"
# 		echo $VS_db
# 		if [ ! -d $VS_db ]
# 		then
# 			curl -OL https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
# 			mkdir -p {params.virSorter_db}
# 			tar -xvzf virsorter-data-v2.tar.gz -C {params.virSorter_db}
# 			rm virsorter-data-v2.tar.gz
# 		fi
# 		"""
rule downloadVcontact2Files:
	output:
		gene2genome=(os.path.join(workflow.basedir,"db/vcontact2/gene-to-genome.30May2020.csv")),
		vcontact_format=(os.path.join(workflow.basedir,"db/vcontact2/vcontact_format_30May2020.faa")),
	message:
		"Downloading vContact2 formatting database"
	threads: 1
	params:
	shell:
		"""
		wget http://s3.climb.ac.uk/ADM_share/crap/website/gene-to-genome.30May2020.csv.gz
		gunzip -c gene-to-genome.30May2020.csv.gz > {output.gene2genome}
		wget http://s3.climb.ac.uk/ADM_share/crap/website/vcontact_format_30May2020.faa.gz
		gunzip -c vcontact_format_30May2020.faa.gz > {output.vcontact_format}
		dos2unix {output.gene2genome}
		dos2unix {output.vcontact_format}
		#sed -i "s/None/None_provided/g" {output.gene2genome}
		"""

rule downloadBLASTviralProteins:
	output:
		blast=(os.path.join(workflow.basedir,"db/ncbi/NCBI_viral_proteins.faa")),
	message:
		"Downloading RefSeq viral proteins for blast annotation"
	threads: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml",
	shell:
		"""
		esearch -db "protein" -query "txid10239[Organism:exp] AND (viruses[filter] AND refseq[filter])" \
			| efetch -format fasta > {output.blast}
		makeblastdb -in {output.blast} -dbtype prot
		"""
rule getClusterONE:
	output:
		clusterONE_dir=directory(config["clusterONE_dir"]),
	message:
		"Downloading clusterONE"
	threads: 1
	shell:
		"""
		mkdir -p {output.clusterONE_dir}
		curl -OL  http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar
		mv cluster_one-1.0.jar {output.clusterONE_dir}
		chmod 744 {output.clusterONE_dir}/cluster_one-1.0.jar
		"""
rule downloadCanu:
	output:
		canu_dir=directory(config['canu_dir']),
	message:
		"Installing Canu assembler"
	threads: 1
	shell:
		"""
		if [ ! -d {output.canu_dir} ]
		then
			if [ {config[operating_system]} == "macOs" ]
			then
				mkdir -p tools
				curl -OL https://github.com/marbl/canu/releases/download/v2.0/canu-2.0.Darwin-amd64.tar.xz
			else
				mkdir -p tools
				curl -OL https://github.com/marbl/canu/releases/download/v2.0/canu-2.0.Linux-amd64.tar.xz
			fi
		fi
		tar -xJf canu-2.0.*.tar.xz -C tools
		"""
rule get_WIsH:
	input:
		representative_list="db/PATRIC/representatives_referece_bacteria_archaea_acc.txt",
	output:
		wish_dir=directory(os.path.join(workflow.basedir, (config['wish_dir']))),
		FNA=directory("db/PATRIC/FNA"),
	params:
		patric_dir="db/PATRIC"
	message:
		"Downloading WIsH"
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 16
	shell:
		"""
		wish_dir={output.wish_dir}
		echo $wish_dir
		if [ ! -d $wish_dir ]
		then
			mkdir -p tools
			cd tools
			git clone https://github.com/soedinglab/WIsH.git
			cd WIsH
			cmake .
			make -j {threads}
		fi
		cd ../..
		mkdir {output.FNA}
		cd {output.FNA}
		#COUNTER=1
		#for i in $(cat < ../.{input.representative_list}); do acc="${{i%.*}}"; echo $acc; COUNTER=$[COUNTER + 1]; echo $COUNTER; wget -qN "ftp://ftp.patricbrc.org/genomes/$i/$i.fna" & done;
		cat ../../../{input.representative_list} | xargs -I {{}} -n 1 -t -P {threads} wget -qN "ftp://ftp.patricbrc.org/genomes/{{}}/{{}}.fna"
		"""

rule get_WTP:
	output:
		WTP_dir=directory(os.path.join(workflow.basedir, config['WTP_dir'])),
	message:
		"Downloading What the Phage"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	threads: 1
	shell:
		"""c
		mkdir -p tools
		cd tools
		mkdir {output.WTP_dir}
		cd {output.WTP_dir}
		singularity pull  --name nanozoo-sourmash-3.4.1--16a8db7.img docker://nanozoo/sourmash:3.4.1--16a8db7
		nextflow run replikation/What_the_Phage -r v1.0.0 --setup -profile local,singularity --cachedir cache_dir
		"""
