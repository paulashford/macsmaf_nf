compiled_module_datasets from NextFlow outputs -
	manually versioned copies of sub-directories in git/macsmaf/macsr_nf/macsr_nf_dev/output
13 02 2025
p.ashford@ucl.ac.uk

NextFlow runs:
	git/macsmaf/macsr_nf/macsr_nf_dev/workflows/run_workflows.sh

git/macsmaf/macsr_nf/macsr_nf_dev/output
	/pre_processed_networks - generated with NextFlow / params.preproc_net_modules = true
		executes optional macs_nf_dev.nf/REPROC_NET_MODULES to gather .dat files from supplied path containing DREAM/Monet module files, 1 for each method, db and cut-off, e.g.
			nf_network_modules_dir  = "${params.user_root_dir}/Dropbox/bioinf/MACSMAF/datasets/d024" 
			
			GeneIDs are either original IDs (K1) or ENSG throughout (M1 and R1) 
				# original gene id type in source database (id_conversion/id_conversion_ad_hoc.nf; params.input_network_file_type = 'original')
				NET_DB_ORIGINAL_GENE_ID_TYPE_DICT = {'cpdb': 'ensembl_gene_id', 'humanbase': 'entrez_id', 'string': 'Protein_stable_ID'}
				# pre-processed gene id type (note: for K1 humanbase and string, gene ID conversion was done via id_conversion_ad_hoc.nf
				NET_DB_PRE_PROCESSED_GENE_ID_TYPE_DICT = {'cpdb': 'ensembl_gene_id', 'humanbase': 'ensembl_gene_id', 'string': 'ensembl_gene_id'}

		If params.preproc_net_modules = false, the above process is skipped and network_inputs channel is expected to contain pre_processed .dat files (1 for each method) in {params.nf_network_modules_dir}/network_processed
	
	/converted_networks - generated with NextFlow workflows/id_conversion_ad_hoc.nf
		ad-hoc id conversion of K1 humanbase and string networks from original Entrez IDs and ENSPs respectively to ENSG IDs, in order to match the M1 and R1 networks
		/batch_processed_networks/
			/hgnc - contains .dat files with geneIDs converted using HGNC map file
			/biomart - contains .dat files with geneIDs converted using biomart map file
			

	v01/	13 02 2025 pre-processed networks and ad-hoc id conversion of K1 humanbase and string networks to ENSG IDs to match the M1 and R1 networks
		pre_processed_networks/
			cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/output
			cp pre_processed_networks/*.dat compiled_module_datasets/v01/pre_processed_networks/
		converted_networks/ensembl_gene_id/
			cp pre_processed_networks/*M1_*.dat compiled_module_datasets/v01/converted_networks/ensembl_gene_id/
			cp pre_processed_networks/*R1_*.dat compiled_module_datasets/v01/converted_networks/ensembl_gene_id/
			cp pre_processed_networks/network_modules_K1_cpdb.dat compiled_module_datasets/v01/converted_networks/ensembl_gene_id/
			cp converted_networks/ad_hoc_id_conversion/network_modules_K1_humanbase_converted_original_to_ensembl_gene_id.dat compiled_module_datasets/v01/converted_networks/ensembl_gene_id/network_modules_K1_humanbase.dat
			cp converted_networks/ad_hoc_id_conversion/network_modules_K1_string_converted_original_to_Gene_stable_ID.dat compiled_module_datasets/v01/converted_networks/ensembl_gene_id/network_modules_K1_string.dat

			
			
