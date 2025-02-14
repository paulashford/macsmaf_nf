id_conversion - notes

Manual convert pre_processed modules for K1
	output/pre_processed_networks/network_modules_K1_humanbase.dat (entrez_id -> ENSG)
	output/pre_processed_networks/network_modules_K1_string.dat (ESP -> ENSG)
reason? explained here: workflows/workflow_README.txt


cut [2] from  
	/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/id_conversion/main.nf
// nf_mod_files            =   "/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/pascal_gwas/output/d024_network_modules_K1_cpdb.dat"
    // nf_mod_files            =   "/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/pascal_gwas/output/d024_network_modules_K1_humanbase.dat"
    // nf_mod_files            =   "/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/pascal_gwas/output/d024_network_modules_K1_string.dat"
    // For HUGO conv...
    // nf_mod_files    =   "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e023/pascal/cpdb/network_modules_K1_cpdb_ENSG.dat"
    // nf_mod_files    =   "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e023/pascal/humanbase/network_modules_K1_humanbase_ENSG.dat"
    // nf_mod_files        =   "/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e023/pascal/humanbase/network_modules_K1_humanbase_ENSG_del_empty_modules.dat"


cut [1] from  
	/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/id_conversion/main.nf
// subworkflow test
// cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/subworkflows/id_convert
// export NF_CONFIG=/Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev/conf/base.config
// nextflow run main.nf -c "${NF_CONFIG}"

// Hi Ash, I looked at the file I got from you before and it was the same format for the GWAS data. 
// The only thing was that it contained all the modules from all methods. I'm attaching it below for your reference. Could you please generate the same file but with all three methods and networks (and the same naming system)?

// Name    chi2Pvalue  empPvalue
// ENSG-R1-humanbase-coval-0.3.dat_158 3.34001336E-4   3.42E-4
// ENSG-R1-humanbase-coval-0.5.dat_122 3.34001336E-4   3.36E-4
// ENSG-M1-cpdb-coval-0.8.dat_202  5.1915018E-4    1.88E-3
// ENSG-R1-humanbase-coval-0.4.dat_239 5.73738664E-4   1.72E-4
// ENSG-R1-cpdb-coval-0.5.dat_49   8.42244723E-4   2.63E-3
// ENSG-M1-cpdb-coval-0.8.dat_89   8.78298669E-4   1.158E-2
// ENSG-M1-cpdb-coval-0.1.dat_98   1.20799046E-3   8.31E-3
// ENSG-R1-humanbase-coval-0.3.dat_67  1.56715245E-3   1.12E-3
// ENSG-M1-cpdb-coval-0.7.dat_128  1.81988014E-3   4.72E-3

// combine netwokr module files in dir (eg for individual module files with different edge threshold cut-offs)
// create a 'labelled' net module file with filename in 1st col, so that file now contains metadata info from path name
// eg: grep -R "${dataset_dir}" --include="*.dat" -e '' > "${exper_dir}/k1_modules_grep.txt" 2>/dev/null
