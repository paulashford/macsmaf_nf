e019 Original network clusters g:profiler and slims analysis 

NOTE 26 04 2024: copied from Dropbox/MACSMAF/experiments/e019

07 12 2022
NOTE: equivalent to e017, but that was for the CSR-specific subnetwork. This is for ALL modules in original files.

NOTE: ** Requires programmatic use of g:Profiler, as the full module queries are causing time-outs on web-page **

NOTE: e017 notes below - see relevant scripts (named below) for up-to-date!

1. CSR net (see ey005) clustered using DREAM/Monet methods R1 and M1 and K1 (see e017.sh).

2. g:Profiler enrichment [Functional enrichment of resulting modules against multiple databases (KEGG, GO, HP etc) ]
    a) Use Monet module files eg:
        experiments/e017/monet_results_outputs/2022-12-07-191927__K1_n200_result-modules__csr-edge-data_simple_no_HDR.tsv
    b) Use VSCode to amend format from:
            265	1	IFI30	RAB31	HCK	PECAM1	CTSS
            266	1	AIF1	LAPTM5	CCL5	GMFG	BCL2A1
            267	1	ARHGDIB	NCF2	NCF4	RAC2
        to g:Profiler 'multi-query' format eg
            >e017_R1_mod_7
            JAK1	IFNLR1	IL10	IL10RB	IL10RA	IFNL3	IFNL1	IL21R
            >e017_R1_mod_8
            JAK2	IFNG	IFNGR1	IFNGR2	LEPR	LEP	CDCP1	ST14
            >e017_R1_mod_9
    c) Save as:
        gprofiler_multi/mod_annot_R1.csv etc 
    d) run g:Profiler with settings as per experiments/e017/gprofiler_multi/gp_params.txt
        (Note have been usng options to save gprofiler run parameters with the output file as prob v useful in future!)

3. Network annotation - clustering by module
    This is simply done by pivoting the Monet module files (pre gprofiler) using:
        experiments/e017/process_modules.R  [note need for 'dummy' header on module files see R script] 
    save output in: 
        annot/mod_annot_R1.csv etc 

4. Functional scoring
    a) Pivot the g:Profiler output as per: 
    experiments/e017/proccess_gprofiler.R  [gives more info than (3) which just has gene/module number in the pivot] [note need for 'dummy' header on module files see R script] 

    b) gprofiler analysis - see:
        proccess_gprofiler.R 
        GOslims: Dropbox/bioinf/MACSMAF/datasets/d011
        which uses owltools: /Users/ash/git/bioinf/owltools
        x-scrivener-item:///Users/ash/Dropbox/bioinf/MACSMAF/scriv/macsmaf_research.scriv?id=60A8403E-C3BB-4203-9282-14712EFCD513

        **updated to use GOAtools
        GO/GOslims: Dropbox/bioinf/MACSMAF/datasets/d014
        datasets/d014/goa_map_slim_test.sh

    c) Check process_GOslims.R to give (in simplest forms) comb_summary_rank_1_* for top ranked enrichments for KEGG, REAC and GO:BP per method   
