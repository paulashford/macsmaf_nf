experiments/e019/functional_enrichment_tables/README.txt
May-July 2023

Functional annotation of modules to capture primary enrichment terms from GO:BP, KEGG and Reactome.
Based on g:Profiler enrichments

Files are for:
	network: {cpdb, humanbase, string}
	method:	{M1, R1}

func_enrichments_*			Full enrichment tables with all enriched terms per module for GO/ KEGG / REAC
							[experiments/e019/create_final_output_tables.R]

func_rank1_agg_detail*		Top ranked functions only, using Rank Aggregation of p-value and MCC (Matthews correlation) rankings
							with tied-1st defaulting to the MCC rank1 (figured would likely give more specific terms)
							1-3 rows per module (GO/KEGG/REAC) for modules with sig. enriched function 
							[experiments/e019/create_final_output_tables_rank_aggregation.R]

func_rank1_agg_summary*		Top ranked functions only as above, simplified to overall function name from combining GO, KEGG and Reactome
							1 row per module reporting function_summary for each module having any sig. enriched function 
							[experiments/e019/create_final_output_tables_rank_aggregation.R]

func_rank1_agg_summary_hugo* These are just the summary files with 2 extra cols for the HGNC converted gene IDs alonside the HGNC descriptions
							Gene name fields concatenated per module with ";"
							Gene descriptions concat with " ~ "
							Note: Literal 'None' may be found in some genes_descriptions or genes_hugo meaning the gene ID used in 
							humanbase/cpdb/string could not be found/converted (eg may be old/expired)

experiments/e019/create_final_output_tables_rank_aggregation.R

uploaded to:
https://drive.google.com/drive/folders/14OPJzYnu-PwcKr7Dl8hL9ni0wzk76t5v
27 07 2023