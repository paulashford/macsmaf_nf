HGNC_VALID_TYPES = ['entrez_id', 'ensembl_gene_id', 'hgnc_id', 'symbol', 'refseq_accession' ]
HGNC_ESSENTIAL_COLS = ['hgnc_id', 'symbol', 'locus_group', 'name', 'status', 'entrez_id', 'ensembl_gene_id', 'refseq_accession', 'alias_symbol', 'prev_symbol']
BIOMART_VALID_TYPES = ['Gene_stable_ID', 'Gene_stable_ID_version', 'Transcript_stable_ID', 'Transcript_stable_ID_version', 'Protein_stable_ID', 'Protein_stable_ID_version', 'NCBI_gene_ID', 'Gene_description', 'HGNC_ID', 'UniProtKB_ID']
BIOMART_ESSENTIAL_COLS = ['Gene_stable_ID', 'Transcript_stable_ID', 'Protein_stable_ID', 'NCBI_gene_ID', 'Gene_description', 'HGNC_ID', 'UniProtKB_ID']
UNIPROTKB_VALID_TYPES = ['UniProtKB-AC', 'UniProtKB-ID', 'EntrezGeneID', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL-CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional_PubMed']
ALL_VALID_TYPES = HGNC_VALID_TYPES + BIOMART_VALID_TYPES + UNIPROTKB_VALID_TYPES
MAPPING_FILE_TYPES = ['hgnc', 'biomart', 'uniprotkb']
MAPPING_TYPE_ESSENTIAL_COLS_DICT = {'hgnc': HGNC_ESSENTIAL_COLS, 'biomart': BIOMART_ESSENTIAL_COLS}