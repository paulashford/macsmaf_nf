# Convert ID formats for genes / proteins used in network dbs
# 04 12 2024 [original version in experiments/e023 - moved to NextFlow project dir 09 12 2024]
# import pandas as pd
# import pathlib
from macsr_nf_dev.constants import (
    ALL_VALID_TYPES,
    HGNC_ESSENTIAL_COLS,
    BIOMART_ESSENTIAL_COLS,
    MAPPING_TYPE_ESSENTIAL_COLS_DICT
    # BIOMART_VALID_TYPES
)
from macsr_nf_dev.errors import (
    ParseMapFileError,
    DelimeterError
)
import pandas as pd
import logging
LOG = logging.getLogger(__name__)

def conv_ids(map_file, map_file_type, id_list, id_type, approved_only, col_filter):
    # check thd id_type
    if not id_type in ALL_VALID_TYPES:
        return(print("id_type must be one of: " + ", ".join(str(x) for x in ALL_VALID_TYPES)))
    
    # convert id_list to DataFrame
    df_id = pd.DataFrame(data = id_list, columns = [id_type], dtype=object)
    
    # Parse mapping file
    try:
        df_map = pd.read_csv(map_file, sep = '\t', dtype=object)
    except:
        raise ParseMapFileError(f"Failed to parse map file ({map_file}) as Pandas.DataFrame.")
    if df_map.columns.size == 1:
        raise DelimeterError(f"Parsing of map file ({map_file}) resulted in only 1 column - check map file is tab-delimeted.")
    
    # HGNC status filter
    if approved_only & (map_file_type == 'hgnc'):
        df_map[df_map['status'] == 'Approved']
        
    # Mapping: merge id_list and HGNC - left join will return all the passed IDs even if no matches in HGNC
    df_merge = pd.merge(df_id, df_map, on = id_type, how = "left")
    # Apply col filter
    if col_filter:
        df_merge = df_merge[MAPPING_TYPE_ESSENTIAL_COLS_DICT[map_file_type]]
    
    return(df_merge)
    
# def conv_uniprot(id_list = [], id_type = '', approved_only = True, col_filter = []):
#     # check thd id_type
#     if not id_type in HGNC_VALID_TYPES:
#         return(print("id_type must be one of: " + ", ".join(str(x) for x in HGNC_VALID_TYPES)))
    
#     # convert id list to DataFrame
#     df_id = pd.DataFrame(data = id_list, columns = [id_type])
    
#     # HGNC file
#     df_hgnc = pd.read_csv(HGNC_MAP_TABLE, sep = '\t')
#     if approved_only:
#         df_hgnc
    
#     # merge id_list and HGNC - left join will return all the passed IDs even if no matches in HGNC
#     df_merge = pd.merge(df_id, df_hgnc, on = id_type, how = "left")
#     if bool(col_filter):
#         df_merge = df_merge[col_filter]
    
#     return(df_merge)
    
    
# # <ipython-input-4-874003d7c70f>:1: DtypeWarning: Columns (32,34,38,40,50) have mixed types. Specify dtype option on import or set low_memory=False.
# <class 'pandas.core.frame.DataFrame'>
# RangeIndex: 43839 entries, 0 to 43838
# Data columns (total 54 columns):
#  #   Column                    Non-Null Count  Dtype  
# ---  ------                    --------------  -----  
#  0   hgnc_id                   43839 non-null  object 
#  1   symbol                    43839 non-null  object 
#  2   name                      43839 non-null  object 
#  3   locus_group               43839 non-null  object 
#  4   locus_type                43839 non-null  object 
#  5   status                    43839 non-null  object 
#  6   location                  43832 non-null  object 
#  7   location_sortable         0 non-null      float64
#  8   alias_symbol              22403 non-null  object 
#  9   alias_name                7620 non-null   object 
#  10  prev_symbol               12551 non-null  object 
#  11  prev_name                 22922 non-null  object 
#  12  gene_group                25863 non-null  object 
#  13  gene_group_id             25863 non-null  object 
#  14  date_approved_reserved    43826 non-null  object 
#  15  date_symbol_changed       10076 non-null  object 
#  16  date_name_changed         24488 non-null  object 
#  17  date_modified             43827 non-null  object 
#  18  entrez_id                 43753 non-null  float64
#  19  ensembl_gene_id           41201 non-null  object 
#  20  vega_id                   32325 non-null  object 
#  21  ucsc_id                   24301 non-null  object 
#  22  ena                       20351 non-null  object 
#  23  refseq_accession          42234 non-null  object 
#  24  ccds_id                   13355 non-null  object 
#  25  uniprot_ids               20294 non-null  object 
#  26  pubmed_id                 23408 non-null  object 
#  27  mgd_id                    19398 non-null  object 
#  28  rgd_id                    18941 non-null  object 
#  29  lsdb                      2151 non-null   object 
#  30  cosmic                    569 non-null    object 
#  31  omim_id                   17336 non-null  object 
#  32  mirbase                   1912 non-null   object 
#  33  homeodb                   312 non-null    float64
#  34  snornabase                400 non-null    object 
#  35  bioparadigms_slc          495 non-null    object 
#  36  orphanet                  4418 non-null   float64
#  37  pseudogene.org            8507 non-null   object 
#  38  horde_id                  856 non-null    object 
#  39  merops                    720 non-null    object 
#  40  imgt                      675 non-null    object 
#  41  iuphar                    3558 non-null   object 
#  42  kznf_gene_catalog         0 non-null      float64
#  43  mamit-trnadb              22 non-null     float64
#  44  cd                        371 non-null    object 
#  45  lncrnadb                  149 non-null    object 
#  46  enzyme_id                 1979 non-null   object 
#  47  intermediate_filament_db  0 non-null      float64
#  48  rna_central_id            8721 non-null   object 
#  49  lncipedia                 2433 non-null   object 
#  50  gtrnadb                   587 non-null    object 
#  51  agr                       40252 non-null  object 
#  52  mane_select               19235 non-null  object 
#  53  gencc                     5165 non-null   object 
# dtypes: float64(7), object(47)
# memory usage: 18.1+ MB
