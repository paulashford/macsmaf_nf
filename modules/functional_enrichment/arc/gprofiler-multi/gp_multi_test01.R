# experiments/e019/gprofiler-multi/gp_multi_test01.R
# Programmatic use of g:Profiler as the full module queries are causing time-outs on web-page
# 05 05 2023

# https://www.tidyverse.org/
# install.packages("tidyverse")

# Note there is also a data.table backend to tidyverse
# https://dtplyr.tidyverse.org/
# library(data.table)
# library(dtplyr)
#  as.data.table(modlist) 
# etc

library(gprofiler2)
library(data.table)  # see note above - don't need to use here but nice to know
library(readr)
library(tidyr)
library(dplyr)
library(purrr)

mod_file <- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/modules-files/modules-R1-cpdb-coval-0.4_TEST.dat'

modules <- read_lines(mod_file)
# test2 <- parse_character(test)

# tidyr - was bit confused how to get character vector into list cols, turns out it's tidyr's:
modules_tib <- tibble( modules )

# A tidyr command for splitting cells
# As only 3 cols specified, the gene list is left alone in single col with  extra="merge"
modules_tib <- separate( modules_tib, col=1, into=c( "mod", "junk", "genes" ), sep="\t", extra="merge" )
# junk col not needed
modules_tib <- select(modules_tib, -junk)

# fixed string to add
sfixed <- "e019_M1_mod_"
# Add module string pre-number and tidy cols
modules_tib <- modules_tib %>%
    mutate( module = map_chr( modules_tib$mod, function(x) paste0(sfixed, x) ) ) %>%
    select(-mod) %>%
    relocate( genes, .after = last_col() )

# tib3 space sep genes
modules_tib <- modules_tib %>%
    rowwise() %>%
    mutate( gl = paste0( gsub( "\t", " ", genes ) )  ) %>%
    select(-genes)

# pivot the genes in gl col to rows
modlong <- modules_tib %>%
    separate_rows( gl, sep=" " ) %>%
    group_by( module )

# Nest genelists in module groups
modnest <- modlong %>%
    # group_by( module ) %>%
    nest( genes = gl)

# modnest <- modnest %>%
#     nest( genes = gl)

# for each row move nested tibbles into lists
modlist <- modnest %>%
    rowwise() %>%
    mutate( genelist = list( tibble::deframe( genes ) ) )

# return list with names of each genelist to be the module
modqry <- setNames(modlist$genelist, modlist$module)

# gprofiler
gostres4 <- gost(query = modqry,
                organism = "hsapiens", 
                ordered_query = FALSE, 
                # multi_query = TRUE, doesn't quite do what expected [see https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html]
                multi_query = FALSE, 
                significant = TRUE, 
                exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, 
                evcodes = FALSE, 
                user_threshold = 0.01, 
                correction_method = "g_SCS", 
                domain_scope = "annotated", 
                custom_bg = NULL, 
                numeric_ns = "", 
                # sources = c( "GO:MF","GO:CC","GO:BP","KEGG","REAC","HPA","CORUM","HP","WP" ), 
                sources = c( "GO:BP","KEGG","REAC" ), 
                as_short_link = FALSE)








# much scrathch!!


a=tibble::deframe(lst1$e019_M1_mod_1)

lst1 <- tibble::deframe(modnest)


#
# test2 <- test %>%
#     rowwise() %>%
#     mutate( glist = list(data))

as.list(setNames(test$data[[1]], test$module[[1]]))
# $e019_M1_mod_1
#  [1] "ENSG00000129103" "ENSG00000106392" "ENSG00000185615" "ENSG00000079150"
#  [5] "ENSG00000120942" "ENSG00000167123" "ENSG00000079459" "ENSG00000171861"
#  [9] "ENSG00000117862" "ENSG00000065057" "ENSG00000149380" "ENSG00000023318"
# [13] "ENSG00000131730" "ENSG00000162433" "ENSG00000001630" "ENSG00000141337"
# [17] "ENSG00000144455"

f <- function(x, y) setNames(y,x)
> sapply(test$data, f, y=test$module) ??


test3 <- map2( test$module, test$data, ~ setNames( .y, .x ) )


test3 <- sapply( modtib, setNames( modtib$data[[1]], modtib$module))



# nest
tib4n <- nest(tib4, gln = gl)

# tibble list
tib5 <- tib4n %>%
    group_by(module) %>%
    list( mod = module, glist = gln)


# create tibble with list cols
tib4 <- tib3 %>%
    rowwise() %>%
    tibble (mod = module, glist = gl)



tib2 <- tib2 %>%
    rowwise() %>%
    mutate(n = list(dim(genelist)))

# Convert the tabs in the genes csv style
tib2$genelist <- map( tib2$genes, function(x) paste0( "'", gsub("\t", "', '", x) , "'" ) ) 
tib2$genelist2 <- map( tib2$genes, function(x) paste0( "as.list( c('", gsub("\t", "', '", x) , "') )" ) )  
tib2$genelistev <- map( tib2$genelist2, function(x) exec(x))
# tib2$genelist2 <- map( tib2$genes, function(x) gsub("\t", " ", x) ) 
tib3 <- tib2 %>%
    select(module, genelist)

# return a df of the genelists with the module id saved in rownames
# tibbles don't like rownames so implict df conversion
df4 <- tibble::column_to_rownames( tib4, var = 'module' )
dt4 <- data.table(df4)

# lst1 <- list( list( df4$genelist ) ) - nope, just wraps whole col!
lst1 <- lapply( df4$genelist, function(x) as.list(x) ) 
lst1 <- lapply( df4$genelist, function(x) as.list(setNames()) ) 

# https://stackoverflow.com/questions/3492379/data-frame-rows-to-a-list
# lst2 <- lapply( df4$genelist, function(x) as.list( transpose(unclass( x )) ) )
lst2 <- lapply( df4$genelist, function(x) ( transpose( as.list(unclass(x)) ) ) )

names(lst1) <- rownames(df4)

# rowwise convert the "csv" style genelist to list
# note: without rowwise will put entire column into a list in each cell!
tib4 <- tib3


tib4 <- tib3 %>%
    rowwise() %>%
    mutate( genes = list( genelist ) )




tst1 <- as.list(tib3$genelist)
tst2 <- lapply(tst1, function(x) as.list(x))

names(tst1) <- tib3$module


tib4 <- tib3 %>%
    rowwise() %>%
    mutate( gl = list(genelist))

tst1 <- list(as.list(tib3)$genelist[[1]])
names(tst1) = tib3$module[[1]]

# note this converts tibble to df as tibbles don't like rownames!
# df3 <- tibble::column_to_rownames( tib3, var="module" )
df3 <- data.frame(tib3)

lst3 <- tibble::deframe(tib3)

# gprofiler
gostres3 <- gost(query = lst1,
                organism = "hsapiens", 
                ordered_query = FALSE, 
                multi_query = TRUE, 
                significant = FALSE, 
                exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, 
                evcodes = FALSE, 
                user_threshold = 0.05, 
                correction_method = "g_SCS", 
                domain_scope = "annotated", 
                custom_bg = NULL, 
                numeric_ns = "", 
                sources = c( "GO:MF","GO:CC","GO:BP","KEGG","REAC","HPA","CORUM","HP","WP" ), 
                as_short_link = FALSE)




# list cols
# tib3 <- tib2 %>%
#     rowwise() %>%
#     mutate( genelist = list(genes) )

# create list col for genes
# tib2$genes <- list(tib2$genes)


# test
# tib2 %>%
#     slice( 1 ) %>%

# map2_chr(tib2, module, genes, paste, collapse='~')

tib3 <- tib2 %>%
        rowwise() %>%
        mutate( test = list(paste0(module,genes) ) )


multiq1 <- list( "e019_M1_mod_1" = c("ENSG00000215251","ENSG00000243927","ENSG00000177646"), "e019_M1_mod_2" = c("ENSG00000174015","ENSG00000113595") )

multiq1 <- list(
		"e019_M1_mod_1" = c("ENSG00000174015","ENSG00000113595","ENSG00000127870","ENSG00000132938","ENSG00000175602","ENSG00000165630","ENSG00000178104","ENSG00000133103","ENSG00000168646","ENSG00000167110","ENSG00000135604","ENSG00000134152","ENSG00000113790","ENSG00000180938","ENSG00000146007","ENSG00000176903","ENSG00000145945","ENSG00000005238","ENSG00000077279","ENSG00000095066","ENSG00000119408","ENSG00000166407","ENSG00000140382","ENSG00000136436","ENSG00000120071","ENSG00000197256","ENSG00000160207","ENSG00000100109","ENSG00000100151","ENSG00000184507","ENSG00000171425","ENSG00000088727","ENSG00000103154","ENSG00000196150","ENSG00000114107","ENSG00000163156","ENSG00000176896","ENSG00000215271","ENSG00000171703","ENSG00000168309","ENSG00000147421","ENSG00000019995","ENSG00000076604","ENSG00000140859","ENSG00000187778","ENSG00000156795","ENSG00000270765","ENSG00000120075","ENSG00000182195","ENSG00000150636","ENSG00000123143","ENSG00000068394","ENSG00000185811","ENSG00000112578","ENSG00000198466","ENSG00000061337","ENSG00000140743","ENSG00000170264","ENSG00000103528","ENSG00000108395","ENSG00000141013","ENSG00000197905","ENSG00000128596","ENSG00000117625","ENSG00000124074","ENSG00000138100","ENSG00000138101","ENSG00000137575","ENSG00000161405","ENSG00000197114","ENSG00000080561","ENSG00000164949","ENSG00000173480","ENSG00000171847","ENSG00000104611","ENSG00000072201","ENSG00000136003","ENSG00000088808","ENSG00000066117","ENSG00000187796","ENSG00000196391","ENSG00000065491","ENSG00000188191","ENSG00000240694","ENSG00000116922","ENSG00000175513","ENSG00000087338"),
		"e019_M1_mod_2" = c("ENSG00000215251","ENSG00000243927","ENSG00000177646","ENSG00000196365","ENSG00000182180","ENSG00000164347","ENSG00000138035","ENSG00000119705","ENSG00000103202","ENSG00000061794","ENSG00000048544","ENSG00000137133","ENSG00000107815","ENSG00000131368","ENSG00000144381","ENSG00000239789","ENSG00000120662","ENSG00000101181","ENSG00000105379","ENSG00000108179","ENSG00000147586","ENSG00000156502","ENSG00000124279","ENSG00000148090","ENSG00000122033","ENSG00000085760","ENSG00000132591","ENSG00000151806","ENSG00000134905","ENSG00000122140","ENSG00000132300","ENSG00000119689","ENSG00000148187","ENSG00000174173","ENSG00000104980","ENSG00000115419","ENSG00000132463","ENSG00000172171","ENSG00000112031","ENSG00000113048","ENSG00000036473","ENSG00000139428","ENSG00000175110","ENSG00000181991","ENSG00000182199","ENSG00000244005","ENSG00000140374","ENSG00000181610","ENSG00000118246","ENSG00000125901","ENSG00000090263","ENSG00000102763","ENSG00000182810","ENSG00000125445","ENSG00000160194","ENSG00000123472","ENSG00000062582","ENSG00000116898","ENSG00000163319","ENSG00000172586","ENSG00000183010","ENSG00000155906","ENSG00000165526","ENSG00000175756","ENSG00000168924","ENSG00000156469","ENSG00000138095","ENSG00000136463","ENSG00000144029","ENSG00000137500","ENSG00000119673","ENSG00000067704","ENSG00000081177","ENSG00000108064","ENSG00000148672","ENSG00000135972","ENSG00000099821","ENSG00000156026","ENSG00000103707","ENSG00000120333","ENSG00000132676","ENSG00000132153","ENSG00000074071","ENSG00000164182","ENSG00000006744","ENSG00000168827","ENSG00000074582","ENSG00000023330","ENSG00000102738","ENSG00000050393")
	)


# test read_lines
# tr '\t' ',' < query-modules-M1-cpdb-coval-0.0-simple.dat > query-modules-M1-cpdb-coval-0.0-simple.csv
fname <- '/Users/ash/Dropbox/bioinf/MACSMAF/experiments/e019/gprofiler-multi/query-modules-M1-cpdb-coval-0.0-simple.csv'
vec <- readLines(fname)
