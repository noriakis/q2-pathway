

## Users should install Tax4Fun2 beforehand
## Tutorial: https://github.com/songweizhi/Tax4Fun2_short_tutorial
## And downloads the reference database
library(Tax4Fun2)

argv <- commandArgs(trailingOnly = TRUE)
pwd_op_folder <- argv[1] ## temp dir

## Need to specify current dir
setwd(pwd_op_folder)

query_otu_seq <- argv[2] ## rep-seqs
query_otu_table <- argv[3] ## table
pwd_ref_data <- argv[4] ## reference data
iden <- as.numeric(argv[5]) ## identity threshold
num_of_threads <- as.numeric(argv[6]) ## thread

cat(pwd_op_folder, query_otu_seq, query_otu_table, pwd_ref_data, iden, num_of_threads)

## From the tutorial code
norm_by_cn      = TRUE                          # normalize_by_copy_number (TRUE or FALSE)

## The plugin do not use pathway output though
norm_path       = TRUE                          # normalize_pathways (TRUE or FALSE)

# predict functions
runRefBlast(path_to_otus = query_otu_seq, path_to_reference_data = pwd_ref_data, path_to_temp_folder = ".", database_mode = "Ref99NR", use_force = T, num_threads = num_of_threads)
makeFunctionalPrediction(path_to_otu_table = query_otu_table, path_to_reference_data = pwd_ref_data, path_to_temp_folder = ".", database_mode = "Ref99NR", normalize_by_copy_number = norm_by_cn, min_identity_to_reference = iden, normalize_pathways = norm_path)
