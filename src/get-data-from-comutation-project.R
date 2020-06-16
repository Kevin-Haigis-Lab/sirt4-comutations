# Get the data for this analysis from the KRAS comutation project.

source(file.path("lib", 'genesets.R'))

# Get the file path from the perspective of the comutation project home dir.
comutation_path <- function(...) { file.path("..", "comutation", ...) }

# Save a file to the data directory for this analysis.
save_data <- function(f, ...) { write_tsv(f, file.path("data", ...)) }


load(comutation_path("cache", "cancer_full_coding_muts_df.RData"))
cancer_full_coding_muts_df %>%
    filter(cancer == "COAD") %>%
    save_data("cancer_full_coding_muts_df.tsv")


load(file.path("cache", "cancer_coding_av_muts_df.RData"))
cancer_coding_av_muts_df %>%
    filter(cancer == "COAD") %>%
    save_data("mutation_data_annotated.tsv")


# NEED TO TRY TO GET CNA DATA FROM GM ORIGINAL BEASTLY OBJECT.


gm_data_object <- readRDS(
    file.path("data", 
              "cancer-data", 
              "full_data_paad_luad_coadread_skcm_mm_Jun2019.rds")
)

gm_data_object$coadread$copynumber$data %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    rename(hugo_symbol = gene_symbol,
           tumor_sample_barcode = case_id,
           cancer = tumor_type,
           cna_value = cn_avalue) %>%
    mutate(cancer = "COAD") %>%
    save_data("coad_cna_df.tsv")

coad_rna <- gm_data_object$coadread$expression$data %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    rename(hugo_symbol = gene_symbol,
           tumor_sample_barcode = case_id,
           cancer = tumor_type) %>%
    mutate(cancer = "COAD") %>%
    save_data("coad_rna_df.tsv")
