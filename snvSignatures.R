# ------------------ Pacotes ------------------
library(data.table)
library(gridExtra)
library(NMF)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(MutationalPatterns)
library(SomaticSignatures)
library(GenomicRanges)
library(VariantAnnotation)
library(BSgenome)
library(bigmemory)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(yaml)
library(rstudioapi)
library(patchwork)


#----- Roda no RStudio e no Rscript
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  # rodando no RStudio
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  # rodando via Rscript
  script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
}

config <- yaml::read_yaml(file.path(script_dir, "config.yaml"))

vcf_dir <- file.path(script_dir, config$vcf_dir)
group_file <- file.path(vcf_dir, config$group_file)

cat("VCF dir:", vcf_dir, "\n")
cat("Group file:", group_file, "\n")

# ------------------ Diretórios de saída ------------------
outdir <- file.path(script_dir, "results")
plots_dir <- file.path(outdir, "plots")
tables_dir <- file.path(outdir, "tables")
rds_dir <- file.path(outdir, "rds")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------ Funções de salvamento ------------------
save_table <- function(df, name) {
  fwrite(df, file.path(tables_dir, paste0(name, ".tsv")), sep = "\t")
}
save_rds <- function(obj, name) {
  saveRDS(obj, file.path(rds_dir, paste0(name, ".rds")))
}
save_plot <- function(plot, name, width = 10, height = 10) {
  ggsave(filename = file.path(plots_dir, paste0(name, ".pdf")),
         plot = plot, width = width, height = height)
}
save_png <- function(plot, name, width = 10, height = 10) {
  ggsave(filename = file.path(plots_dir, paste0(name, ".png")),
         plot = plot, width = width, height = height)
} #fiz essa função para caso eu queira salvar algo em png

# ------------------ Referência ------------------
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
seqinfo_hg38 <- seqinfo(get(ref_genome))

# ------------------ Leitura do groups.txt ------------------
group_file_name <- fread(group_file, header = FALSE)
setnames(group_file_name, c("groups", "sample"))
print(group_file_name)

samples_per_group <- group_file_name %>%
  group_by(groups) %>%
  summarise(count = n_distinct(sample), .groups = "drop")

sample_names <- group_file_name$sample
sample_id <- sub("_snps.hg38_multianno.txt", "", sample_names)
groups <- group_file_name$groups

myfiles <- lapply(file.path(vcf_dir, sample_names), function(x) {
  fread(x, header = TRUE, sep = "\t", colClasses = c("Ref" = "character", "Alt" = "character"))
})

# ------------------ PER SAMPLE ------------------
sample_gr <- list()
for (i in seq_along(sample_names)) {
  file <- as.data.frame(myfiles[[i]][, 1:5])
  colnames(file) <- c("chr", "startvar", "endvar", "ref", "alt")
  
  file$sample_names <- sample_names[i]
  file$sampleID <- sample_id[i]
  file$group <- groups[i]
  
  gr <- GRanges(
    seqnames = file$chr,
    ranges = IRanges(file$startvar, file$endvar),
    ref = file$ref,
    alt = file$alt,
    seqinfo = seqinfo_hg38
  )
  
  mcols(gr)$group <- file$group
  mcols(gr)$sampleID <- file$sampleID
  mcols(gr)$sample_names <- file$sample_names
  
  sample_gr[[i]] <- gr
  names(sample_gr)[i] <- sample_id[i]
}

cat("Objetos GRanges criados para todas as amostras.\n")
sample_grList <- GRangesList(sample_gr)
names(sample_grList) <- sample_id
save_rds(sample_grList, "sample_grList")

# ------------------ PER GROUP ------------------
group_gr <- list()
unique_groups <- unique(groups)
for (g in unique_groups) {
  idx <- which(groups == g)
  file <- do.call(rbind, lapply(idx, function(i) {
    df <- as.data.frame(myfiles[[i]][, 1:5])
    colnames(df) <- c("chr", "startvar", "endvar", "ref", "alt")
    df$sampleID <- sample_id[i]
    df$group <- groups[i]
    df
  }))
  
  gr <- GRanges(
    seqnames = file$chr,
    ranges = IRanges(file$startvar, file$endvar),
    ref = file$ref,
    alt = file$alt,
    seqinfo = seqinfo_hg38
  )
  
  mcols(gr)$sampleID <- file$sampleID
  mcols(gr)$group <- file$group
  
  group_gr[[g]] <- gr
}

groups_grList <- GRangesList(group_gr)
names(groups_grList) <- as.vector(unique(groups))
save_rds(groups_grList, "groups_grList")

# ------------------ CHECK TIME ------------------
resumo_amostras <- tibble(
  sampleID = names(sample_grList),
  grupo = sapply(sample_grList, function(gr) unique(mcols(gr)$group)),
  registros = lengths(sample_grList)
)
resumo_grupos <- tibble(
  grupo = names(groups_grList),
  registros = lengths(groups_grList)  
)
soma_por_grupo <- resumo_amostras %>%
  group_by(grupo) %>%
  summarise(registros_calculados = sum(registros))
conferencia <- resumo_grupos %>%
  left_join(soma_por_grupo, by = "grupo") %>%
  mutate(
    igual = registros == registros_calculados,
    status = ifelse(igual, "OK", "Erro")
  )

save_table(resumo_amostras, "snv_resumo_amostras")
save_table(resumo_grupos, "snv_resumo_grupos")
save_table(conferencia, "snv_conferencia")

cat("\nConferência final (GRangesList x soma de amostras):", conferencia$status, "\n")

# ------------------ PER GROUPS ------------------
type_occurrences_groups <- mut_type_occurrences(sample_grList, ref_genome)
save_rds(type_occurrences_groups, "snv_type_occurrences_groups")

plot_type_occurrences_groups <- plot_spectrum(
  type_occurrences_groups, 
  CT = TRUE, by = groups, indv_points = TRUE, 
)

save_plot(plot_type_occurrences_groups, "snv_spectrum_por_grupo")

# ------------------ PER SAMPLES ------------------
type_occurrences_samples <- mut_type_occurrences(sample_grList, ref_genome)
save_rds(type_occurrences_samples, "snv_type_occurrences_samples")

plot_type_occurrences_samples <- plot_spectrum(type_occurrences_samples, CT = TRUE, by = sample_id)
save_plot(plot_type_occurrences_samples, "snv_spectrum_por_amostra")

mut_mat <- mut_matrix(vcf_list = sample_grList, ref_genome = ref_genome)
save_rds(mut_mat, "snv_mut_mat")

plot96 <- plot_96_profile(mut_mat, ymax = 0.05)
save_plot(plot96, "snv_profile_96_ct")

mut_mat_ext_context <- mut_matrix(sample_grList, ref_genome, extension = 2)
save_rds(mut_mat_ext_context, "snv_mut_mat_ext_context")

heatmap <- plot_profile_heatmap(mut_mat_ext_context, by= groups)
river <- plot_river(mut_mat_ext_context)

save_plot(heatmap, "snv_heatmap_extended")
save_plot(river, "snv_riverplot_extended")

############################## 192 contextos e Strand Bias

# Carregar pacote TxDb
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
mut_mat_strand <- mut_matrix_stranded(
  vcf_list = sample_grList,
  ref_genome = BSgenome.Hsapiens.UCSC.hg38,
  genes_hg38,
  mode = "transcription"
)
save_rds(mut_mat_strand, "snv_mut_mat_strand")

# Plot perfil 192
plot192 <- plot_192_profile(mut_mat_strand, condensed = FALSE, ymax = 0.08)
save_plot(plot192, "snv_plot_192_profile", width = 30, height = 20)

if(exists("plot96")){
  library(gridExtra)
  plot_96TR <- grid.arrange(plot96, plot192)
  save_plot(plot_96TR, "snv_plot_96e192", width = 20, height = 20)
}

# ------------------ Strand Bias ------------------
strand_counts <- strand_occurrences(mut_mat_strand, by = groups)
save_rds(strand_counts, "snv_strand_counts")

strand_bias <- strand_bias_test(strand_counts)
write.table(strand_bias, file = file.path(tables_dir, "snv_strand_bias.txt"), 
            row.names = FALSE, sep = "\t")

# Mais permissivo
strand_bias_notstrict <- strand_bias_test(strand_counts,
                                          p_cutoffs = 0.05,
                                          fdr_cutoffs = 0.1)
write.table(strand_bias_notstrict, file = file.path(tables_dir, "snv_strand_bias_notstrict.txt"),
            row.names = FALSE, sep = "\t")

# Plots de strand bias
ps1 <- plot_strand(strand_counts, mode = "relative")
ps2 <- plot_strand(strand_counts, mode = "absolute")
ps3 <- plot_strand_bias(strand_bias_notstrict, sig_type = "p")

plot_strandbias <- grid.arrange(ps1, ps2, ps3)
save_plot(plot_strandbias, "snv_plot_strandbias", width = 20, height = 25)

# ------------------ Fit para assinaturas COSMIC ------------------
cosmic_signatures <- get_known_signatures(muttype = "snv", source="COSMIC", genome = "GRCh38")

fit_res <- fit_to_signatures(mut_mat, cosmic_signatures)
save_rds(fit_res, "snv_fit_res")

plot_contribution(fit_res$contribution, coord_flip = FALSE)

contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
                                               cosmic_signatures,
                                               n_boots = 100,
                                               method = "strict")
save_rds(contri_boots, "snv_contri_boots")

plot_bootstrapped1 <- plot_bootstrapped_contribution(contri_boots)
plot_bootstrapped2 <- plot_bootstrapped_contribution(contri_boots, 
                                                     mode = "relative", 
                                                     plot_type = "dotplot")
save_plot(plot_bootstrapped1, "snv_bootstraped_absolute", width = 30, height = 25)
save_plot(plot_bootstrapped2, "snv_bootstraped_relative", width = 30, height = 25)

strict_refit <- fit_to_signatures_strict(mut_mat, cosmic_signatures, max_delta = 0.0016)
fit_res_strict <- strict_refit$fit_res
save_rds(fit_res_strict, "snv_fit_res_strict")

# Selecionar assinaturas com contribuição > 0
selected_sig2 <- which(rowSums(fit_res_strict$contribution) > 0)
prelative <- plot_contribution(fit_res_strict$contribution[selected_sig2,],
                               coord_flip = FALSE,
                               mode = "relative")
pabsolute <- plot_contribution(fit_res_strict$contribution[selected_sig2,],
                               coord_flip = FALSE,
                               mode = "absolute")
save_plot(prelative, "snv_SBS_relative", width = 15, height = 17)
save_plot(pabsolute, "snv_SBS_absolute", width = 15, height = 17)

# Plot combinado
combined_plot <- prelative | pabsolute
save_plot(combined_plot, "snv_SBS_combined", width = 30, height = 17)

# Salvar contribuição absoluta selecionada
write.table(fit_res_strict$contribution[selected_sig2,],
            file.path(tables_dir, "snv_mutationalSignaturesABSOLUTE.txt"),
            sep = "\t", row.names = TRUE)


