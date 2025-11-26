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
group_file <- file.path(vcf_dir, config$indel_file)

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
  ggsave(filename = file.path(plots_dir, paste0(name, ".png")),
         plot = plot, width = width, height = height)
}

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
sample_id <- sub("_unique_indels.hg38_multianno.txt", "", sample_names)
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
save_rds(sample_grList, "indel_sample_grList")

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
save_rds(groups_grList, "indel_groups_grList")

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

save_table(resumo_amostras, "indel_resumo_amostras")
save_table(resumo_grupos, "indel_resumo_grupos")
save_table(conferencia, "indel_conferencia")

cat("\nConferência final (GRangesList x soma de amostras):", conferencia$status, "\n")

indel_grl <- get_mut_type(sample_grList, type = "indel")
indel_grl <- get_indel_context(indel_grl, ref_genome)
head(indel_grl[[1]], n = 5)

indel_counts <- count_indel_contexts(indel_grl)
head(indel_counts)

plot_indel <- plot_indel_contexts(indel_counts, condensed = FALSE, same_y = TRUE)
plot_indel <- plot_indel +
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 3, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    strip.text.x = element_text(size = 20, face = "plain", color = "black"), # tipo de sub
    strip.text.y = element_text(size = 20, face = "plain", color = "black") #samples
  )
plot_indel
save_plot(plot_indel, "indel_condensed", width = 25, height = 7)

plot_main_contexts <- plot_main_indel_contexts(indel_counts)
plot_main_contexts <- plot_main_contexts +
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 3, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    strip.text.x = element_text(size = 20, face = "bold", color = "black"), # tipo de sub
    strip.text.y = element_text(size = 20, face = "bold", color = "black") #samples
  )

save_plot(plot_main_contexts, "indel_main_context", height = 10, width = 15)

signatures_indel = get_known_signatures(muttype = "indel")
fit_res_indel <- fit_to_signatures(indel_counts, signatures_indel)
save_rds(fit_res_indel, "fit_res_indel")

selected_sigs_indel <- which(rowSums(fit_res_indel$contribution) > 0)
signatures_indel_reduced <- signatures_indel[, selected_sigs_indel]

strict_refit_indel <- fit_to_signatures_strict(
  mut_mat = indel_counts,
  signatures = signatures_indel_reduced,
  max_delta = 0.03,
  method = "best_subset"
)

fit_res_indel_strict <- strict_refit_indel$fit_res
indel_absolute <- plot_contribution(fit_res_indel_strict$contribution, mode = "absolute")
indel_absolute <- indel_absolute +
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 3, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    strip.text.x = element_text(size = 20, face = "bold", color = "black"), # tipo de sub
    strip.text.y = element_text(size = 20, face = "bold", color = "black") #samples
  )

indel_relative <-plot_contribution(fit_res_indel_strict$contribution, mode = "relative")
indel_relative <- indel_relative +
  theme(
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 3, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    strip.text.x = element_text(size = 20, face = "bold", color = "black"), # tipo de sub
    strip.text.y = element_text(size = 20, face = "bold", color = "black") #samples
  )


save_plot(indel_absolute, "indel_DB_absolute", width = 7, height = 15)
save_plot(indel_relative, "indel_DB_relative", width = 7, height = 15)

# Plot combinado
combined_plot <- indel_relative | indel_absolute
save_plot(combined_plot, "indel_DB_combined", width = 30, height = 17)
