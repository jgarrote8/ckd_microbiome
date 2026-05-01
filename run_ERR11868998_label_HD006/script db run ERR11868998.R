# =========================
# 1) CONFIGURACIÓN DE LA MUESTRA Y RUTAS
# =========================

sample_id <- "ERR11868998" # <--- 
sample_label <- "HD006"   # <--- 
base_dir <- "C:/Users/egt00/Desktop/TFG CKD" # <--- 

project_dir <- file.path(base_dir, paste0("metagenoma_run_", sample_id, "_label_", sample_label))

raw_dir       <- file.path(project_dir, "data_raw")
processed_dir <- file.path(project_dir, "data_processed")
results_dir   <- file.path(project_dir, "results")
scripts_dir   <- file.path(project_dir, "scripts")
docs_dir      <- file.path(project_dir, "docs")

metaphlan_dir <- file.path(results_dir, "MetaPhlAn")
humann_dir    <- file.path(results_dir, "HUMAnN")
figures_dir   <- file.path(results_dir, "figures")

# crear carpetas si no existen
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(scripts_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(docs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(metaphlan_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(humann_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 2) ARCHIVO DE ENTRADA (Análisis inicial de los datos de MetaPhlAn)
# =========================

meta_file <- file.path(metaphlan_dir, paste0(sample_id, "_metaphlan_profile.txt"))

if (!file.exists(meta_file)) {
  stop(paste("ERROR: No se encuentra el archivo", meta_file, "."))
}

meta <- read.delim(
  meta_file,
  comment.char = "#",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(meta) <- c("clade_name", "NCBI_tax_id", "relative_abundance", "additional_species")

#Añadir el nivel taxonómico
meta$tax_level <- ifelse(grepl("\\|s__", meta$clade_name), "species",
                         ifelse(grepl("\\|g__", meta$clade_name), "genus",
                                ifelse(grepl("\\|f__", meta$clade_name), "family",
                                       ifelse(grepl("\\|o__", meta$clade_name), "order",
                                              ifelse(grepl("\\|c__", meta$clade_name), "class",
                                                     ifelse(grepl("\\|p__", meta$clade_name), "phylum",
                                                            ifelse(grepl("^k__", meta$clade_name), "kingdom", "other")))))))

#Separar por nivel
phylum_df  <- subset(meta, tax_level == "phylum")
genus_df   <- subset(meta, tax_level == "genus")
species_df <- subset(meta, tax_level == "species")

phylum_df  <- phylum_df[order(-phylum_df$relative_abundance), ]
genus_df   <- genus_df[order(-genus_df$relative_abundance), ]
species_df <- species_df[order(-species_df$relative_abundance), ]

#Sacar nombre limpio del taxón
get_last_taxon <- function(x) {
  sapply(strsplit(x, "\\|"), function(y) tail(y, 1))
}

phylum_df$taxon  <- sub("p__", "", get_last_taxon(phylum_df$clade_name))
genus_df$taxon   <- sub("g__", "", get_last_taxon(genus_df$clade_name))
species_df$taxon <- sub("s__", "", get_last_taxon(species_df$clade_name))

# Limpiar SGBs (IDs de clúster genético)
species_df <- species_df[!grepl("^t__", species_df$taxon), ]

#Guardar los datos separados por filo, género, especie en CSV usando sample_id
write.csv(phylum_df,  file.path(metaphlan_dir, paste0("phylum_", sample_id, ".csv")), row.names = FALSE)
write.csv(genus_df,   file.path(metaphlan_dir, paste0("genus_", sample_id, ".csv")), row.names = FALSE)
write.csv(species_df, file.path(metaphlan_dir, paste0("species_", sample_id, ".csv")), row.names = FALSE)

#Los nombres de especie o filo pueden ser muy largos para las gráficas. Esta función ajusta las gráficas para estos
short_taxon <- function(x, max_chars = 28) {
  ifelse(nchar(x) > max_chars,
         paste0(substr(x, 1, max_chars - 3), "..."),
         x)
}

#Gráfico barplot de filos
top_phyla <- head(phylum_df, 10)

png(file.path(figures_dir, paste0("top_phyla_", sample_id, ".png")),
    width = 1400, height = 900, res = 150)

par(mar = c(5, 18, 4, 2))  

barplot(
  rev(top_phyla$relative_abundance),
  names.arg = rev(top_phyla$taxon),
  horiz = TRUE,
  las = 1,
  cex.names = 0.9,
  xlab = "Relative abundance (%)",
  main = paste("Top phyla -", sample_id)
)

dev.off()

# =========================
# 4) Resto de datos estándar y gráficas
# =========================

#-----------Top 20 especies----------
top_species <- head(species_df, 20)
top_species$label <- short_taxon(top_species$taxon, max_chars = 30)

png(file.path(figures_dir, paste0("top_species_", sample_id, ".png")), width = 1800, height = 1200, res = 180)
par(mar = c(14, 5, 4, 2))  
barplot(
  top_species$relative_abundance,
  names.arg = top_species$label,
  las = 2, cex.names = 0.8, cex.axis = 0.9, cex.lab = 1,
  ylab = "Relative abundance (%)", main = paste("Top 20 bacterial species -", sample_id)
)
dev.off()

#-------------Calcular ratio Firmicutes/Bacteroidota---------------
firmicutes <- sum(phylum_df$relative_abundance[phylum_df$taxon == "Firmicutes"])
bacteroidota <- sum(phylum_df$relative_abundance[phylum_df$taxon %in% c("Bacteroidota", "Bacteroidetes")])
FB_ratio <- firmicutes / bacteroidota

#-----------Ratio Proteobacteria (marcador de disbiosis)---------
proteobacteria <- sum(phylum_df$relative_abundance[phylum_df$taxon == "Pseudomonadota"])

#----------Diversidad taxonómica simple y Shannon---------
n_species <- nrow(species_df)
if (!require("vegan")) install.packages("vegan")
library(vegan)
shannon <- diversity(species_df$relative_abundance, index = "shannon")

#------------Guardar resumen taxonómico----
summary_table <- data.frame(
  metric = c("Number of species", "Firmicutes/Bacteroidota ratio", "Proteobacteria abundance (%)", "Shannon diversity"),
  value = c(n_species, FB_ratio, proteobacteria, shannon)
)

write.csv(summary_table, file.path(metaphlan_dir, paste0("taxonomic_summary_", sample_id, ".csv")), row.names = FALSE)

# =====================================================================================================
#---------                  ANÁLISIS CON HUMAnN (ACTUALIZADO A KOs)        ----------------------------
# -----------------------------------------------------------------------------------------------------
# =====================================================================================================

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# -------- Archivos de entrada HUMAnN (NUEVOS) --------
path_relab_file <- file.path(humann_dir, paste0(sample_id, "_pathabundance_relab.tsv"))
ko_relab_file   <- file.path(humann_dir, paste0(sample_id, "_genefamilies_ko_relab.tsv"))

if (!file.exists(path_relab_file) | !file.exists(ko_relab_file)) {
  stop("ERROR: Faltan archivos de HUMAnN. Comprueba que las rutas son correctas.")
}

# -------- Cargar tablas --------
path_relab <- read_tsv(path_relab_file, show_col_types = FALSE)
ko_relab   <- read_tsv(ko_relab_file, show_col_types = FALSE)

colnames(path_relab)[1] <- "feature"
colnames(ko_relab)[1]   <- "feature"

# -------- Función para separar filas globales y estratificadas --------
split_humann_table <- function(df) {
  list(
    global = df %>% filter(!str_detect(feature, "\\|")),
    stratified = df %>% filter(str_detect(feature, "\\|"))
  )
}

path_split <- split_humann_table(path_relab)
ko_split   <- split_humann_table(ko_relab)

# -------- Quitar columnas sin identificar --------
path_global <- path_split$global %>%
  filter(!feature %in% c("UNMAPPED", "UNINTEGRATED"))

ko_global <- ko_split$global %>%
  filter(!feature %in% c("UNMAPPED", "UNGROUPED"))

# -------- Columnas de abundancia --------
path_value_col <- colnames(path_global)[2]
ko_value_col   <- colnames(ko_global)[2]

# -------- Top rutas --------
top_pathways <- path_global %>%
  arrange(desc(.data[[path_value_col]])) %>%
  slice_head(n = 20)

write.csv(top_pathways, file.path(humann_dir, paste0("top_pathways_", sample_id, ".csv")), row.names = FALSE)

# -------- Top KOs (KEGG Orthology) --------
top_kos <- ko_global %>%
  arrange(desc(.data[[ko_value_col]])) %>%
  slice_head(n = 30)

write.csv(top_kos, file.path(humann_dir, paste0("top_KOs_", sample_id, ".csv")), row.names = FALSE)

# -------- Panel de funciones candidatas relacionadas con CKD --------
ckd_pattern <- paste(
  c("tryptophan", "tyrosine", "phenylalanine",
    "arginine", "ornithine", "urea", "urease", "nitrogen",
    "ammonia", "polyamine", "choline", "carnitine",
    "putrescine", "spermidine", "indole", "cresol", "tma", "trimethylamine"),
  collapse = "|"
)

ckd_pathways <- path_global %>%
  filter(str_detect(str_to_lower(feature), ckd_pattern)) %>%
  arrange(desc(.data[[path_value_col]]))

ckd_kos <- ko_global %>%
  filter(str_detect(str_to_lower(feature), ckd_pattern)) %>%
  arrange(desc(.data[[ko_value_col]]))

write.csv(ckd_pathways, file.path(humann_dir, paste0("ckd_candidate_pathways_", sample_id, ".csv")), row.names = FALSE)
write.csv(ckd_kos, file.path(humann_dir, paste0("ckd_candidate_KOs_", sample_id, ".csv")), row.names = FALSE)

# -------- Tabla resumen conjunta --------
ckd_pathways_summary <- ckd_pathways %>%
  transmute(tipo = "pathway", funcion = feature, abundancia = .data[[path_value_col]])

ckd_kos_summary <- ckd_kos %>%
  transmute(tipo = "KO", funcion = feature, abundancia = .data[[ko_value_col]])

ckd_summary <- bind_rows(ckd_pathways_summary, ckd_kos_summary) %>%
  arrange(desc(abundancia))

write.csv(ckd_summary, file.path(humann_dir, paste0("ckd_function_summary_", sample_id, ".csv")), row.names = FALSE)

# -------- Gráfico top rutas funcionales --------
shorten_names <- function(x, max_length = 50) {
  ifelse(nchar(x) > max_length, paste0(substr(x, 1, max_length - 3), "..."), x)
}

top_pathways$short_feature <- shorten_names(top_pathways$feature)

p_top_paths <- ggplot(top_pathways, 
                      aes(x = reorder(short_feature, .data[[path_value_col]]), 
                          y = .data[[path_value_col]])) +
  geom_col(fill = "#8DA0CB") +
  coord_flip() +
  labs(x = "Pathway", y = "Relative abundance", 
       title = paste("Top 20 functional pathways -", sample_id)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 9))

ggsave(
  filename = file.path(figures_dir, paste0("top_pathways_", sample_id, ".png")),
  plot = p_top_paths, width = 11, height = 8, dpi = 300
)

# -------- Gráfico panel CKD: rutas --------
if (nrow(ckd_pathways) > 0) {
  p1 <- ggplot(ckd_pathways %>% slice_head(n = 15),
               aes(x = reorder(feature, .data[[path_value_col]]), y = .data[[path_value_col]])) +
    geom_col(fill = "#4C72B0") +
    coord_flip() +
    labs(x = "Pathway", y = "Relative abundance", title = paste("CKD-related candidate pathways -", sample_id)) +
    theme_minimal(base_size = 11)
  
  ggsave(filename = file.path(figures_dir, paste0("ckd_candidate_pathways_", sample_id, ".png")), plot = p1, width = 10, height = 7, dpi = 300)
}

# -------- Gráfico panel CKD: KOs --------
if (nrow(ckd_kos) > 0) {
  p2 <- ggplot(ckd_kos %>% slice_head(n = 15),
               aes(x = reorder(feature, .data[[ko_value_col]]), y = .data[[ko_value_col]])) +
    geom_col(fill = "#55A868") +
    coord_flip() +
    labs(x = "KEGG Orthology (KO)", y = "Relative abundance", title = paste("CKD-related candidate KOs -", sample_id)) +
    theme_minimal(base_size = 11)
  
  ggsave(filename = file.path(figures_dir, paste0("ckd_candidate_KOs_", sample_id, ".png")), plot = p2, width = 10, height = 7, dpi = 300)
}

# =========================
# --------- Heatmaps HUMAnN ---------
# =========================

if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

if (nrow(ckd_pathways) > 0) {
  mat_ckd_path <- as.data.frame(ckd_pathways)
  rownames(mat_ckd_path) <- mat_ckd_path$feature
  mat_ckd_path$feature <- NULL
  mat_ckd_path <- as.matrix(mat_ckd_path)
  colnames(mat_ckd_path) <- "Abundancia relativa"
  
  png(file.path(figures_dir, paste0("heatmap_ckd_candidate_pathways_", sample_id, ".png")), width = 1400, height = 1800, res = 180)
  pheatmap(
    mat_ckd_path, scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
    border_color = NA, main = paste("Rutas candidatas (CKD) -", sample_id), fontsize_row = 9, fontsize_col = 11
  )
  dev.off()
}

if (nrow(ckd_kos) > 0) {
  mat_ckd_ko <- as.data.frame(ckd_kos)
  rownames(mat_ckd_ko) <- mat_ckd_ko$feature
  mat_ckd_ko$feature <- NULL
  mat_ckd_ko <- as.matrix(mat_ckd_ko)
  colnames(mat_ckd_ko) <- "Abundancia relativa"
  
  png(file.path(figures_dir, paste0("heatmap_ckd_candidate_KOs_", sample_id, ".png")), width = 1400, height = 2200, res = 180)
  pheatmap(
    mat_ckd_ko, scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
    border_color = NA, main = paste("KOs candidatos (CKD) -", sample_id), fontsize_row = 8, fontsize_col = 11
  )
  dev.off()
}

top50_kos <- ko_global %>%
  arrange(desc(.data[[ko_value_col]])) %>%
  slice_head(n = 50)

if(nrow(top50_kos) > 0) {
  mat_top50_ko <- as.data.frame(top50_kos)
  rownames(mat_top50_ko) <- mat_top50_ko$feature
  mat_top50_ko$feature <- NULL
  mat_top50_ko <- as.matrix(mat_top50_ko)
  colnames(mat_top50_ko) <- "Abundancia relativa"
  
  png(file.path(figures_dir, paste0("heatmap_top50_KOs_", sample_id, ".png")), width = 1600, height = 2400, res = 200)
  pheatmap(
    mat_top50_ko, scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
    border_color = NA, main = paste("Top 50 KOs -", sample_id), fontsize_row = 7, fontsize_col = 11
  )
  dev.off()
}