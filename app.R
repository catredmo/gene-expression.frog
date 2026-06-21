# Xenopus developmental proteomics browser (Van Itallie 2025 reanalysis)
#
# Protein-level data: FragPipe combined output (Rep A + Rep B, TMT 11-plex)
#   ../cluster_reanalysis/fragpipe/tmt-report/abundance_gene_MD.tsv
#     14,540 proteins x (RepA: 11 stages, RepB: 11 stages)
#     MSstatsTMT-normalized log2 abundances (relative to channel-median ref)
#
# Phosphosite data: Rep A+B re-search against standard v10.1 DB
#   ../cluster_reanalysis/results_phospho/ -> all_phosphosites_perstage.tsv
#     25,574 phosphosites across the detectable proteome (delta-score localised)
#
# Run locally:
#   shiny::runApp("shiny_keratin_browser")

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(readr)
# Cross-species alignment tab dependencies
library(Biostrings)
library(msa)
library(ggnewscale)

# ---- File locations ----
# The app supports two layouts:
#   1) Development:        data files live in ../cluster_reanalysis/...
#   2) Self-contained:     all data files sit next to app.R (preferred for deployment)

resolve_data_file <- function(filename, dev_subpath) {
    # Self-contained first (next to app.R), then development location
    if (file.exists(filename)) return(filename)
    dev_path <- file.path("..", "cluster_reanalysis", dev_subpath, filename)
    if (file.exists(dev_path)) return(dev_path)
    stop("Cannot find data file '", filename, "'. ",
         "Place it next to app.R, or place app.R in shiny_keratin_browser/ ",
         "alongside the original cluster_reanalysis/ directory.")
}

ABUND_TSV       <- resolve_data_file("abundance_gene_MD.tsv",
                                     "fragpipe/tmt-report")
PEPS_TSV        <- resolve_data_file("keratin_peptides_all_perstage.tsv", ".")
ALL_PHOS_TSV    <- resolve_data_file("all_phosphosites_perstage.tsv", ".")
RNA_DEV_TSV     <- resolve_data_file("session2016_dev_tpm.tsv", ".")
RNA_TISSUE_TSV  <- resolve_data_file("session2016_tissue_tpm.tsv", ".")

DEV_STAGES_RNA <- c("egg", "st08", "st09", "st10", "st12", "st15",
                    "st20", "st25", "st30", "st35", "st40")
TISSUES_RNA <- c("brain", "heart", "intestine", "kidney", "liver", "lung",
                 "muscle", "ovary", "pancreas", "skin", "spleen", "stomach",
                 "testis")

# Cross-modality stage mapping for the combined view.
# Only stages where the proteomics TMT channel and the Session 2016 RNA-seq
# stage are at the same (or near-identical) Nieuwkoop-Faber stage. Stages
# without a confident pair (Oocyte_VI, St18/St22/St24/St26, St46) are omitted
# rather than wrongly aligned.
MATCHED_STAGES <- data.frame(
    pair_label = c("Egg",  "St9",   "St12",  "St30",  "St40/41"),
    prot_stage = c("Egg",  "St9",   "St12",  "St30",  "St41"),
    rna_stage  = c("egg",  "st09",  "st12",  "st30",  "st40"),
    stringsAsFactors = FALSE
)

STAGES <- c("Oocyte_VI", "Egg", "St9", "St12", "St18",
            "St22", "St24", "St26", "St30", "St41", "St46")
CHANNELS <- c("126", "127N", "127C", "128N", "128C",
              "129N", "129C", "130N", "130C", "131N", "131C")
INT_COLS <- paste0("intensity_", CHANNELS, "_", STAGES)

# ---- Helper: pull gene name from RefSeq:NP_XXXX|gene|... index strings ----
extract_gene <- function(s) {
    parts <- strsplit(s, "|", fixed = TRUE)
    vapply(parts, function(p) if (length(p) >= 2) p[2] else p[1], character(1))
}

# ---- Load protein-level abundance ----
abund_raw <- read_tsv(ABUND_TSV, show_col_types = FALSE)

# The TSV repeats stage column names (RepA stages then RepB stages). Rename
# columns positionally: first 5 = metadata, next 11 = RepA, next 11 = RepB.
stopifnot(ncol(abund_raw) == 5 + length(STAGES) * 2)
names(abund_raw)[6:16]  <- paste0("RepA_", STAGES)
names(abund_raw)[17:27] <- paste0("RepB_", STAGES)
abund_raw$gene <- extract_gene(abund_raw$Index)
abund_raw <- abund_raw %>%
    relocate(gene, .before = Index) %>%
    filter(!is.na(gene) & gene != "") %>%
    # Some gene symbols appear on >1 protein accession (a curated RefSeq plus a
    # redundant predicted model mapping to the same gene). Keep only the
    # best-supported row per symbol (most PSMs, ties broken by ReferenceIntensity)
    # so each gene plots as a single trajectory instead of two overlaid lines.
    # ~267 of ~14,240 symbols are affected; e.g. pkp3.S had a 30-PSM XP_ model
    # and a 3-PSM NP_ model under the same label.
    arrange(desc(NumberPSM), desc(ReferenceIntensity)) %>%
    distinct(gene, .keep_all = TRUE)

genes_all <- sort(unique(abund_raw$gene))

abund_long <- abund_raw %>%
    pivot_longer(matches("^Rep[AB]_"),
                 names_to = c("rep", "stage"),
                 names_pattern = "^(Rep[AB])_(.+)$",
                 values_to = "log2_abundance") %>%
    mutate(stage = factor(stage, levels = STAGES))

# ---- Load phosphosite data (Rep A+B against standard v10.1 DB) ----
# This is the single phosphosite dataset in the app. The former Rep C keratin-only
# tab (corrected krt12.4.S DB) was retired: it contained no unique phosphosites
# (its krt12.4.S S409 is the same residue as this dataset's S349, just in the
# 435-aa corrected vs 375-aa standard coordinate frame), and the restored
# krt12.4.S head is unphosphorylated, so the corrected DB exposes no new sites.
all_phos_raw <- read_tsv(ALL_PHOS_TSV, show_col_types = FALSE)
all_phos_raw <- all_phos_raw %>%
    # Renumber krt12.4.S to the corrected 435-aa frame (standard +60) so its
    # phosphosite positions match the corrected gene model used elsewhere.
    mutate(site_position = if_else(gene == "krt12.4.S",
                                   site_position + 60, site_position)) %>%
    arrange(desc(best_delta_score)) %>%
    distinct(gene, site_position, residue, .keep_all = TRUE) %>%
    mutate(confidence_tier = case_when(
        best_delta_score >= 0.3 ~ "high",
        best_delta_score >= 0.1 ~ "medium",
        TRUE ~ "low"
    ))

all_phos_long <- all_phos_raw %>%
    mutate(label = paste0(gene, " ", residue, site_position)) %>%
    pivot_longer(all_of(INT_COLS),
                 names_to = "stage_col", values_to = "intensity") %>%
    mutate(stage = factor(sub("intensity_[^_]+_", "", stage_col), levels = STAGES))

# ---- Load Session 2016 RNA-seq TPM (pre-computed) ----
rna_dev <- read_tsv(RNA_DEV_TSV, show_col_types = FALSE)
rna_tissue <- read_tsv(RNA_TISSUE_TSV, show_col_types = FALSE)

rna_dev_long <- rna_dev %>%
    pivot_longer(all_of(DEV_STAGES_RNA), names_to = "stage", values_to = "tpm") %>%
    mutate(stage = factor(stage, levels = DEV_STAGES_RNA))
rna_tissue_long <- rna_tissue %>%
    pivot_longer(all_of(TISSUES_RNA), names_to = "tissue", values_to = "tpm") %>%
    mutate(tissue = factor(tissue, levels = TISSUES_RNA))

# Union of gene names available in any data source (for the picker)
rna_genes <- unique(c(rna_dev$Geneid, rna_tissue$Geneid))
all_phos_genes <- unique(all_phos_raw$gene)
genes_union <- sort(unique(c(genes_all, rna_genes, all_phos_genes)))

# ---- Homeolog merging helpers ----
# X. laevis genes carry a .L / .S subgenome suffix. Stripping it gives the base
# symbol so a homeolog pair (e.g. krt8.1.L + krt8.1.S) can be collapsed to one.
homeolog_base <- function(g) sub("[.][LS]$", "", g)
# Combine MSstatsTMT log2 abundances by summing in linear space (total protein
# output), returning to log2. NA if every homeolog is missing in the group.
log2_sum <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else log2(sum(2^x))
}

# ---- Cross-species alignment tab (defines alignment_tab_ui / alignment_tab_server) ----
source("alignment_tab.R", local = FALSE)

# ---- UI ----
ui <- fluidPage(
    tags$head(tags$title("Xenopus laevis gene expression explorer")),
    titlePanel(HTML("<em>Xenopus laevis</em> gene expression explorer")),
    tabsetPanel(
      id = "top_tabs",
      tabPanel("Gene expression",
    p(em("Mass spectrometry: Van Itallie 2025 reanalysis (FragPipe Rep A+B) | ",
         "RNA-seq: Session et al. 2016 reanalysis (TPM, developmental + adult tissues)")),
    p(actionLink("show_docs",
                 HTML("&#x1F4D6; <strong>Data sources &amp; methods</strong> (click to expand)"),
                 style = "font-size: 14px;")),
    sidebarLayout(
        sidebarPanel(
            width = 3,
            h4("Gene selection"),
            selectizeInput("genes", "Gene(s):",
                           choices = NULL,
                           multiple = TRUE),
            helpText(sprintf("Choose from %s genes (proteomics + RNA-seq union). Start typing the gene name.",
                             format(length(genes_union), big.mark=","))),
            actionButton("preset_krt", "Keratin family", class = "btn-sm"),
            actionButton("preset_dsm", "Desmosome", class = "btn-sm"),
            actionButton("preset_ajc", "AJC", class = "btn-sm"),
            br(), br(),
            textAreaInput("genes_paste", "Paste a gene list:",
                          placeholder = "krt8.1.L, krt18.1.L, dsp.L  (comma, space, or newline separated)",
                          rows = 2),
            actionButton("genes_paste_add", "Add pasted genes", class = "btn-sm"),
            br(), br(),
            radioButtons("rep_mode", "Replicate display:",
                         choices = c("Both reps overlaid" = "both",
                                     "Mean +/- range" = "mean",
                                     "RepA only" = "RepA",
                                     "RepB only" = "RepB"),
                         selected = "both"),
            radioButtons("homeolog_mode", "Homeologs (.L / .S):",
                         choices = c("Show both homeologs" = "both",
                                     "Merge (sum .L + .S)" = "merge"),
                         selected = "both"),
            helpText(em("Merge collapses each .L/.S pair to its base symbol, ",
                        "summing both subgenomes (total gene output). Pulls both ",
                        "homeologs even if only one is selected. Applies to plots ",
                        "and heatmaps; tables stay per-homeolog.")),
            hr(),
            h4("Phosphosite tab filters"),
            helpText(em("Applies to the Phosphosites trajectory and table tabs ",
                        "(Rep A+B against the standard v10.1 DB).")),
            checkboxGroupInput("tiers", "Confidence tier:",
                               choices = c("high", "medium", "low"),
                               selected = c("high", "medium")),
            sliderInput("min_psms", "Minimum num_psms:",
                        min = 1, max = max(all_phos_raw$num_psms), value = 1),
            hr(),
            h4("Plot customization"),
            selectInput("palette", "Categorical palette (lines / bars):",
                        choices = c("Default (ggplot)" = "default",
                                    "Set1 (RColorBrewer)" = "Set1",
                                    "Set2 (RColorBrewer)" = "Set2",
                                    "Dark2 (RColorBrewer)" = "Dark2",
                                    "Paired" = "Paired",
                                    "Viridis (discrete)" = "viridis",
                                    "Plasma (discrete)" = "plasma"),
                        selected = "default"),
            selectInput("point_shape", "Point shape:",
                        choices = c("Filled circle" = "16",
                                    "Open circle" = "1",
                                    "Filled square" = "15",
                                    "Open square" = "0",
                                    "Filled triangle" = "17",
                                    "Open triangle" = "2",
                                    "Filled diamond" = "18",
                                    "Open diamond" = "5",
                                    "Plus" = "3",
                                    "Cross" = "4",
                                    "Asterisk" = "8"),
                        selected = "16"),
            sliderInput("point_size", "Point size:",
                        min = 0.5, max = 6, value = 2.5, step = 0.5),
            sliderInput("line_width", "Line width:",
                        min = 0.3, max = 3, value = 0.8, step = 0.1),
            selectInput("heatmap_palette", "Heatmap palette:",
                        choices = c("Red-blue (diverging)" = "RdBu",
                                    "Red-yellow-blue" = "RdYlBu",
                                    "Purple-orange" = "PuOr",
                                    "Viridis" = "viridis",
                                    "Magma" = "magma",
                                    "Plasma" = "plasma",
                                    "Cividis" = "cividis"),
                        selected = "RdBu")
        ),
        mainPanel(
            width = 9,
            tabsetPanel(
                id = "top_tabs", type = "tabs",
                tabPanel("Mass spectrometry",
                         br(),
                         tabsetPanel(
                             id = "ms_tabs", type = "pills",
                             tabPanel("Protein abundance",
                                      br(),
                                      plotlyOutput("plot_protein", height = "520px"),
                                      br(),
                                      p("Log2 normalized abundance per developmental stage. ",
                                        "Two TMT 11-plex replicates (RepA, RepB) overlaid by default. ",
                                        "Values are MSstatsTMT-normalized log2 ratios relative to the channel median.")),
                             tabPanel("Protein heatmaps",
                                      br(),
                                      h4("Z-scored (within gene, across stages)"),
                                      plotlyOutput("plot_protein_hm", height = "500px"),
                                      br(),
                                      h4("Raw log2 abundance (MSstatsTMT-normalized)"),
                                      plotlyOutput("plot_protein_hm_raw", height = "500px"),
                                      br(),
                                      p("Top: row-Z-scored so temporal patterns are comparable across proteins of different absolute abundance. ",
                                        "Bottom: raw log2-normalized abundance (replicate-averaged) preserving absolute level. ",
                                        "Diverging palettes work for both; switch to a sequential palette in the sidebar for raw values if preferred.")),
                             tabPanel("Phosphosites trajectory",
                                      br(),
                                      plotlyOutput("plot_all_phos", height = "520px"),
                                      br(),
                                      p("Phosphosite TMT intensity across stages. ",
                                        strong("Source:"), " Rep A+B Comet search against the standard Xenbase v10.1 database. ",
                                        "25,574 unique phosphosites across the detectable proteome. ",
                                        "Filterable by confidence tier and PSM count. ",
                                        em("Hover a point for the AScore site-localization and the peptide motif "),
                                        em("(the modified residue is lower-case; see the table for full details)."))),
                             tabPanel("Protein table",
                                      br(),
                                      DTOutput("table_protein")),
                             tabPanel("Phosphosite table",
                                      br(),
                                      DTOutput("table_all_phos"),
                                      br(),
                                      p(em("Phosphosite table (Rep A+B). Filtered by sidebar gene selection plus the phosphosite tier and PSM filters. "),
                                        em("'motif (+/-7)' shows the sequence window with the phosphoresidue in lower-case; "),
                                        em("'candidate S/T/Y' is how many residues in the detected peptide could carry the phosphate; "),
                                        em("'AScore'/'localization' give the site-localization confidence (confident >=19, likely 13-19, ambiguous <13). "),
                                        em("Sites without AScore show 'not assessed'.")))
                         )),
                tabPanel("RNA-seq",
                         br(),
                         tabsetPanel(
                             id = "rna_tabs", type = "pills",
                             tabPanel("Developmental time course",
                                      br(),
                                      plotlyOutput("plot_rna_dev", height = "420px"),
                                      br(),
                                      h4("Z-scored (within gene, across stages)"),
                                      plotlyOutput("plot_rna_dev_hm", height = "340px"),
                                      br(),
                                      h4("Raw log10(TPM+1)"),
                                      plotlyOutput("plot_rna_dev_hm_raw", height = "340px"),
                                      br(),
                                      p("Top: line plot of log10(TPM+1) across stages. ",
                                        "Middle: row-Z-scored heatmap so temporal shape is ",
                                        "comparable across genes of different absolute expression. ",
                                        "Bottom: raw log10(TPM+1) heatmap preserving absolute level. ",
                                        "Session et al. 2016 reanalysis.")),
                             tabPanel("Adult tissues",
                                      br(),
                                      plotlyOutput("plot_rna_tissue", height = "420px"),
                                      br(),
                                      h4("Z-scored (within gene, across tissues)"),
                                      plotlyOutput("plot_rna_tissue_hm", height = "340px"),
                                      br(),
                                      h4("Raw log10(TPM+1)"),
                                      plotlyOutput("plot_rna_tissue_hm_raw", height = "340px"),
                                      br(),
                                      p("Top: bar chart of raw TPM. ",
                                        "Middle: Z-scored tissue specificity. ",
                                        "Bottom: raw log10(TPM+1) for absolute comparison across tissues.")),
                             tabPanel("Developmental table",
                                      br(),
                                      DTOutput("table_rna_dev"),
                                      br(),
                                      p("TPM per gene across 11 developmental stages. ",
                                        "Rows filtered by the gene selection in the sidebar. ",
                                        "Use the Copy / CSV / Excel buttons to download.")),
                             tabPanel("Adult tissue table",
                                      br(),
                                      DTOutput("table_rna_tissue"),
                                      br(),
                                      p("TPM per gene across 13 adult tissues. ",
                                        "Rows filtered by the gene selection in the sidebar."))
                         )),
                tabPanel("Combined (matched stages)",
                         br(),
                         p(em("Cross-modality view at developmental stages where the ",
                              "proteomics TMT channel and Session 2016 RNA-seq sample ",
                              "are at the same (or near-identical) stage:"), br(),
                           strong("Egg (Egg / egg), St9 (St9 / st09), St12 (St12 / st12), ",
                                  "St30 (St30 / st30), St40/41 (St41 / st40)"), br(),
                           "Both modalities are Z-scored across the matched stages, ",
                           "so trajectory shape is comparable but absolute level is not."),
                         br(),
                         plotlyOutput("plot_combined", height = "520px"),
                         br(),
                         h4("Z-scored (within gene, within modality)"),
                         plotlyOutput("plot_combined_hm", height = "360px"),
                         br(),
                         h4("Raw values (protein: log2 abundance; RNA: log10 TPM+1)"),
                         plotlyOutput("plot_combined_hm_raw", height = "360px"),
                         br(),
                         p("Top: paired line plot (solid = protein log2 abundance, ",
                           "dashed = RNA log10 TPM, both Z-scored per gene). ",
                           "Middle: side-by-side Z-score heatmaps. ",
                           "Bottom: side-by-side raw heatmaps; modalities use different units so each facet has its own color scale."))
            )
        )
    )
    ),
    tabPanel("Cross-species alignment", alignment_tab_ui)
    )
)

# ---- Server ----
server <- function(input, output, session) {

    # ---- Cross-species alignment tab ----
    alignment_tab_server(input, output, session)

    # ---- Documentation modal ----
    observeEvent(input$show_docs, {
        showModal(modalDialog(
            title = "Data sources & methods",
            size = "l",
            easyClose = TRUE,
            footer = modalButton("Close"),

            tags$h4("Mass spectrometry"),
            tags$p(strong("Source:"), " Van Itallie ES ", em("et al."),
                   " (2025) developmental proteomics of ", em("Xenopus laevis"), ". ",
                   "Raw data: PRIDE accession ", strong("PXD060481"), "."),
            tags$p(strong("Design:"), " three TMT 11-plex SPS-MS3 replicates ",
                   "(Rep A, Rep B, Rep C) sampling 11 Nieuwkoop-Faber stages: ",
                   "Oocyte_VI, Egg, St9, St12, St18, St22, St24, St26, St30, St41, St46. ",
                   "Static mods: TMT11-plex on N-term/K, NEM on Cys. ",
                   "Variable mods: Met oxidation, Asn deamidation, ",
                   "phospho on S/T/Y (phospho searches only)."),
            tags$p(strong("Reanalysis pipeline (this app uses both outputs):")),
            tags$ul(
                tags$li(strong("Protein abundance"), " (Rep A + Rep B): ",
                        "FragPipe / MSFragger search against ",
                        em("Xenopus laevis"), " Xenbase v10.1 with corrected keratin ",
                        "FASTA, MSstatsTMT normalisation. ",
                        "Source file: ", tags$code("abundance_gene_MD.tsv"),
                        " (14,241 proteins, two replicates, log2 normalised)."),
                tags$li(strong("Phosphosites (Rep A+B): "),
                        "Comet search against the standard Xenbase v10.1 ",
                        "database. ",
                        tags$strong("25,574 unique phosphosites"),
                        " across the detectable proteome. FDR 1% PSM-level, ",
                        "phospho positions localised by delta-score ",
                        "(rank1 - rank2 xcorr). Use the 'Phosphosites trajectory' ",
                        "and 'Phosphosite table' tabs. ",
                        em("Note:"), " a separate Rep C search against the ",
                        "krt12.4.S head-restored corrected FASTA was retired - it ",
                        "yielded no unique phosphosites (the restored head is ",
                        "unphosphorylated). The one krt12.4.S phosphosite has been ",
                        "renumbered to the corrected 435-aa frame (+60) here for ",
                        "consistency with the rest of the app."),
                tags$li(strong("Site localization (AScore): "),
                        "phosphosites are annotated with ", em("AScore"),
                        " site-localization where available (", tags$code("ascore"),
                        ", ", tags$code("loc_tier"), ": confident >=19 ~p<0.01, ",
                        "likely 13-19 ~p<0.05, ambiguous <13), plus the number of ",
                        "candidate S/T/Y in the detected peptide and a +/-7 sequence ",
                        "motif (phosphoresidue lower-case). ",
                        tags$strong("Coverage caveat:"), " AScore was run on a ",
                        "subset of phospho-PSMs, so ~24% of sites (6,025) carry a ",
                        "localization score; the rest show 'not assessed'. Many sites ",
                        "- especially in the serine-rich keratin head/tail domains - ",
                        "are genuinely ", em("ambiguous"), ": the phosphate cannot be ",
                        "pinned to a single residue among several adjacent S/T.")
            ),
            tags$p(strong("Key correction applied:"), " the Xenbase v10.1 ",
                   "annotation of krt12.4.S (NP_001079456.1, 375 aa) was found to ",
                   "be truncated. The full-length 435-aa form (60-aa N-terminal ",
                   "head restored from the alternative ATG at Chr9_10S:1162844) ",
                   "was used in the reanalysis. MS data directly confirms the ",
                   "full-length form (15+ PSMs across the head N-term peptide; ",
                   "zero PSMs for the truncated-form N-term)."),

            tags$hr(),

            tags$h4("RNA-seq"),
            tags$p(strong("Source:"), " Session AM ", em("et al."),
                   " (2016) ", em("Genome evolution in the allotetraploid frog "),
                   em("Xenopus laevis"), ". ", em("Nature"), " 538, 336-343. ",
                   "GEO accession ", strong("GSE73430"),
                   " / BioProject PRJNA313977."),
            tags$p(strong("Developmental time course (11 samples):"),
                   " egg, st08, st09, st10, st12, st15, st20, st25, st30, st35, st40. ",
                   "Single-replicate per stage."),
            tags$p(strong("Adult tissue panel (13 samples):"),
                   " brain, heart, intestine, kidney, liver, lung, muscle, ovary, ",
                   "pancreas, skin, spleen, stomach, testis."),
            tags$p(strong("Reanalysis pipeline:"),
                   " HISAT2 alignment to Xenbase v10.1 ", em("X. laevis"),
                   " genome, featureCounts quantification, ",
                   "TPM normalisation per sample. ",
                   "Source files: ", tags$code("session2016_dev_tpm.tsv"),
                   " (42,675 genes × 11 stages), ",
                   tags$code("session2016_tissue_tpm.tsv"),
                   " (42,675 genes × 13 tissues)."),

            tags$hr(),

            tags$h4("Notes / caveats"),
            tags$ul(
                tags$li(strong("Nomenclature:"),
                        " The proteomics keratin gene symbols have been harmonised ",
                        "to match the post-2018 Xenbase RNA-seq nomenclature: ",
                        tags$code("krt8.L"), " -> ", tags$code("krt8.1.L"), ", ",
                        tags$code("krt8.S"), " -> ", tags$code("krt8.1.S"), ", ",
                        tags$code("krt18.L"), " -> ", tags$code("krt18.1.L"), ", ",
                        tags$code("krt18.S"), " -> ", tags$code("krt18.1.S"), ", ",
                        tags$code("krt15.S"), " -> ", tags$code("krt15.1.S"), ", ",
                        tags$code("krt78.L"), " -> ", tags$code("krt78.1.L"),
                        ". The entry labelled ", tags$code("krt12.L"),
                        " is deliberately kept under its legacy symbol because the ",
                        "Xenbase v10.1 gene model is a corrupted prediction (a ",
                        "truncated krt23.L missing the N-terminal 38 aa) - treat with caution."),
                tags$li(strong("Stage matching:"),
                        " the Combined tab pairs only stages with confident ",
                        "Nieuwkoop-Faber correspondence (Egg, St9, St12, St30, St40/41). ",
                        "Intermediate proteomics stages (St18, St22, St24, St26, St46) ",
                        "have no exact RNA-seq counterpart."),
                tags$li(strong("Phosphosite source:"),
                        " phosphosites come from the Rep A+B Comet search against ",
                        "the standard v10.1 FASTA (25,574 sites). The earlier ",
                        "Rep C keratin-only search (krt12.4.S-corrected FASTA) was ",
                        "retired as redundant - it added no unique sites because the ",
                        "restored krt12.4.S head is unphosphorylated. krt12.4.S ",
                        "positions are reported in the corrected 435-aa frame."),
                tags$li(strong("Redundant protein models:"),
                        " ~267 gene symbols in the abundance data appeared on more ",
                        "than one protein accession (a curated RefSeq plus a redundant ",
                        "predicted model mapping to the same gene). For each symbol the ",
                        "app shows only the best-supported accession (most PSMs), so ",
                        "every gene plots as a single trajectory. e.g. ",
                        tags$code("pkp3.S"), " uses the 30-PSM ", tags$code("XP_018115036.1"),
                        " model rather than the 3-PSM ", tags$code("NP_001084424.1"), " one."),
                tags$li(strong("Homeolog merging:"),
                        " the sidebar 'Merge (sum .L + .S)' option collapses each ",
                        "homeolog pair to its base symbol, summing the two subgenomes ",
                        "in linear space (protein: ", tags$code("log2(2^L + 2^S)"),
                        "; RNA: summed TPM) to represent total gene output. It pulls ",
                        "both homeologs from the full data even if only one is selected, ",
                        "and applies to the plots and heatmaps only - the data tables ",
                        "always stay per-homeolog.")
            ),

            tags$hr(),

            tags$h4("Source code"),
            tags$p("The app and its underlying data TSVs are version-controlled at ",
                   tags$a(href = "https://github.com/catredmo/gene-expression.frog",
                          target = "_blank",
                          tags$code("github.com/catredmo/gene-expression.frog")),
                   ". Pull requests and issues welcome."),

            tags$p(em("Last updated: 2026-06-16. Data is subject to revision as ",
                      "the underlying gene model corrections progress."))
        ))
    })

    # Populate gene selector server-side (efficient with 40k+ genes)
    updateSelectizeInput(session, "genes",
                         choices = genes_union,
                         selected = c("krt19.L", "krt19.S"),
                         server = TRUE)

    observeEvent(input$preset_krt, {
        krt <- grep("^krt", genes_union, value = TRUE)
        updateSelectizeInput(session, "genes", selected = krt[1:min(20, length(krt))],
                             choices = genes_union, server = TRUE)
    })
    observeEvent(input$preset_dsm, {
        dsm <- intersect(genes_union,
                         c("dsg2.L","dsg2.S","dsg3.L","dsg3.S","dsc1.L","dsc1.S",
                           "dsc2.L","dsc2.S","dsp.L","dsp.S","pkp1.L","pkp1.S",
                           "pkp3.L","pkp3.S","jup.L","jup.S","pgr.L","pgr.S"))
        updateSelectizeInput(session, "genes", selected = dsm,
                             choices = genes_union, server = TRUE)
    })
    observeEvent(input$preset_ajc, {
        ajc <- intersect(genes_union,
                         c("cdh1.L","cdh1.S","cdh3.L","cdh3.S","ctnna1.L","ctnna1.S",
                           "ctnnb1.L","ctnnb1.S","ctnnd1.L","ctnnd1.S","tjp1.L","tjp1.S",
                           "tjp2.L","tjp2.S","jam2.L","jam2.S","cldn1.L","cldn1.S",
                           "ocln.L","ocln.S","afdn.L","afdn.S"))
        updateSelectizeInput(session, "genes", selected = ajc,
                             choices = genes_union, server = TRUE)
    })
    observeEvent(input$genes_paste_add, {
        req(input$genes_paste)
        # Split on commas, semicolons, or any whitespace (incl. newlines/tabs).
        toks <- trimws(strsplit(input$genes_paste, "[,;[:space:]]+")[[1]])
        toks <- toks[toks != ""]
        req(length(toks) > 0)
        # Case-insensitive match so e.g. "KRT8.1.L" resolves to "krt8.1.L".
        matched <- genes_union[tolower(genes_union) %in% tolower(toks)]
        unmatched <- toks[!tolower(toks) %in% tolower(genes_union)]
        # Append to the current selection rather than replacing it.
        sel <- union(input$genes, matched)
        updateSelectizeInput(session, "genes", selected = sel,
                             choices = genes_union, server = TRUE)
        if (length(unmatched) > 0) {
            showNotification(
                paste0("Not found (", length(unmatched), "): ",
                       paste(unmatched, collapse = ", ")),
                type = "warning", duration = 8)
        }
    })

    sel <- reactive({
        req(input$genes)
        if (input$homeolog_mode == "merge") {
            bases <- unique(homeolog_base(input$genes))
            abund_long %>%
                mutate(gene = homeolog_base(gene)) %>%
                filter(gene %in% bases) %>%
                group_by(gene, rep, stage) %>%
                summarise(log2_abundance = log2_sum(log2_abundance),
                          .groups = "drop")
        } else {
            abund_long %>% filter(gene %in% input$genes)
        }
    })

    # RNA-seq selection reactives, homeolog-aware (sum TPM across .L + .S).
    rna_dev_sel <- reactive({
        req(input$genes)
        if (input$homeolog_mode == "merge") {
            bases <- unique(homeolog_base(input$genes))
            rna_dev_long %>%
                mutate(Geneid = homeolog_base(Geneid)) %>%
                filter(Geneid %in% bases) %>%
                group_by(Geneid, stage) %>%
                summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop")
        } else {
            rna_dev_long %>% filter(Geneid %in% input$genes)
        }
    })
    rna_tissue_sel <- reactive({
        req(input$genes)
        if (input$homeolog_mode == "merge") {
            bases <- unique(homeolog_base(input$genes))
            rna_tissue_long %>%
                mutate(Geneid = homeolog_base(Geneid)) %>%
                filter(Geneid %in% bases) %>%
                group_by(Geneid, tissue) %>%
                summarise(tpm = sum(tpm, na.rm = TRUE), .groups = "drop")
        } else {
            rna_tissue_long %>% filter(Geneid %in% input$genes)
        }
    })

    output$plot_protein <- renderPlotly({
        df <- sel()
        validate(need(nrow(df) > 0, "No data for the selected gene(s)."))

        if (input$rep_mode == "mean") {
            df2 <- df %>%
                group_by(gene, stage) %>%
                summarise(mean_abund = mean(log2_abundance, na.rm = TRUE),
                          lo = min(log2_abundance, na.rm = TRUE),
                          hi = max(log2_abundance, na.rm = TRUE),
                          .groups = "drop")
            p <- ggplot(df2, aes(x = stage, y = mean_abund, group = gene, color = gene,
                                 text = paste0("gene: ", gene,
                                               "<br>stage: ", stage,
                                               "<br>mean log2: ", round(mean_abund, 2),
                                               "<br>range: ", round(lo, 2), " - ", round(hi, 2)))) +
                geom_ribbon(aes(ymin = lo, ymax = hi, fill = gene),
                            alpha = 0.15, color = NA) +
                geom_line(linewidth = ln_width()) + geom_point(size = pt_size(), shape = pt_shape())
        } else {
            if (input$rep_mode %in% c("RepA", "RepB")) {
                df <- df %>% filter(rep == input$rep_mode)
            }
            p <- ggplot(df, aes(x = stage, y = log2_abundance,
                                group = interaction(gene, rep),
                                color = gene, linetype = rep,
                                text = paste0("gene: ", gene,
                                              "<br>rep: ", rep,
                                              "<br>stage: ", stage,
                                              "<br>log2 abundance: ", round(log2_abundance, 2)))) +
                geom_line(linewidth = ln_width()) + geom_point(size = pt_size(), shape = pt_shape())
        }
        p <- p + theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = "log2 abundance (MSstatsTMT-normalized)",
                 color = NULL, linetype = NULL)
        ggplotly(apply_palette(p), tooltip = "text")
    })

    output$plot_protein_hm <- renderPlotly({
        df <- sel()
        validate(need(nrow(df) > 0, "No data for the selected gene(s)."))
        mat <- df %>%
            group_by(gene, stage) %>%
            summarise(mean_abund = mean(log2_abundance, na.rm = TRUE), .groups = "drop") %>%
            group_by(gene) %>%
            mutate(z = if (sd(mean_abund, na.rm = TRUE) > 0)
                       (mean_abund - mean(mean_abund, na.rm = TRUE)) /
                           sd(mean_abund, na.rm = TRUE)
                   else 0) %>%
            ungroup()
        p <- ggplot(mat, aes(x = stage, y = gene, fill = z,
                             text = paste0(gene,
                                           "<br>stage: ", stage,
                                           "<br>z-score: ", round(z, 2),
                                           "<br>mean log2 abundance: ", round(mean_abund, 2)))) +
            geom_tile(color = "white") +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL)
        ggplotly(apply_heatmap_palette(p), tooltip = "text")
    })

    # ---- Plot styling helpers (driven by sidebar customization inputs) ----
    pt_shape <- reactive(as.numeric(input$point_shape))
    pt_size  <- reactive(input$point_size)
    ln_width <- reactive(input$line_width)

    apply_palette <- function(p) {
        pal <- input$palette
        if (is.null(pal) || pal == "default") return(p)
        if (pal %in% c("Set1", "Set2", "Dark2", "Paired")) {
            return(p + scale_color_brewer(palette = pal, na.value = "grey50") +
                   scale_fill_brewer(palette = pal, na.value = "grey50"))
        }
        opt <- switch(pal, "viridis" = "viridis", "plasma" = "plasma")
        p + scale_color_viridis_d(option = opt) +
            scale_fill_viridis_d(option = opt)
    }

    # Z-scored heatmaps use the defaults (legend "z-score", diverging about 0).
    # Raw-value heatmaps pass their own legend label and a midpoint centred on the
    # data range so the diverging palette keeps its contrast (raw values sit
    # entirely on one side of 0). Helper centres on the value range:
    mid_of <- function(x) mean(range(x, na.rm = TRUE))
    apply_heatmap_palette <- function(p, name = "z-score", midpoint = 0) {
        pal <- input$heatmap_palette
        if (is.null(pal)) pal <- "RdBu"
        diverging_endpoints <- list(
            RdBu  = c(low = "#2166AC", mid = "white", high = "#B2182B"),
            RdYlBu = c(low = "#4575B4", mid = "#FFFFBF", high = "#D73027"),
            PuOr  = c(low = "#542788", mid = "white", high = "#B35806")
        )
        if (pal %in% names(diverging_endpoints)) {
            ep <- diverging_endpoints[[pal]]
            return(p + scale_fill_gradient2(low = ep["low"], mid = ep["mid"],
                                            high = ep["high"], midpoint = midpoint,
                                            name = name))
        }
        opt <- switch(pal,
                      "viridis" = "viridis",
                      "magma" = "magma",
                      "plasma" = "plasma",
                      "cividis" = "cividis",
                      "viridis")
        p + scale_fill_viridis_c(option = opt, name = name)
    }

    zscore_by_group <- function(df, group_col, value_col) {
        df %>%
            group_by(.data[[group_col]]) %>%
            mutate(z = if (sd(.data[[value_col]], na.rm = TRUE) > 0)
                       (.data[[value_col]] - mean(.data[[value_col]], na.rm = TRUE)) /
                           sd(.data[[value_col]], na.rm = TRUE)
                   else 0) %>%
            ungroup()
    }

    output$plot_protein_hm_raw <- renderPlotly({
        df <- sel()
        validate(need(nrow(df) > 0, "No data for the selected gene(s)."))
        mat <- df %>%
            group_by(gene, stage) %>%
            summarise(mean_abund = mean(log2_abundance, na.rm = TRUE), .groups = "drop")
        p <- ggplot(mat, aes(x = stage, y = gene, fill = mean_abund,
                             text = paste0(gene,
                                           "<br>stage: ", stage,
                                           "<br>log2 abundance: ", round(mean_abund, 2)))) +
            geom_tile(color = "white") +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL)
        ggplotly(apply_heatmap_palette(p, name = "log2 abundance",
                                       midpoint = mid_of(mat$mean_abund)),
                 tooltip = "text")
    })

    output$plot_all_phos <- renderPlotly({
        req(input$genes)
        df <- all_phos_long %>%
            filter(gene %in% input$genes,
                   confidence_tier %in% input$tiers,
                   num_psms >= input$min_psms) %>%
            mutate(loc_display = if_else(
                       is.na(ascore), "not assessed",
                       paste0(loc_tier, " (AScore ", round(ascore, 1), ", ",
                              n_candidate_sites, " candidate S/T/Y in peptide)")),
                   motif_display = if_else(is.na(peptide_motif), "n/a", peptide_motif))
        validate(need(nrow(df) > 0,
                      "No phosphosites for the selected genes + filters. ",
                      "Try lowering the PSM cutoff or including the 'low' tier."))
        p <- ggplot(df, aes(x = stage, y = intensity, group = label, color = label,
                            text = paste0(label,
                                          "<br>tier: ", confidence_tier,
                                          "<br>num_psms: ", num_psms,
                                          "<br>localization: ", loc_display,
                                          "<br>peptide motif: ", motif_display,
                                          "<br>stage: ", stage,
                                          "<br>intensity: ",
                                          formatC(intensity, format = "g", digits = 4)))) +
            geom_line(alpha = 0.8, linewidth = ln_width()) +
            geom_point(size = pt_size(), shape = pt_shape()) +
            scale_y_log10(labels = scales::label_number(big.mark = ",")) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = "TMT intensity", color = NULL)
        ggplotly(apply_palette(p), tooltip = "text")
    })

    output$table_all_phos <- renderDT({
        req(input$genes)
        all_phos_raw %>%
            filter(gene %in% input$genes,
                   confidence_tier %in% input$tiers,
                   num_psms >= input$min_psms) %>%
            mutate(localization = coalesce(loc_tier, "not assessed"),
                   peptide_motif = coalesce(peptide_motif, "n/a")) %>%
            select(gene, site_position, residue, peptide_motif,
                   n_candidate_sites, ascore, localization,
                   confidence_tier, num_psms, num_peptides, best_delta_score,
                   all_of(INT_COLS), protein) %>%
            datatable(
                extensions = "Buttons",
                options = list(pageLength = 25, scrollX = TRUE,
                               dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
                rownames = FALSE,
                colnames = c("AScore" = "ascore",
                             "candidate S/T/Y" = "n_candidate_sites",
                             "motif (+/-7)" = "peptide_motif")
            ) %>%
            formatRound(INT_COLS, 1) %>%
            formatRound("best_delta_score", 3)
    })

    output$plot_rna_dev <- renderPlotly({
        req(input$genes)
        df <- rna_dev_sel()
        validate(need(nrow(df) > 0,
                      "No Session 2016 developmental RNA-seq for the selected gene(s). ",
                      "(Nomenclature mismatch? krt8.L (protein) -> krt8.1.L (RNA), etc.)"))
        p <- ggplot(df, aes(x = stage, y = tpm + 1, group = Geneid, color = Geneid,
                            text = paste0("gene: ", Geneid,
                                          "<br>stage: ", stage,
                                          "<br>TPM: ", round(tpm, 2)))) +
            geom_line(linewidth = ln_width()) + geom_point(size = pt_size(), shape = pt_shape()) +
            scale_y_log10(labels = scales::label_number(big.mark = ",")) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = "TPM + 1 (log10)", color = NULL)
        ggplotly(apply_palette(p), tooltip = "text")
    })

    output$plot_rna_dev_hm <- renderPlotly({
        req(input$genes)
        df <- rna_dev_sel() %>%
            zscore_by_group("Geneid", "tpm")
        validate(need(nrow(df) > 0, "No Session 2016 dev RNA-seq for selection."))
        p <- ggplot(df, aes(x = stage, y = Geneid, fill = z,
                            text = paste0(Geneid,
                                          "<br>stage: ", stage,
                                          "<br>z-score: ", round(z, 2),
                                          "<br>TPM: ", round(tpm, 2)))) +
            geom_tile(color = "white") +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL)
        ggplotly(apply_heatmap_palette(p), tooltip = "text")
    })

    output$plot_rna_dev_hm_raw <- renderPlotly({
        req(input$genes)
        df <- rna_dev_sel() %>%
            mutate(log10_tpm = log10(tpm + 1))
        validate(need(nrow(df) > 0, "No Session 2016 dev RNA-seq for selection."))
        p <- ggplot(df, aes(x = stage, y = Geneid, fill = log10_tpm,
                            text = paste0(Geneid,
                                          "<br>stage: ", stage,
                                          "<br>log10(TPM+1): ", round(log10_tpm, 2),
                                          "<br>TPM: ", round(tpm, 2)))) +
            geom_tile(color = "white") +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL)
        ggplotly(apply_heatmap_palette(p, name = "log10(TPM+1)",
                                       midpoint = mid_of(df$log10_tpm)),
                 tooltip = "text")
    })

    output$plot_rna_tissue <- renderPlotly({
        req(input$genes)
        df <- rna_tissue_sel()
        validate(need(nrow(df) > 0,
                      "No Session 2016 adult-tissue RNA-seq for the selected gene(s)."))
        p <- ggplot(df, aes(x = tissue, y = tpm, fill = Geneid,
                            text = paste0("gene: ", Geneid,
                                          "<br>tissue: ", tissue,
                                          "<br>TPM: ", round(tpm, 2)))) +
            geom_col(position = position_dodge(width = 0.85)) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = "TPM", fill = NULL)
        ggplotly(apply_palette(p), tooltip = "text")
    })

    output$plot_rna_tissue_hm <- renderPlotly({
        req(input$genes)
        df <- rna_tissue_sel() %>%
            zscore_by_group("Geneid", "tpm")
        validate(need(nrow(df) > 0, "No Session 2016 tissue RNA-seq for selection."))
        p <- ggplot(df, aes(x = tissue, y = Geneid, fill = z,
                            text = paste0(Geneid,
                                          "<br>tissue: ", tissue,
                                          "<br>z-score: ", round(z, 2),
                                          "<br>TPM: ", round(tpm, 2)))) +
            geom_tile(color = "white") +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL)
        ggplotly(apply_heatmap_palette(p), tooltip = "text")
    })

    output$plot_rna_tissue_hm_raw <- renderPlotly({
        req(input$genes)
        df <- rna_tissue_sel() %>%
            mutate(log10_tpm = log10(tpm + 1))
        validate(need(nrow(df) > 0, "No Session 2016 tissue RNA-seq for selection."))
        p <- ggplot(df, aes(x = tissue, y = Geneid, fill = log10_tpm,
                            text = paste0(Geneid,
                                          "<br>tissue: ", tissue,
                                          "<br>log10(TPM+1): ", round(log10_tpm, 2),
                                          "<br>TPM: ", round(tpm, 2)))) +
            geom_tile(color = "white") +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL)
        ggplotly(apply_heatmap_palette(p, name = "log10(TPM+1)",
                                       midpoint = mid_of(df$log10_tpm)),
                 tooltip = "text")
    })

    # ---- Combined matched-stages view ----
    combined_data <- reactive({
        req(input$genes)
        # Protein: average over RepA/RepB, restrict to matched proteomics stages.
        # sel() / rna_dev_sel() already apply the homeolog merge toggle.
        prot <- sel() %>%
            filter(as.character(stage) %in% MATCHED_STAGES$prot_stage) %>%
            group_by(gene, stage) %>%
            summarise(value = mean(log2_abundance, na.rm = TRUE), .groups = "drop") %>%
            rename(Geneid = gene) %>%
            left_join(MATCHED_STAGES %>%
                          mutate(stage = factor(prot_stage, levels = STAGES)),
                      by = "stage") %>%
            mutate(modality = "Protein") %>%
            select(Geneid, pair_label, value, modality)

        # RNA: log10(TPM+1), restrict to matched RNA stages
        rna <- rna_dev_sel() %>%
            filter(as.character(stage) %in% MATCHED_STAGES$rna_stage) %>%
            mutate(value = log10(tpm + 1)) %>%
            left_join(MATCHED_STAGES %>%
                          mutate(stage = factor(rna_stage, levels = DEV_STAGES_RNA)),
                      by = "stage") %>%
            mutate(modality = "RNA") %>%
            select(Geneid, pair_label, value, modality)

        both <- bind_rows(prot, rna) %>%
            mutate(pair_label = factor(pair_label, levels = MATCHED_STAGES$pair_label),
                   modality = factor(modality, levels = c("Protein", "RNA")),
                   series = paste0(Geneid, " — ", modality)) %>%
            group_by(Geneid, modality) %>%
            mutate(z = if (sd(value, na.rm = TRUE) > 0)
                       (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
                   else 0) %>%
            ungroup()
        both
    })

    # Facet labeller so heatmap headers show units even though data uses
    # short modality names ("Protein", "RNA") to keep the line-plot legend clean.
    modality_labeller <- as_labeller(c(
        "Protein" = "Protein (log2 abundance)",
        "RNA"     = "RNA (log10 TPM+1)"
    ))

    output$plot_combined <- renderPlotly({
        df <- combined_data()
        validate(need(nrow(df) > 0,
                      "No matched-stage data for selection. ",
                      "Check that the gene exists in both modalities."))
        # One legend entry per (gene, modality) pair, with explicit label.
        # Shape encodes modality (circle = Protein, triangle = RNA);
        # linetype reinforces visually (solid vs longdash).
        p <- ggplot(df, aes(x = pair_label, y = z,
                            group = series, color = series,
                            linetype = modality, shape = modality,
                            text = paste0("gene: ", Geneid,
                                          "<br>modality: ", modality,
                                          "<br>stage: ", pair_label,
                                          "<br>z-score: ", round(z, 2),
                                          "<br>raw value: ", round(value, 2)))) +
            geom_line(linewidth = ln_width()) +
            geom_point(size = pt_size()) +
            scale_linetype_manual(values = c(Protein = "solid", RNA = "longdash"),
                                  guide = "none") +
            scale_shape_manual(values = c(Protein = 16, RNA = 17),
                               guide = "none") +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = "Z-score (within gene, within modality)",
                 color = "Gene — modality")
        ggplotly(apply_palette(p), tooltip = "text")
    })

    output$plot_combined_hm <- renderPlotly({
        df <- combined_data()
        validate(need(nrow(df) > 0, "No matched-stage data for selection."))
        # Faceted heatmap: protein left, RNA right
        p <- ggplot(df, aes(x = pair_label, y = Geneid, fill = z,
                            text = paste0(Geneid,
                                          "<br>modality: ", modality,
                                          "<br>stage: ", pair_label,
                                          "<br>z-score: ", round(z, 2),
                                          "<br>raw value: ", round(value, 2)))) +
            geom_tile(color = "white") +
            facet_wrap(~ modality, nrow = 1, labeller = modality_labeller) +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL)
        ggplotly(apply_heatmap_palette(p), tooltip = "text")
    })

    output$plot_combined_hm_raw <- renderPlotly({
        df <- combined_data()
        validate(need(nrow(df) > 0, "No matched-stage data for selection."))
        # Each modality on its own color scale (different units)
        p <- ggplot(df, aes(x = pair_label, y = Geneid, fill = value,
                            text = paste0(Geneid,
                                          "<br>modality: ", modality,
                                          "<br>stage: ", pair_label,
                                          "<br>value: ", round(value, 2)))) +
            geom_tile(color = "white") +
            facet_wrap(~ modality, nrow = 1, scales = "free",
                       labeller = modality_labeller) +
            theme_minimal(base_size = 11) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL)
        ggplotly(apply_heatmap_palette(p, name = "raw value",
                                       midpoint = mid_of(df$value)),
                 tooltip = "text")
    })

    output$table_protein <- renderDT({
        req(input$genes)
        abund_raw %>%
            filter(gene %in% input$genes) %>%
            select(gene, Index, NumberPSM, MaxPepProb, ReferenceIntensity,
                   matches("^Rep[AB]_")) %>%
            datatable(
                extensions = "Buttons",
                options = list(pageLength = 25, scrollX = TRUE,
                               dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
                rownames = FALSE
            ) %>%
            formatRound(c("ReferenceIntensity",
                          grep("^Rep[AB]_", names(abund_raw), value = TRUE)),
                        digits = 2)
    })

    output$table_rna_dev <- renderDT({
        req(input$genes)
        df <- rna_dev %>% filter(Geneid %in% input$genes)
        validate(need(nrow(df) > 0, "No Session 2016 dev RNA-seq for selection."))
        datatable(
            df,
            extensions = "Buttons",
            options = list(pageLength = 25, scrollX = TRUE,
                           dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
            rownames = FALSE
        ) %>%
            formatRound(DEV_STAGES_RNA, digits = 2)
    })

    output$table_rna_tissue <- renderDT({
        req(input$genes)
        df <- rna_tissue %>% filter(Geneid %in% input$genes)
        validate(need(nrow(df) > 0, "No Session 2016 tissue RNA-seq for selection."))
        datatable(
            df,
            extensions = "Buttons",
            options = list(pageLength = 25, scrollX = TRUE,
                           dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
            rownames = FALSE
        ) %>%
            formatRound(TISSUES_RNA, digits = 2)
    })

}

shinyApp(ui, server)
