# ============================================================================
# Cross-species protein alignment - sourced by app.R as a tab.
# Defines: alignment_tab_ui (a sidebarLayout) and
#          alignment_tab_server(input, output, session).
# Users pick X. laevis / human / mouse genes -> MSA (msa, no system binary) ->
# conservation + residue-equivalency views.  Libraries are loaded in app.R.
# ============================================================================

DATA_FILE <- if (file.exists("sequences.tsv")) "sequences.tsv" else "dev_sequences.tsv"
seqs <- readr::read_tsv(DATA_FILE, show_col_types = FALSE)

choices_for <- function(sp) {
    s <- seqs %>% filter(species == sp) %>% arrange(gene)
    setNames(s$accession, s$gene)   # value = accession (unique), label = gene
}
SP_SHORT <- c("X. laevis" = "Xl", "human" = "Hs", "mouse" = "Mm")

# Amino-acid chemistry categories (Clustal-style) -> one colour per category,
# so the alignment legend is compact and meaningful.
AA_CATEGORY <- c(
    A="Hydrophobic", I="Hydrophobic", L="Hydrophobic", M="Hydrophobic",
    F="Hydrophobic", W="Hydrophobic", V="Hydrophobic", C="Hydrophobic",
    K="Positive", R="Positive",
    D="Negative", E="Negative",
    N="Polar", Q="Polar", S="Polar", T="Polar",
    G="Glycine", P="Proline",
    H="Aromatic", Y="Aromatic", "-"="Gap")
CAT_LEVELS <- c("Hydrophobic","Positive","Negative","Polar",
                "Aromatic","Glycine","Proline","Gap")
CAT_COLORS <- c(Hydrophobic="#80a0f0", Positive="#f01505", Negative="#c048c0",
                Polar="#15c015", Aromatic="#15a4a4", Glycine="#f09048",
                Proline="#c0c000", Gap="#eeeeee")

# Alternative per-residue scheme (Taylor 1997, as used by Jalview).
TAYLOR <- c(A="#ccff00", R="#0000ff", N="#cc00ff", D="#ff0000", C="#ffff00",
            Q="#ff00cc", E="#ff0066", G="#ff9900", H="#0066ff", I="#66ff00",
            L="#33ff00", K="#6600ff", M="#00ff00", F="#00ff66", P="#ffcc00",
            S="#ff3300", T="#ff6600", W="#00ccff", Y="#00ffcc", V="#99ff00",
            "-"="#eeeeee")

# Friendly colour choices (R named colours) + ready-made category palettes,
# so users pick from a menu instead of typing hex.
COLOR_CHOICES <- c(
    "royalblue","steelblue","dodgerblue","skyblue","navy","blue",
    "red","firebrick","red3","salmon",
    "darkorange","orange","gold","yellow",
    "limegreen","green3","forestgreen","seagreen","olivedrab",
    "turquoise","cyan3",
    "purple","mediumorchid","magenta","violet","deeppink",
    "brown","tan","pink",
    "grey30","grey40","grey45","grey50","grey60","grey70","grey80","grey90",
    "black","white")
DEFAULT_NAMED <- c(Hydrophobic = "royalblue", Positive = "red", Negative = "purple",
                   Polar = "limegreen", Aromatic = "turquoise", Glycine = "darkorange",
                   Proline = "gold", Gap = "grey80")
PALETTES <- list(
    "Default"          = DEFAULT_NAMED,
    "Colourblind-safe" = c(Hydrophobic="skyblue", Positive="darkorange",
                           Negative="mediumorchid", Polar="seagreen",
                           Aromatic="dodgerblue", Glycine="orange",
                           Proline="gold", Gap="grey80"),
    "Vivid"            = c(Hydrophobic="blue", Positive="red", Negative="magenta",
                           Polar="green3", Aromatic="cyan3", Glycine="darkorange",
                           Proline="yellow", Gap="grey80"),
    "Grayscale"        = c(Hydrophobic="grey30", Positive="grey50", Negative="grey40",
                           Polar="grey60", Aromatic="grey45", Glycine="grey70",
                           Proline="grey80", Gap="grey90"))

# A colour is valid if 6-digit hex (optionally '#') OR an R colour name.
valid_color <- function(x) grepl("^#?[0-9A-Fa-f]{6}$", x) ||
    tolower(x) %in% tolower(grDevices::colors())
norm_color  <- function(x) if (grepl("^[0-9A-Fa-f]{6}$", x)) paste0("#", x) else x

# Parse "KEY=colour" pairs (comma/newline-separated; colour = hex or name).
parse_kv <- function(txt) {
    if (is.null(txt) || !nzchar(trimws(txt))) return(character(0))
    toks <- trimws(unlist(strsplit(txt, "[,\n]+"))); toks <- toks[nzchar(toks)]
    out <- character(0)
    for (t in toks) {
        kv <- trimws(strsplit(t, "=")[[1]])
        if (length(kv) == 2 && valid_color(kv[2])) out[kv[1]] <- norm_color(kv[2])
    }
    out
}

wrap60 <- function(x) gsub("(.{60})", "\\1\n", x)

# Clustal consensus groups: a column scores ':' if all its residues fall in one
# "strong" group, '.' if in one "weak" group, '*' if identical, else blank.
STRONG <- strsplit(c("STA","NEQK","NHQK","NDEQ","QHRK","MILV","MILF","HY","FYW"), "")
WEAK   <- strsplit(c("CSA","ATV","SAG","STNK","STPA","SGND","SNDEQK","NDEQHK",
                     "NEQHRK","FVLIM","HFY"), "")
consensus_symbol <- function(col) {
    if (any(col == "-")) return(" ")
    u <- unique(col)
    if (length(u) == 1) return("*")
    if (any(vapply(STRONG, function(g) all(u %in% g), logical(1)))) return(":")
    if (any(vapply(WEAK,   function(g) all(u %in% g), logical(1)))) return(".")
    " "
}

alignment_tab_ui <- sidebarLayout(
        sidebarPanel(
            width = 3,
            helpText("Pick >=2 genes across the three species (each box is ",
                     "multi-select), or paste a list below, then Align."),
            selectizeInput("xl", "X. laevis gene(s):", choices = NULL, multiple = TRUE),
            selectizeInput("hs", "Human gene(s):",     choices = NULL, multiple = TRUE),
            selectizeInput("mm", "Mouse gene(s):",      choices = NULL, multiple = TRUE),
            textAreaInput("gene_paste", "Paste gene list (any species):",
                          placeholder = "krt8.L, KRT18, dsp.S  (comma / space / newline)",
                          rows = 2),
            actionButton("paste_add", "Add pasted genes", class = "btn-sm"),
            br(), br(),
            selectInput("method", "Alignment method:",
                        c("ClustalOmega", "ClustalW", "Muscle")),
            checkboxInput("show_letters", "Show residue letters on alignment", TRUE),
            actionButton("align", "Align", class = "btn-primary"),
            br(), br(),
            verbatimTextOutput("eta"),
            helpText(em("Seed set = keratins only (93 seqs). Swap in full ",
                        "proteomes later.")),
            hr(),
            h5("Colouring"),
            selectInput("color_scheme", "Scheme:",
                        c("Chemistry (categories)" = "chemistry",
                          "Taylor (per-residue)" = "taylor")),
            selectInput("palette_preset", "Category palette:",
                        choices = names(PALETTES), selected = "Default"),
            checkboxInput("custom_cats", "Fine-tune individual category colours", FALSE),
            conditionalPanel(
                "input.custom_cats == true",
                do.call(tagList, lapply(CAT_LEVELS, function(ct)
                    selectInput(paste0("cat_", ct), paste0(ct, ":"),
                                choices = COLOR_CHOICES, selected = DEFAULT_NAMED[[ct]])))
            ),
            helpText(em("Pick a category palette (Chemistry scheme) or switch to the ",
                        "Taylor per-residue scheme; tick fine-tune to set individual ",
                        "category colours.")),
            hr(),
            h5("Selected sequences (FASTA)"),
            tags$div(style = "max-height: 320px; overflow-y: auto; font-size: 11px;",
                     verbatimTextOutput("seq_fasta"))
        ),
        mainPanel(
            width = 9,
            verbatimTextOutput("summary_line"),
            h4("Pairwise % identity (all vs all)"),
            helpText(em("Identical residues / aligned non-gap positions, per pair.")),
            DTOutput("pid_tbl"),
            h4("Conservation across alignment"),
            helpText(em("Drag left-right here OR on the alignment below to select a "),
                     em("region of interest; its conservation is summarised below. "),
                     em("Use this bar for regions spanning multiple wrapped blocks. "),
                     em("Double-click here (or use Reset selection) to clear.")),
            plotOutput("cons_plot", height = "200px",
                       brush = brushOpts(id = "cons_brush", direction = "x"),
                       dblclick = "cons_dblclick"),
            downloadButton("dl_cons_plot", "Download figure", class = "btn-sm"),
            actionButton("clear_region", "Reset selection", class = "btn-warning btn-sm"),
            verbatimTextOutput("region_stats"),
            h5("Residue numbering across species (selected region)"),
            helpText(em("Native (per-sequence) residue numbers for the selected "),
                     em("columns. Select a single column - click the same residue "),
                     em("twice - to compare one residue across species.")),
            DTOutput("residue_map"),
            h4("Alignment"),
            helpText(em("Wrapped at 60 columns; sequence tiles coloured by AA ",
                        "chemistry. Per column, top to bottom: a ", strong("Consensus"),
                        " symbol (", strong("*"), " identical, ", strong(":"),
                        " strongly similar chemistry, ", strong("."), " weakly similar) ",
                        "and a ", strong("Conservation"), " square coloured by % conserved. ",
                        "Block headers show the block mean % conserved. ",
                        strong("Drag within a block, or click a start then an end "),
                        strong("residue (works across blocks), to select a region "),
                        strong("(outlined in orange)."))),
            plotOutput("aln_plot", height = "auto",
                       brush = brushOpts(id = "aln_brush", direction = "x"),
                       click = clickOpts(id = "aln_click")),
            downloadButton("dl_aln_plot", "Download figure", class = "btn-sm"),
            actionButton("clear_region_btm", "Reset selection", class = "btn-warning btn-sm")
        )
    )

alignment_tab_server <- function(input, output, session, figs = reactiveValues()) {
    updateSelectizeInput(session, "xl", choices = choices_for("X. laevis"),
                         selected = choices_for("X. laevis")[["krt8.L"]], server = TRUE)
    updateSelectizeInput(session, "hs", choices = choices_for("human"),
                         selected = choices_for("human")[["KRT8"]], server = TRUE)
    updateSelectizeInput(session, "mm", choices = choices_for("mouse"),
                         selected = choices_for("mouse")[["Krt8"]], server = TRUE)

    # Paste a gene list -> distribute matches to the right species picker(s).
    SP_BY_ID <- c(xl = "X. laevis", hs = "human", mm = "mouse")
    # Normalise a symbol for cross-species matching: lower-case, strip the .L/.S
    # homeolog suffix and the Xenopus .N sub-model number, so KRT9 == krt9.1.L.
    norm_gene <- function(x) sub("\\.[0-9]+$", "", sub("\\.[ls]$", "", tolower(x)))
    # Per-species lookup tables (built once).
    SP_TAB <- lapply(SP_BY_ID, function(sp) {
        d <- seqs[seqs$species == sp, c("gene", "accession")]
        d$lg <- tolower(d$gene); d$core <- norm_gene(d$gene); d
    })

    observeEvent(input$paste_add, {
        req(input$gene_paste)
        toks <- trimws(unlist(strsplit(input$gene_paste, "[,;\n[:space:]]+")))
        toks <- toks[nzchar(toks)]
        req(length(toks) > 0)

        add <- list(xl = character(0), hs = character(0), mm = character(0))
        report <- character(0); none <- character(0)
        for (tok in toks) {
            tl <- tolower(tok); tc <- norm_gene(tok); parts <- character(0)
            for (id in names(SP_BY_ID)) {
                d <- SP_TAB[[id]]; tag <- SP_SHORT[[SP_BY_ID[[id]]]]
                ex <- d$lg == tl
                if (any(ex)) {                              # exact (case-insensitive)
                    add[[id]] <- c(add[[id]], d$accession[ex])
                    parts <- c(parts, paste0(tag, " ", paste(d$gene[ex], collapse = "/")))
                } else {
                    nm <- d$core == tc                      # cross-species normalised
                    if (any(nm)) {
                        add[[id]] <- c(add[[id]], d$accession[nm])
                        parts <- c(parts, paste0(tag, " ", paste(d$gene[nm], collapse = "/"), "*"))
                    }
                }
            }
            if (length(parts)) report <- c(report, paste0("<b>", tok, "</b> -> ",
                                                          paste(parts, collapse = ", ")))
            else none <- c(none, tok)
        }
        for (id in names(SP_BY_ID))
            if (length(add[[id]]))
                updateSelectizeInput(session, id, choices = choices_for(SP_BY_ID[[id]]),
                                     selected = union(input[[id]], unique(add[[id]])),
                                     server = TRUE)

        # Suggestions for tokens that matched nothing at all.
        if (length(none)) {
            allg <- unique(seqs$gene); gl <- tolower(allg)
            for (tok in head(none, 20)) {
                tl <- tolower(tok)
                cand <- allg[norm_gene(allg) == norm_gene(tok)]
                if (!length(cand)) cand <- allg[startsWith(gl, tl)]
                if (!length(cand)) cand <- allg[adist(tl, gl)[1, ] <= 1]
                report <- c(report, if (length(cand))
                    paste0("<b>", tok, "</b> - not found; did you mean: ",
                           paste(head(unique(cand), 6), collapse = ", "))
                    else paste0("<b>", tok, "</b> - not found"))
            }
        }
        if (length(report))
            showNotification(
                HTML(paste0(paste(report, collapse = "<br>"),
                            "<br><i>* added by cross-species name match</i>")),
                type = if (length(none)) "warning" else "message", duration = 14)
    })

    # Live: whatever is currently selected, independent of the Align button.
    selected_seqs <- reactive({
        accs <- c(input$xl, input$hs, input$mm)
        seqs %>% filter(accession %in% accs) %>%
            mutate(label = paste0(SP_SHORT[species], " ", gene))
    })

    output$seq_fasta <- renderText({
        s <- selected_seqs()
        if (!nrow(s)) return("(no genes selected)")
        paste(sprintf(">%s | %s | %d aa\n%s",
                      s$label, s$accession, s$length, wrap60(s$sequence)),
              collapse = "\n\n")
    })

    # Live rough ETA from the current selection, shown before aligning.
    output$eta <- renderText({
        s <- selected_seqs(); n <- nrow(s)
        if (n < 2) return("Select >=2 genes, then click Align.")
        total_aa <- sum(s$length)
        est <- ceiling(1 + 0.25 * n + total_aa / 4000 + 0.05 * n * (n - 1))
        sprintf("%d sequences, %d aa total.\nRough estimate: ~%d s to align + render.",
                n, total_aa, est)
    })

    # Gather selected sequences into a named AAStringSet (only on Align click).
    picked <- eventReactive(input$align, {
        s <- selected_seqs()
        validate(need(nrow(s) >= 2, "Select at least two genes, then Align."))
        aa <- AAStringSet(s$sequence); names(aa) <- s$label
        aa
    })

    align_secs <- reactiveVal(NA_real_)
    aligned <- reactive({
        aa <- picked()
        t0 <- Sys.time()
        m <- withProgress(message = "Aligning sequences...", value = 0.2, {
            aln <- msa(aa, method = input$method)
            incProgress(0.6, message = "Building figures...")
            al <- as(aln, "AAStringSet")
            mm <- do.call(rbind, strsplit(as.character(al), ""))
            rownames(mm) <- names(al)
            mm
        })
        align_secs(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
        m
    })

    output$summary_line <- renderText({
        m <- aligned()
        sprintf("%d sequences | alignment length %d | method %s | done in %s s",
                nrow(m), ncol(m), input$method, align_secs())
    })

    # All-vs-all pairwise % identity over the ungapped overlap of each pair.
    output$pid_tbl <- renderDT({
        m <- aligned(); n <- nrow(m)
        pid <- function(a, b) {
            keep <- a != "-" & b != "-"
            if (!any(keep)) return(NA_real_)
            round(100 * mean(a[keep] == b[keep]), 1)
        }
        M <- matrix(100, n, n, dimnames = list(rownames(m), rownames(m)))
        if (n >= 2) for (i in 1:(n - 1)) for (j in (i + 1):n) {
            v <- pid(m[i, ], m[j, ]); M[i, j] <- v; M[j, i] <- v
        }
        datatable(as.data.frame(M), options = list(dom = "t", scrollX = TRUE),
                  rownames = TRUE) %>%
            formatStyle(colnames(as.data.frame(M)),
                        backgroundColor = styleInterval(c(40, 70, 90),
                            c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5")))
    })

    cons <- reactive({
        m <- aligned()
        fr <- vapply(seq_len(ncol(m)), function(j) {
            col <- m[, j][m[, j] != "-"]
            if (length(col)) max(table(col)) / nrow(m) else 0
        }, numeric(1))
        consensus <- vapply(seq_len(ncol(m)), function(j) {
            tb <- sort(table(m[, j][m[, j] != "-"]), decreasing = TRUE)
            if (length(tb)) names(tb)[1] else "-"
        }, character(1))
        symbol <- vapply(seq_len(ncol(m)), function(j) consensus_symbol(m[, j]),
                         character(1))
        data.frame(pos = seq_len(ncol(m)), frac = fr,
                   consensus = consensus, symbol = symbol)
    })

    BLOCK <- 60   # alignment columns per wrapped row (used by plot + brushes)

    # Region of interest (alignment columns). NULL = whole alignment. Updated by
    # brushing EITHER the conservation bar OR the amino-acid alignment itself.
    sel_region <- reactiveVal(NULL)
    anchor <- reactiveVal(NULL)   # first click of a two-click alignment selection
    region <- reactive({
        r <- sel_region(); N <- nrow(cons())
        if (is.null(r)) c(1L, N) else c(max(1L, r[1]), min(N, r[2]))
    })
    observeEvent(input$cons_brush, {
        N <- nrow(cons())
        a <- max(1L, round(input$cons_brush$xmin)); b <- min(N, round(input$cons_brush$xmax))
        sel_region(sort(c(a, b))); anchor(NULL)
    })
    observeEvent(input$aln_brush, {
        # Alignment is faceted by block; panelvar1 is the block index, x is the
        # within-block column (1..BLOCK). Map both back to a global column range.
        bi <- suppressWarnings(as.integer(input$aln_brush$panelvar1))
        if (is.na(bi)) return()
        x0 <- max(1, floor(input$aln_brush$xmin)); x1 <- min(BLOCK, ceiling(input$aln_brush$xmax))
        sel_region(sort(c((bi - 1) * BLOCK + x0, (bi - 1) * BLOCK + x1))); anchor(NULL)
    })
    # Click two residues (even in different blocks) to select a region that wraps.
    observeEvent(input$aln_click, {
        bi <- suppressWarnings(as.integer(input$aln_click$panelvar1))
        if (is.na(bi)) return()
        N <- nrow(cons())
        pos <- min(max((bi - 1) * BLOCK + round(input$aln_click$x), 1L), N)
        if (is.null(anchor())) anchor(pos)              # first click: set start
        else { sel_region(sort(c(anchor(), pos))); anchor(NULL) }  # second: complete
    })
    reset_sel <- function() {
        sel_region(NULL); anchor(NULL)
        session$resetBrush("cons_brush"); session$resetBrush("aln_brush")
    }
    observeEvent(input$clear_region, reset_sel())       # button above alignment
    observeEvent(input$clear_region_btm, reset_sel())   # button below alignment
    observeEvent(input$cons_dblclick, reset_sel())      # double-click conservation bar
    observeEvent(input$align, reset_sel())              # new alignment -> clear stale region

    # Conservation / identity stats for the current region (shared by the stats
    # box and the alignment subtitle).
    region_summary <- reactive({
        cn <- cons(); m <- aligned(); rg <- region(); reg <- rg[1]:rg[2]
        sym <- cn$symbol[reg]; sub <- m[, reg, drop = FALSE]; nseq <- nrow(m)
        ids <- c()
        if (nseq >= 2) for (i in 1:(nseq - 1)) for (j in (i + 1):nseq) {
            keep <- sub[i, ] != "-" & sub[j, ] != "-"
            if (any(keep)) ids <- c(ids, mean(sub[i, keep] == sub[j, keep]))
        }
        list(a = rg[1], b = rg[2], n = length(reg),
             meancons = 100 * mean(cn$frac[reg]),
             nstar = sum(sym == "*"), ncolon = sum(sym == ":"), ndot = sum(sym == "."),
             pid = if (length(ids)) 100 * mean(ids) else NA_real_,
             scope = if (is.null(sel_region())) "whole alignment" else "selected region")
    })

    output$cons_plot <- renderPlot({
        d <- cons(); rg <- region()
        p <- ggplot(d, aes(pos, frac))
        if (!is.null(sel_region()))
            p <- p + annotate("rect", xmin = rg[1] - 0.5, xmax = rg[2] + 0.5,
                              ymin = 0, ymax = 1, fill = "#ffd54f", alpha = 0.4)
        g <- p + geom_col(width = 1, fill = "#2166AC") +
            scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
            scale_x_continuous(expand = c(0, 0)) +
            labs(x = "alignment position", y = "% conserved") +
            theme_minimal(base_size = 11)
        figs[["cons_plot"]] <- g
        figs[["__dim_cons_plot"]] <- c(w = 11, h = 2.6)   # wide, short track
        g
    })

    output$region_stats <- renderText({
        s <- region_summary()
        sprintf(paste0("Region (%s): positions %d-%d  (%d columns)\n",
                       "Mean %% conserved: %.0f%%\n",
                       "Identical (*): %d  |  strongly similar (:): %d  |  ",
                       "weakly similar (.): %d\n",
                       "Mean pairwise %% identity: %s"),
                s$scope, s$a, s$b, s$n, s$meancons,
                s$nstar, s$ncolon, s$ndot,
                if (is.na(s$pid)) "n/a" else sprintf("%.1f%%", s$pid))
    })

    # Map the selected alignment columns to each sequence's own residue numbering
    # (gaps excluded). For a single-column selection this gives the equivalent
    # residue in every species, e.g. Xl krt8.S S418 = Hs KRT8 S420.
    output$residue_map <- renderDT({
        m <- aligned(); rg <- region(); a <- rg[1]; b <- rg[2]
        rows <- lapply(seq_len(nrow(m)), function(i) {
            row <- m[i, ]; nong <- row != "-"; cumn <- cumsum(nong)
            idx <- which(nong[a:b]) + (a - 1)            # non-gap columns in region
            if (!length(idx))
                return(data.frame(sequence = rownames(m)[i],
                                  residue_range = "- (all gap here)", n = 0L))
            f <- idx[1]; l <- idx[length(idx)]
            rng <- if (f == l) sprintf("%s%d", row[f], cumn[f])
                   else sprintf("%s%d - %s%d", row[f], cumn[f], row[l], cumn[l])
            data.frame(sequence = rownames(m)[i], residue_range = rng,
                       n = length(idx))
        })
        datatable(do.call(rbind, rows), rownames = FALSE,
                  options = list(dom = "t", scrollX = TRUE),
                  colnames = c("sequence", "native residue range", "residues in region"))
    })

    # ---- Colouring (customisable scheme + per-residue overrides + sites) ----
    # Selecting a palette preset fills in the per-category dropdowns.
    observeEvent(input$palette_preset, {
        pal <- PALETTES[[input$palette_preset]]
        if (is.null(pal)) return()
        for (ct in CAT_LEVELS)
            updateSelectInput(session, paste0("cat_", ct), selected = pal[[ct]])
    }, ignoreInit = TRUE)
    cat_colors_map <- reactive({
        setNames(vapply(CAT_LEVELS, function(ct) {
            v <- input[[paste0("cat_", ct)]]
            if (is.null(v) || !nzchar(v)) DEFAULT_NAMED[[ct]] else v
        }, character(1)), CAT_LEVELS)
    })
    # Per-residue fill: a group label (for the legend) and a group->colour map.
    aa_style <- reactive({
        if (identical(input$color_scheme, "taylor"))
            list(grp = setNames(names(TAYLOR), names(TAYLOR)), col = TAYLOR)
        else
            list(grp = AA_CATEGORY, col = cat_colors_map())
    })

    # Static plot renders server-side in <1s (vs ~15s for interactive plotly on
    # a long alignment). Top to bottom per block: Consensus symbols, Conservation
    # squares, sequences. Brushable to select a region (see aln_brush observer).
    output$aln_plot <- renderPlot({
        m <- aligned(); cn <- cons()
        np <- ncol(m); nseq <- nrow(m)
        nblocks <- ceiling(np / BLOCK)
        blk <- function(pos) ((pos - 1) %/% BLOCK) + 1
        blkmean <- vapply(seq_len(nblocks),
                          function(b) mean(cn$frac[blk(cn$pos) == b]), numeric(1))
        blab <- sprintf("positions %d-%d  (%.0f%% conserved)",
                        (seq_len(nblocks) - 1) * BLOCK + 1,
                        pmin(seq_len(nblocks) * BLOCK, np), 100 * blkmean)
        # Facet by integer block id so the brush's panelvar1 maps straight back
        # to a block index; show the descriptive label via a labeller. When a
        # region is selected, append its summary to every block header so it is
        # visible without scrolling back to the top.
        s <- region_summary()
        strip_lab <- blab
        if (!is.null(sel_region()))
            strip_lab <- paste0(blab, "\nselection ", s$a, "-", s$b, ": ",
                                sprintf("%.0f%% conserved", s$meancons),
                                if (!is.na(s$pid)) sprintf(", %.0f%% identity", s$pid) else "")
        labeller_map <- as_labeller(setNames(strip_lab, as.character(seq_len(nblocks))))
        bfac <- function(pos) factor(blk(pos), levels = seq_len(nblocks))
        ylev <- c(rev(rownames(m)), "Conservation", "Consensus")

        # Per-block x-range of the selected region, for highlighting.
        hb <- NULL
        if (!is.null(sel_region())) {
            rg <- region()
            for (b in seq_len(nblocks)) {
                g0 <- (b - 1) * BLOCK + 1; g1 <- min(b * BLOCK, np)
                lo <- max(rg[1], g0); hi <- min(rg[2], g1)
                if (lo <= hi) hb <- rbind(hb, data.frame(
                    block = factor(b, levels = seq_len(nblocks)),
                    xmin = (lo - g0 + 1) - 0.5, xmax = (hi - g0 + 1) + 0.5))
            }
        }

        sdf <- data.frame(
            seq = factor(rep(rownames(m), times = np), levels = ylev),
            pos = rep(seq_len(np), each = nseq),
            residue = as.vector(m))
        sdf$xin <- ((sdf$pos - 1) %% BLOCK) + 1
        sdf$block <- bfac(sdf$pos)
        st <- aa_style()
        sdf$group <- factor(unname(st$grp[sdf$residue]), levels = names(st$col))

        cons_sq <- data.frame(seq = factor("Conservation", levels = ylev),
                              pos = cn$pos, frac = cn$frac)
        cons_sq$xin <- ((cons_sq$pos - 1) %% BLOCK) + 1
        cons_sq$block <- bfac(cons_sq$pos)

        sym <- data.frame(seq = factor("Consensus", levels = ylev),
                          pos = cn$pos, symbol = cn$symbol)
        sym$xin <- ((sym$pos - 1) %% BLOCK) + 1
        sym$block <- bfac(sym$pos)

        p <- ggplot(mapping = aes(xin, seq)) +
            # conservation squares (continuous fill scale)
            geom_tile(data = cons_sq, aes(fill = frac), colour = "white", linewidth = 0.1) +
            scale_fill_gradient(low = "grey88", high = "#08306b", limits = c(0, 1),
                                labels = scales::percent, name = "% conserved") +
            new_scale_fill() +
            # sequence tiles (discrete residue/category fill scale)
            geom_tile(data = sdf, aes(fill = group), colour = "white", linewidth = 0.1) +
            scale_fill_manual(values = st$col, drop = FALSE, na.value = "#dddddd",
                              name = if (identical(input$color_scheme, "taylor"))
                                         "residue" else "AA colour")
        # Translucent highlight band over the selected columns (full block height,
        # so consensus + conservation + sequences are all highlighted). Drawn
        # under the text so letters/symbols stay crisp.
        if (!is.null(hb))
            p <- p + geom_rect(data = hb, inherit.aes = FALSE,
                               aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                               fill = "#ffe082", alpha = 0.35, colour = NA)
        p <- p +
            # consensus notation symbols on top
            geom_text(data = sym, aes(label = symbol), size = 5, fontface = "bold",
                      colour = "grey20") +
            facet_wrap(~ block, ncol = 1, labeller = labeller_map) +
            scale_y_discrete(limits = ylev) +   # bottom -> top: seqs, Conservation, Consensus
            scale_x_continuous(breaks = seq(0, BLOCK, 10), expand = c(0, 0)) +
            labs(x = NULL, y = NULL,
                 subtitle = sprintf("%s: positions %d-%d  -  %.0f%% conserved  -  %s identity",
                             s$scope, s$a, s$b, s$meancons,
                             if (is.na(s$pid)) "n/a" else sprintf("%.0f%%", s$pid))) +
            theme_minimal(base_size = 11) +
            theme(panel.grid = element_blank(),
                  plot.subtitle = element_text(face = "bold", colour = "#e65100"))
        if (isTRUE(input$show_letters))
            p <- p + geom_text(data = sdf, aes(label = residue), size = 3, colour = "grey15")
        # Dashed marker at the pending first click (waiting for the second click).
        if (!is.null(anchor())) {
            ap <- anchor()
            adf <- data.frame(block = factor(((ap - 1) %/% BLOCK) + 1, levels = seq_len(nblocks)),
                              xin = ((ap - 1) %% BLOCK) + 1)
            p <- p + geom_vline(data = adf, aes(xintercept = xin),
                                colour = "#e65100", linetype = "dashed", linewidth = 0.6)
        }
        # Bold border around the selected region, on top of everything.
        if (!is.null(hb))
            p <- p + geom_rect(data = hb, inherit.aes = FALSE,
                               aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                               fill = NA, colour = "#e65100", linewidth = 1.1)
        figs[["aln_plot"]] <- p
        # Export size matched to the on-screen wrapped layout (px height / 96 dpi)
        # so the downloaded alignment isn't squished into the default 5 in.
        figs[["__dim_aln_plot"]] <- c(w = 11,
            h = max(3, nblocks * ((nseq + 2) * 20 + 80) / 96))
        p
    },
    height = function() {
        m <- aligned()
        nblocks <- ceiling(ncol(m) / BLOCK)
        max(260, nblocks * ((nrow(m) + 2) * 20 + 80))
    })
}

