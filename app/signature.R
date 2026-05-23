# ─── Scoring Engines ─────────────────────────────────────────────────────────
# Each returns a numeric matrix of signatures (rows) × samples (columns).

run_gsva <- function(mat, sigs, kcdf) {
    if (exists("gsvaParam", where=asNamespace("GSVA"), inherits=FALSE)) {
        par <- GSVA::gsvaParam(exprData=mat, geneSets=sigs, kcdf=kcdf)
        GSVA::gsva(par, verbose=FALSE)
    } else {
        GSVA::gsva(mat, sigs, method="gsva", kcdf=kcdf, verbose=FALSE)
    }
}

run_ssgsea <- function(mat, sigs, normalize) {
    if (exists("ssgseaParam", where=asNamespace("GSVA"), inherits=FALSE)) {
        par <- GSVA::ssgseaParam(exprData=mat, geneSets=sigs, normalize=normalize)
        GSVA::gsva(par, verbose=FALSE)
    } else {
        GSVA::gsva(mat, sigs, method="ssgsea", ssgsea.norm=normalize, verbose=FALSE)
    }
}

run_aucell <- function(mat, sigs, thresh) {
    rankings   <- suppressMessages(AUCell::AUCell_buildRankings(mat, plotStats=FALSE, verbose=FALSE))
    aucMaxRank <- max(1, ceiling(thresh * nrow(mat)))
    cells_AUC  <- suppressMessages(AUCell::AUCell_calcAUC(sigs, rankings, aucMaxRank=aucMaxRank, verbose=FALSE))
    AUCell::getAUC(cells_AUC)
}

run_singscore <- function(mat, sigs, direction) {
    ranked <- singscore::rankGenes(mat)
    rows <- lapply(sigs, function(g) {
        if (direction == "both") {
            singscore::simpleScore(ranked, upSet=g, knownDirection=FALSE)$TotalScore
        } else if (direction == "down") {
            -singscore::simpleScore(ranked, upSet=g)$TotalScore
        } else {
            singscore::simpleScore(ranked, upSet=g)$TotalScore
        }
    })
    out <- do.call(rbind, rows)
    rownames(out) <- names(sigs)
    colnames(out) <- colnames(ranked)
    out
}

# ─── Score Matrix Reactive ─────────────────────────────────────────────────────

ScoreMatrix <- reactiveVal(NULL)

observeEvent(input$run_score, {
    mat  <- ProcessedMatrix()
    shiny::validate(need(!is.null(mat), "Load an expression matrix first (Data & Scoring tab)."))

    # Assemble signatures from the chosen source, then expose them for the preview.
    sigs <- assemble_signatures()
    if (is.null(sigs) || length(sigs) == 0) return()
    LoadedSignatures(sigs)

    # Restrict signatures to genes present; drop those with < 2 matched genes.
    present <- rownames(mat)
    sigs_f  <- lapply(sigs, function(g) intersect(unique(as.character(g)), present))
    keep    <- vapply(sigs_f, length, integer(1)) >= 2
    dropped <- names(sigs_f)[!keep]
    if (length(dropped) > 0) {
        showNotification(
            paste("Skipped (fewer than 2 matched genes):", paste(dropped, collapse=", ")),
            type='warning', duration=8)
    }
    sigs_f <- sigs_f[keep]
    shiny::validate(need(length(sigs_f) >= 1,
        "No signatures had at least 2 genes present in the matrix."))

    method <- input$score_method %||% "gsva"
    withProgress(message=paste("Scoring with", toupper(method)), value=0.3, {
        res <- tryCatch({
            if (method == "gsva") {
                run_gsva(mat, sigs_f, input$gsva_kcdf %||% "Gaussian")
            } else if (method == "ssgsea") {
                run_ssgsea(mat, sigs_f, isTRUE(input$ssgsea_norm))
            } else if (method == "aucell") {
                run_aucell(mat, sigs_f, input$aucell_thresh %||% 0.05)
            } else if (method == "singscore") {
                run_singscore(mat, sigs_f, input$singscore_dir %||% "up")
            }
        }, error=function(e) {
            showNotification(paste("Scoring error:", e$message), type='error', duration=NULL); NULL
        })
        incProgress(0.6)
        if (!is.null(res)) {
            res <- as.matrix(res)
            storage.mode(res) <- "numeric"
            ScoreMatrix(res)
            updateTabsetPanel(session, inputId="ResultTabs", selected="Score Matrix")
            showNotification(
                sprintf("Scored %d signatures × %d samples with %s.",
                        nrow(res), ncol(res), toupper(method)),
                type='message', duration=5)
        }
    })
})

# ─── Score Tab Outputs ──────────────────────────────────────────────────────────

output$score_status <- renderUI({
    sm <- ScoreMatrix()
    if (is.null(sm)) return(helpText("No scores yet. Configure a method and click Run Scoring."))
    tagList(
        strong("Last run:"),
        p(sprintf("%d signatures × %d samples", nrow(sm), ncol(sm)),
          style="margin-top:4px;")
    )
})

output$score_matrix_table <- DT::renderDataTable({
    sm <- ScoreMatrix()
    shiny::validate(need(!is.null(sm), "Run scoring to populate the score matrix."))
    df <- data.frame(Signature=rownames(sm), round(as.data.frame(sm), 4), check.names=FALSE)
    rownames(df) <- NULL
    DT::datatable(df, style='bootstrap', rownames=FALSE,
        options=list(pageLength=15, scrollX=TRUE))
})

# Keep the per-signature selector in sync with the current score matrix.
observe({
    sm <- ScoreMatrix()
    req(sm)
    updateSelectInput(session, "dist_signature",
        choices=rownames(sm), selected=rownames(sm)[1])
})

# ═══════════════════════════════════════════════════════════════════════════════
#  SCORE HEATMAP
# ═══════════════════════════════════════════════════════════════════════════════

# ─── Annotation column selector + dynamic color pickers ─────────────────────────

output$hm_anno_col_selector <- renderUI({
    meta <- MetadataReactive()
    req(meta)
    meta_cols <- colnames(meta)[-1]
    checkboxGroupInput("hm_anno_cols", "Metadata columns to annotate:",
        choices=meta_cols, selected=meta_cols[1])
})

output$anno_color_ui <- renderUI({
    req(input$hm_anno_cols)
    meta <- ProcessedMeta()
    req(meta)
    default_colors <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
                        "#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB")
    ui_elements <- lapply(input$hm_anno_cols, function(col) {
        vals <- unique(as.character(meta[[col]])); vals <- vals[!is.na(vals)]
        if (length(vals) <= 12) {
            color_inputs <- lapply(seq_along(vals), function(i) {
                column(4, colourInput(
                    inputId = paste0("anno_color_", col, "_", gsub("[^A-Za-z0-9]", "_", vals[i])),
                    label   = vals[i],
                    value   = default_colors[((i - 1) %% 12) + 1]))
            })
            tagList(strong(paste("Colors for:", col)), fluidRow(color_inputs))
        } else {
            p(paste(col, "has too many levels for manual colors (> 12)."))
        }
    })
    do.call(tagList, ui_elements)
})

# ─── Builders ───────────────────────────────────────────────────────────────────

# Wrap long labels onto multiple lines. Breaks at "_" / space / "-" / ":" so that
# long MSigDB names (e.g. HALLMARK_INTERFERON_ALPHA_RESPONSE, which have no spaces
# for str_wrap to use) wrap sensibly. ComplexHeatmap renders "\n" as line breaks.
wrap_label <- function(x, width) {
    width <- max(4, as.integer(width))
    vapply(as.character(x), function(s) {
        if (is.na(s) || nchar(s) <= width) return(s)
        parts <- strsplit(s, "(?<=[_ :-])", perl=TRUE)[[1]]   # keep separators
        lines <- character(0); cur <- ""
        for (p in parts) {
            if (nzchar(cur) && nchar(cur) + nchar(p) > width) { lines <- c(lines, cur); cur <- p }
            else cur <- paste0(cur, p)
        }
        if (nzchar(cur)) lines <- c(lines, cur)
        paste(lines, collapse="\n")
    }, character(1), USE.NAMES=FALSE)
}

build_score_color_fun <- function(mat_vals, palette, reverse, centered) {
    rng <- range(mat_vals, na.rm=TRUE)
    pal_colors <- switch(palette,
        "RdBu"    = rev(brewer.pal(11, "RdBu")),   # high = red, low = blue
        "viridis" = viridis::viridis(11),
        "magma"   = viridis::magma(11),
        rev(brewer.pal(11, "RdBu")))
    if (isTRUE(reverse)) pal_colors <- rev(pal_colors)

    if (centered) {
        m <- max(abs(rng), na.rm=TRUE)
        if (!is.finite(m) || m == 0) m <- 1
        colorRamp2(seq(-m, m, length.out=length(pal_colors)), pal_colors)
    } else {
        if (!is.finite(rng[1]) || rng[1] == rng[2]) rng <- c(rng[1] - 1, rng[1] + 1)
        colorRamp2(seq(rng[1], rng[2], length.out=length(pal_colors)), pal_colors)
    }
}

build_column_annotation <- function(meta_df, anno_cols, anno_bar_size, anno_font_size, input) {
    if (is.null(anno_cols) || length(anno_cols) == 0) return(NULL)

    anno_list  <- list()
    color_list <- list()
    for (col in anno_cols) {
        vals <- as.character(meta_df[[col]])
        uniq <- unique(vals[!is.na(vals)])
        col_colors <- sapply(uniq, function(v) {
            id  <- paste0("anno_color_", col, "_", gsub("[^A-Za-z0-9]", "_", v))
            clr <- input[[id]]
            if (is.null(clr)) "#999999" else clr
        })
        names(col_colors) <- uniq
        anno_list[[col]]  <- vals
        color_list[[col]] <- col_colors
    }

    # Match the heatmap legend's font sizes on the annotation-track legends.
    leg_title <- input$hm_legend_title_size %||% 12
    leg_label <- input$hm_legend_label_size %||% 10
    legend_params <- lapply(anno_cols, function(col) {
        list(title_gp = gpar(fontsize=leg_title, fontface="bold"),
             labels_gp = gpar(fontsize=leg_label))
    })
    names(legend_params) <- anno_cols

    HeatmapAnnotation(
        df                      = as.data.frame(anno_list),
        col                     = color_list,
        annotation_height       = unit(anno_bar_size, "mm"),
        annotation_name_gp      = gpar(fontsize=anno_font_size),
        annotation_legend_param = legend_params,
        which                   = "column"
    )
}

# ─── Main Heatmap Reactive ──────────────────────────────────────────────────────

ScoreHeatmapPlotter <- reactive({
    sm <- ScoreMatrix()
    req(sm)
    meta <- tryCatch(ProcessedMeta(), error=function(e) NULL)

    mat_sub    <- sm
    scale_mode <- input$hm_scale %||% "row"
    if (scale_mode == "row" && nrow(mat_sub) >= 1) {
        mat_sub <- t(scale(t(mat_sub))); mat_sub[is.nan(mat_sub)] <- 0
    } else if (scale_mode == "col") {
        mat_sub <- scale(mat_sub); mat_sub[is.nan(mat_sub)] <- 0
    }
    centered <- scale_mode %in% c("row", "col")

    col_fun <- build_score_color_fun(
        mat_vals = mat_sub,
        palette  = input$hm_palette %||% "RdBu",
        reverse  = isTRUE(input$hm_palette_reverse),
        centered = centered
    )

    top_anno <- NULL
    if (isTRUE(input$show_anno) && !is.null(input$hm_anno_cols) &&
        length(input$hm_anno_cols) > 0 && !is.null(meta)) {
        top_anno <- build_column_annotation(
            meta_df       = meta,
            anno_cols     = input$hm_anno_cols,
            anno_bar_size = input$anno_bar_size  %||% 5,
            anno_font_size= input$anno_font_size %||% 10,
            input         = input
        )
    }

    nsig <- nrow(mat_sub)
    nsmp <- ncol(mat_sub)

    # Optional label wrapping for long signature / sample names
    row_labels <- rownames(mat_sub)
    col_labels <- colnames(mat_sub)
    if (isTRUE(input$wrap_row_names)) row_labels <- wrap_label(row_labels, input$wrap_row_width %||% 24)
    if (isTRUE(input$wrap_col_names)) col_labels <- wrap_label(col_labels, input$wrap_col_width %||% 16)

    Heatmap(
        mat_sub,
        col               = col_fun,
        name              = if (scale_mode == "none") "score" else "z-score",
        top_annotation    = top_anno,
        cluster_rows      = isTRUE(input$cluster_rows) && nsig >= 2,
        cluster_columns   = isTRUE(input$cluster_cols) && nsmp >= 2,
        show_row_dend     = isTRUE(input$show_dend),
        show_column_dend  = isTRUE(input$show_dend),
        show_row_names    = isTRUE(input$show_row_names),
        show_column_names = isTRUE(input$show_col_names),
        row_labels        = row_labels,
        column_labels     = col_labels,
        row_names_gp      = gpar(fontsize = input$row_font_size %||% 9),
        column_names_gp   = gpar(fontsize = input$col_font_size %||% 8),
        column_names_rot  = input$col_font_angle %||% 45,
        heatmap_legend_param = list(
            title_gp  = gpar(fontsize = input$hm_legend_title_size %||% 12, fontface = "bold"),
            labels_gp = gpar(fontsize = input$hm_legend_label_size %||% 10)
        ),
        border            = TRUE,
        rect_gp           = gpar(col="white", lwd=0.5)
    )
})

output$score_heatmap_out <- renderPlot({
    ht  <- ScoreHeatmapPlotter()
    pos <- input$hm_legend_pos %||% "right"
    draw(ht, merge_legend=TRUE, heatmap_legend_side=pos, annotation_legend_side=pos)
}, height = function() input$hm_height %||% 600,
   width  = function() input$hm_width  %||% 900)

outputOptions(output, "score_heatmap_out", suspendWhenHidden=FALSE)

output$download_heatmap <- downloadHandler(
    filename = function() paste0("signature_heatmap.", input$hm_download_format %||% "png"),
    content  = function(file) {
        ht   <- ScoreHeatmapPlotter()
        h_in <- (input$hm_height %||% 600) / 96
        w_in <- (input$hm_width  %||% 900) / 96
        fmt  <- input$hm_download_format %||% "png"
        if (fmt == "png") {
            png(file, height=input$hm_height %||% 600, width=input$hm_width %||% 900, res=96)
        } else if (fmt == "jpeg") {
            jpeg(file, height=input$hm_height %||% 600, width=input$hm_width %||% 900, res=96)
        } else if (fmt == "tiff") {
            tiff(file, height=input$hm_height %||% 600, width=input$hm_width %||% 900, res=96)
        } else if (fmt == "pdf") {
            pdf(file, height=h_in, width=w_in)
        } else if (fmt == "svg") {
            svg(file, height=h_in, width=w_in)
        } else if (fmt == "eps") {
            setEPS(); postscript(file, height=h_in, width=w_in)
        }
        pos <- input$hm_legend_pos %||% "right"
        draw(ht, merge_legend=TRUE, heatmap_legend_side=pos, annotation_legend_side=pos)
        dev.off()
    }
)

# ═══════════════════════════════════════════════════════════════════════════════
#  PER-SIGNATURE DISTRIBUTIONS
# ═══════════════════════════════════════════════════════════════════════════════

build_fill_scale <- function(palette_name, custom_colors=NULL) {
    switch(palette_name,
        "default" = scale_fill_discrete(),
        "Set1"    = scale_fill_brewer(palette="Set1"),
        "Set2"    = scale_fill_brewer(palette="Set2"),
        "viridis" = viridis::scale_fill_viridis(discrete=TRUE),
        "custom"  = if (!is.null(custom_colors) && length(custom_colors) > 0)
                        scale_fill_manual(values=custom_colors) else scale_fill_discrete(),
        scale_fill_discrete()
    )
}

output$dist_group_selector <- renderUI({
    meta <- MetadataReactive()
    req(meta)
    meta_cols <- colnames(meta)[-1]
    selectInput("dist_group", "Group by (metadata):", choices=meta_cols, selected=meta_cols[1])
})

output$dist_color_ui <- renderUI({
    req(input$dist_group)
    meta <- ProcessedMeta()
    req(meta)
    vals <- unique(as.character(meta[[input$dist_group]])); vals <- vals[!is.na(vals)]
    base_pal <- brewer.pal(max(3, min(9, length(vals))), "Set1")
    pickers <- lapply(seq_along(vals), function(i) {
        colourInput(
            inputId = paste0("dist_color_", gsub("[^A-Za-z0-9]", "_", vals[i])),
            label   = vals[i],
            value   = base_pal[((i - 1) %% length(base_pal)) + 1])
    })
    do.call(tagList, pickers)
})

DistData <- reactive({
    sm <- ScoreMatrix()
    req(sm)
    sig <- input$dist_signature
    req(sig, sig %in% rownames(sm))
    meta <- ProcessedMeta()
    req(meta)
    grp <- input$dist_group
    req(grp, grp %in% colnames(meta))

    samples <- colnames(sm)
    df <- data.frame(
        Sample = samples,
        Score  = as.numeric(sm[sig, ]),
        Group  = as.character(meta[samples, grp]),
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$Group), ]
    df$Group <- factor(df$Group)
    df
})

DistPlotter <- reactive({
    df <- DistData()
    req(df, nrow(df) > 0)
    sig <- input$dist_signature

    p <- ggplot(df, aes(x=Group, y=Score, fill=Group))

    ptype <- input$dist_plot_type %||% "boxplot"
    if (ptype == "boxplot") {
        p <- p + geom_boxplot(outlier.shape = if (isTRUE(input$dist_jitter)) NA else 19)
    } else if (ptype == "violin") {
        p <- p + geom_violin(trim=FALSE)
    } else {
        p <- p + geom_violin(trim=FALSE) +
                 geom_boxplot(width=0.15, fill="white", alpha=0.7, outlier.shape=NA)
    }
    if (isTRUE(input$dist_jitter)) {
        p <- p + geom_jitter(width=0.15, alpha=0.6, size=1.5)
    }

    if (isTRUE(input$dist_show_stats) && nlevels(df$Group) >= 2) {
        if (isTRUE(input$dist_pairwise)) {
            comps <- combn(levels(df$Group), 2, simplify=FALSE)
            p <- p + stat_compare_means(
                method      = input$dist_pairwise_test %||% "wilcox.test",
                comparisons = comps)
        }
        if (isTRUE(input$dist_global)) {
            gm <- if ((input$dist_global_test %||% "kruskal.test") == "anova") "anova" else "kruskal.test"
            p <- p + stat_compare_means(method=gm, label.y.npc="top")
        }
    }

    custom <- NULL
    if ((input$dist_palette %||% "default") == "custom") {
        lv <- levels(df$Group)
        custom <- sapply(lv, function(v) {
            input[[paste0("dist_color_", gsub("[^A-Za-z0-9]", "_", v))]] %||% "#999999"
        })
        names(custom) <- lv
    }
    p <- p + build_fill_scale(input$dist_palette %||% "default", custom)

    p + labs(title=sig, x=input$dist_group %||% "Group", y="Signature score") +
        theme_bw() +
        theme(
            plot.title  = element_text(size=input$dist_title_font %||% 14, face="bold", hjust=0.5),
            axis.title  = element_text(size=input$dist_axis_font %||% 12),
            axis.text    = element_text(size=input$dist_tick_font %||% 10),
            axis.text.x  = element_text(angle=input$dist_x_angle %||% 45, hjust=1),
            legend.title = element_text(size=input$dist_legend_title_size %||% 12),
            legend.text  = element_text(size=input$dist_legend_label_size %||% 10),
            legend.position = input$dist_legend_pos %||% "right"
        )
})

output$dist_out <- renderPlot({
    DistPlotter()
}, height = function() input$dist_height %||% 550,
   width  = function() input$dist_width  %||% 700)

outputOptions(output, "dist_out", suspendWhenHidden=FALSE)

output$download_dist <- downloadHandler(
    filename = function() paste0("signature_distribution.", input$dist_download_format %||% "png"),
    content  = function(file) {
        p    <- DistPlotter()
        hpx  <- input$dist_height %||% 550
        wpx  <- input$dist_width  %||% 700
        h_in <- hpx / 96
        w_in <- wpx / 96
        fmt  <- input$dist_download_format %||% "png"
        if (fmt == "png") {
            png(file, height=hpx, width=wpx, res=96)
        } else if (fmt == "jpeg") {
            jpeg(file, height=hpx, width=wpx, res=96)
        } else if (fmt == "tiff") {
            tiff(file, height=hpx, width=wpx, res=96)
        } else if (fmt == "pdf") {
            pdf(file, height=h_in, width=w_in)
        } else if (fmt == "svg") {
            svg(file, height=h_in, width=w_in)
        } else if (fmt == "eps") {
            setEPS(); postscript(file, height=h_in, width=w_in)
        }
        print(p)
        dev.off()
    }
)

# ═══════════════════════════════════════════════════════════════════════════════
#  DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════════════

ScoreExport <- reactive({
    sm <- ScoreMatrix()
    req(sm)
    if ((input$dl_orientation %||% "sig_by_sample") == "sig_by_sample") {
        df <- data.frame(Signature=rownames(sm), as.data.frame(sm), check.names=FALSE)
    } else {
        tm <- t(sm)
        df <- data.frame(Sample=rownames(tm), as.data.frame(tm), check.names=FALSE)
    }
    rownames(df) <- NULL
    df
})

output$download_scores <- downloadHandler(
    filename = function() {
        orient <- input$dl_orientation %||% "sig_by_sample"
        ext    <- input$dl_format %||% "csv"
        paste0("signature_scores_", orient, ".", ext)
    },
    content = function(file) {
        df <- ScoreExport()
        sep <- if ((input$dl_format %||% "csv") == "tsv") "\t" else ","
        fwrite(df, file, sep=sep)
    }
)

# ═══════════════════════════════════════════════════════════════════════════════
#  HELP MODAL  (triggered by the floating help button)
# ═══════════════════════════════════════════════════════════════════════════════

show_sig_help_ui <- function() {
    showModal(modalDialog(
        title     = tagList(icon("circle-question"), " Signature Scoring Tool Help"),
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        tabsetPanel(
            tabPanel("Overview",
                br(),
                h4("What this tool does"),
                p("It scores every sample against any number of gene-set signatures, producing a
                   signature × sample score matrix that flows directly into the MaGIC Survival and
                   Correlation tools."),
                h4("Workflow"),
                tags$ol(
                    tags$li("Load an expression matrix + metadata (Data Input)."),
                    tags$li("Define gene-set signatures (Gene Sets)."),
                    tags$li("Pick a scoring method and run it (Score)."),
                    tags$li("Visualize and download (Score Heatmap, Distributions, Download).")
                )
            ),
            tabPanel("Scoring Methods",
                br(),
                tags$ul(
                    tags$li(strong("GSVA —"), " rank-based, cohort-relative; best for stable comparisons across samples. kcdf: Gaussian for continuous/log data, Poisson for raw integer counts."),
                    tags$li(strong("ssGSEA —"), " single-sample GSEA enrichment scores; optional normalization rescales across the dataset."),
                    tags$li(strong("AUCell —"), " area-under-the-recovery-curve over per-sample gene rankings; robust to dropouts. The AUC threshold sets the top fraction of the ranking considered."),
                    tags$li(strong("singscore —"), " rank-based directional scoring. Up = signature genes expected highly expressed; Down = expected lowly expressed (score negated); Both = undirected, scoring coordinated extremes.")
                )
            ),
            tabPanel("Gene Sets & Input",
                br(),
                h4("Expression matrix"),
                p("CSV/TSV, genes × samples, gene IDs in the first column, normalized (log-scale preferred). Metadata: samples × variables, sample names in the first column matching the matrix columns."),
                h4("Signature sources"),
                tags$ul(
                    tags$li(strong("GMT upload —"), " standard Broad format (name, description, then genes per line)."),
                    tags$li(strong("MSigDB —"), " pick a collection (Hallmark, C2 KEGG/Reactome, C5 GO:BP, C7) and organism, then choose specific signatures."),
                    tags$li(strong("Paste —"), " name a signature, paste its genes, and add it; repeat for as many as you like.")
                ),
                p("The preview table shows, per signature, how many genes are in the set and how many are present in your matrix. Signatures with fewer than 2 matched genes are skipped at scoring time.")
            ),
            tabPanel("Plots & Output",
                br(),
                h4("Score Heatmap"),
                p("Signatures (rows) × samples (columns). Scale per row, per column, or use raw scores; choose a palette; add metadata annotation tracks with per-level colors; toggle clustering, dendrograms, and label wrapping for long names."),
                h4("Distributions"),
                p("Pick a signature and a metadata grouping to see its score distribution as a boxplot/violin, with optional pairwise/global statistics and per-group colors."),
                h4("Reproducible code"),
                p("Each plot has a ", icon("file-code"), " button that opens the R code to reproduce your current view offline."),
                h4("Output orientation"),
                p("The canonical download is signatures × samples (signature name in the first column), which the Survival and Correlation tools read directly as a feature × sample matrix. A samples × signatures orientation is also available for joining with metadata tables.")
            )
        )
    ))
}

observeEvent(input$show_help_float, { show_sig_help_ui() })

# ═══════════════════════════════════════════════════════════════════════════════
#  CODE MODALS  (reproducible R code for each plot)
# ═══════════════════════════════════════════════════════════════════════════════

show_code_modal_ui <- function(code) {
    showModal(modalDialog(
        title     = tagList(icon("file-code"), " Reproducible R Code"),
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        p("Copy this code to reproduce your current plot in an offline R session.",
          style="color:#555; margin-bottom:12px;"),
        tags$pre(
            style = paste(
                "background:#1e1e1e; color:#d4d4d4; border-radius:6px;",
                "padding:16px; font-size:12px; max-height:520px; overflow-y:auto;",
                "white-space:pre; font-family:'Courier New', monospace;"
            ),
            code
        )
    ))
}

build_hm_code <- function(sm, meta, inp) {
    scale_mode <- inp$hm_scale %||% "row"
    pal        <- inp$hm_palette %||% "RdBu"
    centered   <- scale_mode %in% c("row", "col")

    scale_block <- switch(scale_mode,
        "row" = "mat <- t(scale(t(mat))); mat[is.nan(mat)] <- 0   # row z-score\n",
        "col" = "mat <- scale(mat); mat[is.nan(mat)] <- 0          # column z-score\n",
        "")
    pal_line <- switch(pal,
        "viridis" = "pal <- viridis::viridis(11)",
        "magma"   = "pal <- viridis::magma(11)",
        'pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))')
    rev_line <- if (isTRUE(inp$hm_palette_reverse)) "\npal <- rev(pal)" else ""
    col_block <- if (centered)
        "m <- max(abs(range(mat, na.rm = TRUE)))\ncol_fun <- circlize::colorRamp2(seq(-m, m, length.out = length(pal)), pal)\n"
    else
        "rng <- range(mat, na.rm = TRUE)\ncol_fun <- circlize::colorRamp2(seq(rng[1], rng[2], length.out = length(pal)), pal)\n"

    anno_cols <- inp$hm_anno_cols
    if (isTRUE(inp$show_anno) && !is.null(anno_cols) && length(anno_cols) > 0 && !is.null(meta)) {
        col_defs <- vapply(anno_cols, function(col) {
            vals   <- unique(as.character(meta[[col]])); vals <- vals[!is.na(vals)]
            colors <- vapply(vals, function(v) {
                clr <- inp[[paste0("anno_color_", col, "_", gsub("[^A-Za-z0-9]", "_", v))]]
                if (is.null(clr)) "#999999" else clr
            }, character(1))
            sprintf('        %s = c(%s)', col,
                    paste0('"', vals, '" = "', colors, '"', collapse = ", "))
        }, character(1))
        df_args  <- paste(sprintf('%s = meta[["%s"]]', anno_cols, anno_cols), collapse = ", ")
        anno_block <- paste0(
            "top_anno <- ComplexHeatmap::HeatmapAnnotation(\n",
            "    df  = data.frame(", df_args, "),\n",
            "    col = list(\n", paste(col_defs, collapse = ",\n"), "\n    ),\n",
            "    which = \"column\"\n)\n")
        anno_arg <- "    top_annotation  = top_anno,\n"
    } else { anno_block <- ""; anno_arg <- "" }

    paste0(
        "library(ComplexHeatmap)\nlibrary(circlize)\nlibrary(RColorBrewer)\n\n",
        "# scores: signature x sample matrix produced by this tool\nmat <- scores\n\n",
        if (nzchar(scale_block)) paste0("# Scaling\n", scale_block, "\n") else "",
        "# Colour scale\n", pal_line, rev_line, "\n", col_block, "\n",
        if (nzchar(anno_block)) paste0("# Column annotation (metadata)\n", anno_block, "\n") else "",
        "Heatmap(\n",
        "    mat,\n",
        "    col             = col_fun,\n",
        sprintf('    name            = "%s",\n', if (scale_mode == "none") "score" else "z-score"),
        anno_arg,
        sprintf("    cluster_rows    = %s,\n", tolower(as.character(isTRUE(inp$cluster_rows)))),
        sprintf("    cluster_columns = %s,\n", tolower(as.character(isTRUE(inp$cluster_cols)))),
        sprintf("    show_row_dend    = %s,\n", tolower(as.character(isTRUE(inp$show_dend)))),
        sprintf("    show_column_dend = %s,\n", tolower(as.character(isTRUE(inp$show_dend)))),
        sprintf("    row_names_gp    = gpar(fontsize = %s),\n", as.character(inp$row_font_size %||% 9)),
        sprintf("    column_names_gp = gpar(fontsize = %s),\n", as.character(inp$col_font_size %||% 8)),
        sprintf("    column_names_rot = %s,\n", as.character(inp$col_font_angle %||% 45)),
        "    heatmap_legend_param = list(\n",
        sprintf("        title_gp  = gpar(fontsize = %s, fontface = \"bold\"),\n", as.character(inp$hm_legend_title_size %||% 12)),
        sprintf("        labels_gp = gpar(fontsize = %s)\n", as.character(inp$hm_legend_label_size %||% 10)),
        "    ),\n",
        "    border = TRUE, rect_gp = gpar(col = \"white\", lwd = 0.5)\n",
        sprintf(") |> ComplexHeatmap::draw(merge_legend = TRUE, heatmap_legend_side = \"%s\", annotation_legend_side = \"%s\")\n",
                inp$hm_legend_pos %||% "right", inp$hm_legend_pos %||% "right")
    )
}

build_dist_code <- function(inp) {
    sig   <- inp$dist_signature %||% "SIGNATURE"
    grp   <- inp$dist_group %||% "Group"
    ptype <- inp$dist_plot_type %||% "boxplot"

    geom_block <- switch(ptype,
        "boxplot"    = sprintf("    geom_boxplot(outlier.shape = %s) +\n", if (isTRUE(inp$dist_jitter)) "NA" else "19"),
        "violin"     = "    geom_violin(trim = FALSE) +\n",
        "violin_box" = "    geom_violin(trim = FALSE) +\n    geom_boxplot(width = 0.15, fill = \"white\", alpha = 0.7, outlier.shape = NA) +\n",
        "    geom_boxplot() +\n")
    jitter_line <- if (isTRUE(inp$dist_jitter)) "    geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +\n" else ""

    stats_block <- ""
    if (isTRUE(inp$dist_show_stats)) {
        if (isTRUE(inp$dist_pairwise))
            stats_block <- paste0(stats_block, sprintf(
                '    ggpubr::stat_compare_means(method = "%s",\n        comparisons = combn(levels(df$Group), 2, simplify = FALSE)) +\n',
                inp$dist_pairwise_test %||% "wilcox.test"))
        if (isTRUE(inp$dist_global)) {
            gm <- if ((inp$dist_global_test %||% "kruskal.test") == "anova") "anova" else "kruskal.test"
            stats_block <- paste0(stats_block, sprintf(
                '    ggpubr::stat_compare_means(method = "%s", label.y.npc = "top") +\n', gm))
        }
    }
    fill_block <- switch(inp$dist_palette %||% "default",
        "Set1"    = '    scale_fill_brewer(palette = "Set1") +\n',
        "Set2"    = '    scale_fill_brewer(palette = "Set2") +\n',
        "viridis" = '    viridis::scale_fill_viridis(discrete = TRUE) +\n',
        "custom"  = '    scale_fill_manual(values = my_colors) +   # set per-group colours\n',
        "")

    paste0(
        "library(ggplot2)\nlibrary(ggpubr)\n\n",
        "# scores: signature x sample matrix; metadata: samples x variables\n",
        sprintf('df <- data.frame(\n    Score = scores["%s", ],\n    Group = factor(metadata[["%s"]])\n)\n\n', sig, grp),
        "ggplot(df, aes(x = Group, y = Score, fill = Group)) +\n",
        geom_block, jitter_line, stats_block, fill_block,
        sprintf('    labs(title = "%s", x = "%s", y = "Signature score") +\n', sig, grp),
        "    theme_bw() +\n",
        "    theme(\n",
        sprintf("        legend.title = element_text(size = %s),\n", as.character(inp$dist_legend_title_size %||% 12)),
        sprintf("        legend.text  = element_text(size = %s),\n", as.character(inp$dist_legend_label_size %||% 10)),
        sprintf('        legend.position = "%s"\n', inp$dist_legend_pos %||% "right"),
        "    )\n"
    )
}

observeEvent(input$show_hm_code, {
    sm   <- isolate(tryCatch(ScoreMatrix(),  error=function(e) NULL))
    meta <- isolate(tryCatch(ProcessedMeta(), error=function(e) NULL))
    code <- if (is.null(sm)) "# Run scoring first, then click here for the reproducible code."
            else tryCatch(build_hm_code(sm, meta, input),
                          error=function(e) paste0("# Code generation error: ", conditionMessage(e)))
    show_code_modal_ui(code)
})

observeEvent(input$show_dist_code, {
    sm   <- isolate(tryCatch(ScoreMatrix(), error=function(e) NULL))
    code <- if (is.null(sm)) "# Run scoring first, then click here for the reproducible code."
            else tryCatch(build_dist_code(input),
                          error=function(e) paste0("# Code generation error: ", conditionMessage(e)))
    show_code_modal_ui(code)
})
