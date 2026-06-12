# ─── Utilities ─────────────────────────────────────────────────────────────────

`%||%` <- function(a, b) {
    if (is.null(a)) return(b)
    if (length(a) == 0) return(b)
    if (length(a) == 1 && is.character(a) && !nzchar(a)) return(b)
    a
}

read_delim_auto <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("tsv", "txt")) {
        fread(path, sep="\t")
    } else {
        fread(path, sep=",")
    }
}

# Parse a Broad-format GMT file into a named list of gene-symbol vectors.
parse_gmt <- function(path) {
    lines <- readLines(path, warn=FALSE)
    lines <- lines[nzchar(trimws(lines))]
    parsed <- lapply(lines, function(ln) {
        parts <- strsplit(ln, "\t")[[1]]
        if (length(parts) < 3) return(NULL)
        genes <- trimws(parts[-c(1, 2)])
        genes <- unique(genes[nzchar(genes)])
        if (length(genes) == 0) return(NULL)
        list(name = trimws(parts[1]), genes = genes)
    })
    parsed <- Filter(Negate(is.null), parsed)
    out <- lapply(parsed, `[[`, "genes")
    names(out) <- make.unique(vapply(parsed, `[[`, character(1), "name"))
    out
}

# Map a collection code + species to an msigdbr data.frame (gs_name, gene_symbol).
# msigdbr >= 10 renamed the arguments to collection/subcollection and split KEGG
# into CP:KEGG_LEGACY / CP:KEGG_MEDICUS. We target that API and fall back to the
# legacy category/subcategory arguments so the tool survives msigdbr version drift.
msigdbr_fetch <- function(species, collection) {
    spec <- switch(collection,
        "H"           = list(coll="H",  sub=NULL),
        "C2_KEGG"     = list(coll="C2", sub="CP:KEGG_LEGACY"),
        "C2_REACTOME" = list(coll="C2", sub="CP:REACTOME"),
        "C5_BP"       = list(coll="C5", sub="GO:BP"),
        "C7"          = list(coll="C7", sub=NULL),
        list(coll="H", sub=NULL)
    )
    args_new <- c(list(species=species, collection=spec$coll),
                  if (!is.null(spec$sub)) list(subcollection=spec$sub))
    tryCatch(
        do.call(msigdbr, args_new),
        error = function(e) {
            sub_old <- if (identical(spec$sub, "CP:KEGG_LEGACY")) "CP:KEGG" else spec$sub
            args_old <- c(list(species=species, category=spec$coll),
                          if (!is.null(sub_old)) list(subcategory=sub_old))
            do.call(msigdbr, args_old)
        }
    )
}

# ─── Demo Data Generator ───────────────────────────────────────────────────────

DemoDataCache <- reactiveVal(NULL)

generate_demo_data <- function() {
    set.seed(42)
    n_per  <- 10
    groups <- rep(c("Control", "TreatA", "TreatB"), each=n_per)
    n      <- length(groups)
    sample_names <- paste0("Sample_", sprintf("%02d", 1:n))
    sex    <- rep(c("Male", "Female"), length.out=n)
    batch  <- rep(c("Batch1", "Batch2"), length.out=n)

    modules <- list(
        DEMO_CELL_CYCLE = c("CDK1","CCNB1","AURKA","PLK1","TOP2A","MKI67","BUB1","CCNA2",
                            "CDC20","BIRC5","KIF11","CENPA","MCM2","MCM5","PCNA"),
        DEMO_INFLAMMATORY = c("TNF","IL6","CCL2","CXCL8","TLR4","NFKB1","IL1B","ICAM1",
                              "VCAM1","CCL5","CXCL10","IRF1","STAT1","CD14","CASP1"),
        DEMO_METABOLISM = c("GAPDH","PKM","LDHA","HK2","ENO1","PGK1","PFKL","GPI",
                            "ALDOA","TPI1","PGAM1","MDH2","IDH1","CS","SDHA")
    )
    # Each module is elevated in one group
    elevate <- c(DEMO_CELL_CYCLE="TreatA", DEMO_INFLAMMATORY="TreatB", DEMO_METABOLISM="Control")

    rows <- list()
    gene_names <- character(0)
    for (mod in names(modules)) {
        for (g in modules[[mod]]) {
            base  <- rnorm(n, mean=8, sd=1)
            base[groups == elevate[[mod]]] <- base[groups == elevate[[mod]]] + 2.5
            rows[[length(rows) + 1]] <- round(base, 3)
            gene_names <- c(gene_names, g)
        }
    }
    # Background genes with no module structure
    n_bg <- 150
    for (i in seq_len(n_bg)) {
        rows[[length(rows) + 1]] <- round(rnorm(n, mean=7, sd=1.5), 3)
        gene_names <- c(gene_names, paste0("GENE", i))
    }

    expr_mat <- do.call(rbind, rows)
    expr_dt  <- data.table(Gene = gene_names, as.data.table(expr_mat))
    setnames(expr_dt, c("Gene", sample_names))

    meta_dt <- data.table(
        Sample = sample_names,
        Group  = groups,
        Sex    = sex,
        Batch  = batch
    )

    list(matrix = expr_dt, metadata = meta_dt, signatures = modules)
}

# ─── Data Loading Reactives ─────────────────────────────────────────────────────

MatrixReactive <- reactive({
    if (input$DemoData == FALSE) {
        cached <- DemoDataCache()
        if (is.null(cached)) { cached <- generate_demo_data(); DemoDataCache(cached) }
        cached$matrix
    } else {
        shiny::validate(need(!is.null(input$matrix_file), "Please upload an expression matrix file."))
        tryCatch(
            read_delim_auto(input$matrix_file$datapath),
            error = function(e) { showNotification(paste("Matrix parse error:", e$message), type='error', duration=NULL); NULL }
        )
    }
})

MetadataReactive <- reactive({
    if (input$DemoData == FALSE) {
        cached <- DemoDataCache()
        if (is.null(cached)) { cached <- generate_demo_data(); DemoDataCache(cached) }
        cached$metadata
    } else {
        shiny::validate(need(!is.null(input$metadata_file), "Please upload a metadata file."))
        tryCatch(
            read_delim_auto(input$metadata_file$datapath),
            error = function(e) { showNotification(paste("Metadata parse error:", e$message), type='error', duration=NULL); NULL }
        )
    }
})

# ─── Gene Column Selector (for custom uploads) ──────────────────────────────────

output$gene_col_selector <- renderUI({
    req(input$matrix_file)
    mat <- MatrixReactive()
    req(mat)
    tagList(
        selectInput("gene_col", "Which column contains gene names?",
            choices = colnames(mat), selected = colnames(mat)[1]),
        hr()
    )
})

# ─── Preview Tables ──────────────────────────────────────────────────────────────

output$matrix_table <- DT::renderDataTable({
    mat <- MatrixReactive()
    req(mat)
    DT::datatable(mat, style='bootstrap', options=list(pageLength=15, scrollX=TRUE))
})

output$metadata_table <- DT::renderDataTable({
    meta <- MetadataReactive()
    req(meta)
    DT::datatable(meta, style='bootstrap', options=list(pageLength=15, scrollX=TRUE))
})

# ─── Processed Matrix (numeric, rownames = genes) ───────────────────────────────

ProcessedMatrix <- reactive({
    mat <- MatrixReactive()
    req(mat)

    gene_col <- if (!is.null(input$gene_col) && input$DemoData) input$gene_col else colnames(mat)[1]

    gene_names <- as.character(mat[[gene_col]])
    num_cols   <- setdiff(colnames(mat), gene_col)

    m <- as.matrix(mat[, ..num_cols])
    rownames(m) <- make.unique(gene_names)
    storage.mode(m) <- "numeric"

    # Restrict to samples present in BOTH the matrix and the metadata (shared set,
    # in matrix-column order). Matches the qc/deg idiom and keeps NA-metadata
    # samples out of scoring/annotations instead of showing them as an "NA" group.
    meta_raw <- MetadataReactive()
    if (!is.null(meta_raw) && ncol(meta_raw) >= 1) {
        sample_ids <- as.character(meta_raw[[colnames(meta_raw)[1]]])
        shared <- intersect(colnames(m), sample_ids)
        shiny::validate(need(length(shared) >= 1,
            "No samples are shared between the matrix columns and the metadata sample names."))
        if (length(shared) < ncol(m)) {
            showNotification(sprintf("Using %d of %d matrix samples that are present in the metadata.",
                length(shared), ncol(m)), type='warning', duration=8)
        }
        m <- m[, shared, drop=FALSE]
    }
    m
})

# ─── Aligned Metadata (rows = samples, ordered to match matrix columns) ─────────

ProcessedMeta <- reactive({
    meta  <- MetadataReactive()
    mat_m <- ProcessedMatrix()
    req(meta, mat_m)

    sample_col   <- colnames(meta)[1]
    sample_names <- as.character(meta[[sample_col]])

    col_order <- colnames(mat_m)
    idx       <- match(col_order, sample_names)

    if (any(is.na(idx))) {
        showNotification(
            paste("Warning: Some matrix sample names were not found in metadata:",
                  paste(col_order[is.na(idx)], collapse=", ")),
            type='warning', duration=8
        )
    }

    df <- as.data.frame(meta)[idx, , drop=FALSE]
    rownames(df) <- col_order
    df
})

# ─── Gene Set Source Selector (dynamic: demo option only with demo data) ────────

output$gs_mode_ui <- renderUI({
    if (isTRUE(input$DemoData)) {
        # Custom data uploaded — no built-in demo signatures
        radioButtons("gs_mode", "Signature source:",
            choices=c("GMT file upload"="gmt",
                      "MSigDB collections"="msigdb",
                      "Paste custom signatures"="paste"),
            selected="msigdb")
    } else {
        # Demo data — offer the matching built-in signatures, selected by default
        radioButtons("gs_mode", "Signature source:",
            choices=c("Demo signatures (built-in)"="demo",
                      "GMT file upload"="gmt",
                      "MSigDB collections"="msigdb",
                      "Paste custom signatures"="paste"),
            selected="demo")
    }
})

# ─── MSigDB Signature Picker Population ─────────────────────────────────────────

observe({
    req(input$gs_mode == "msigdb")
    tryCatch({
        species_sel <- input$msigdb_species %||% "Homo sapiens"
        coll        <- input$msigdb_collection %||% "H"
        gs_df       <- msigdbr_fetch(species_sel, coll)
        gene_sets   <- sort(unique(gs_df$gs_name))
        updateSelectizeInput(session, "msigdb_sets", choices=gene_sets, server=TRUE,
            options=list(placeholder='Select one or more gene sets...', maxOptions=2000))
    }, error=function(e) {
        showNotification(paste("MSigDB query error:", e$message), type='error', duration=8)
    })
})

# ─── Paste-in Signatures (staged) ───────────────────────────────────────────────

PastedSignatures <- reactiveVal(list())

observeEvent(input$add_paste, {
    nm  <- trimws(input$paste_name %||% "")
    raw <- input$paste_genes %||% ""
    shiny::validate(need(nchar(nm) > 0, "Please provide a signature name."))
    genes <- trimws(unlist(strsplit(raw, "[,\n\t ]+")))
    genes <- unique(genes[nzchar(genes)])
    if (length(genes) == 0) {
        showNotification("No gene symbols found in the pasted text.", type='warning', duration=5)
        return()
    }
    current <- PastedSignatures()
    current[[nm]] <- genes
    PastedSignatures(current)
    updateTextInput(session, "paste_name", value="")
    updateTextAreaInput(session, "paste_genes", value="")
    showNotification(paste0("Added '", nm, "' (", length(genes), " genes)."), type='message', duration=4)
})

observeEvent(input$clear_paste, {
    PastedSignatures(list())
    showNotification("Cleared staged signatures.", type='message', duration=3)
})

output$paste_staged <- renderText({
    sigs <- PastedSignatures()
    if (length(sigs) == 0) return("(none yet)")
    paste(sprintf("%s (%d genes)", names(sigs), vapply(sigs, length, integer(1))), collapse="\n")
})

# ─── Loaded Signatures ──────────────────────────────────────────────────────────

LoadedSignatures <- reactiveVal(NULL)

# Assemble a named list of gene-symbol vectors from the active source, or return
# NULL (with a notification) if nothing usable is available. Called by the single
# Run button on the Data & Scoring tab.
assemble_signatures <- function() {
    mode <- input$gs_mode
    sigs <- NULL

    if (is.null(mode) || mode == "demo") {
        cached <- DemoDataCache()
        if (is.null(cached)) { cached <- generate_demo_data(); DemoDataCache(cached) }
        sigs <- cached$signatures

    } else if (mode == "gmt") {
        if (is.null(input$gmt_file)) {
            showNotification("Please upload a GMT file.", type='warning', duration=6); return(NULL)
        }
        sigs <- tryCatch(parse_gmt(input$gmt_file$datapath),
            error=function(e) { showNotification(paste("GMT parse error:", e$message), type='error', duration=NULL); NULL })

    } else if (mode == "msigdb") {
        sigs <- tryCatch({
            species_sel <- input$msigdb_species %||% "Homo sapiens"
            coll        <- input$msigdb_collection %||% "H"
            gs_df       <- msigdbr_fetch(species_sel, coll)
            picks       <- input$msigdb_sets
            if (is.null(picks) || length(picks) == 0) {
                if (coll == "H") {
                    picks <- sort(unique(gs_df$gs_name))   # all 50 Hallmark sets
                } else {
                    showNotification("Select at least one gene set for this collection.", type='warning', duration=6)
                    return(NULL)
                }
            }
            gs_df <- gs_df[gs_df$gs_name %in% picks, ]
            lapply(split(gs_df$gene_symbol, gs_df$gs_name), function(g) unique(as.character(g)))
        }, error=function(e) { showNotification(paste("MSigDB error:", e$message), type='error', duration=NULL); NULL })

    } else if (mode == "paste") {
        sigs <- PastedSignatures()
        if (length(sigs) == 0) {
            showNotification("Add at least one pasted signature first.", type='warning', duration=6); return(NULL)
        }
    }

    if (is.null(sigs) || length(sigs) == 0) {
        showNotification("No signatures were loaded.", type='warning', duration=6); return(NULL)
    }
    sigs
}

# ─── Gene Set Preview Table ─────────────────────────────────────────────────────

GenesetPreviewDF <- reactive({
    sigs <- LoadedSignatures()
    req(sigs, length(sigs) > 0)

    mat     <- tryCatch(ProcessedMatrix(), error=function(e) NULL)
    present <- if (!is.null(mat)) rownames(mat) else character(0)

    df <- data.frame(
        Signature = names(sigs),
        `Genes in set` = vapply(sigs, length, integer(1)),
        `Genes found in matrix` = vapply(sigs, function(g) sum(g %in% present), integer(1)),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    df$`Genes missing` <- df$`Genes in set` - df$`Genes found in matrix`
    rownames(df) <- NULL
    df
})

output$geneset_preview <- DT::renderDataTable({
    shiny::validate(need(!is.null(LoadedSignatures()) && length(LoadedSignatures()) > 0,
        "No signatures yet. Choose a source on the left and click 'Run Scoring'."))
    df <- GenesetPreviewDF()
    DT::datatable(df, style='bootstrap', rownames=FALSE,
        options=list(pageLength=15, scrollX=TRUE)) |>
        DT::formatStyle('Genes found in matrix',
            color=DT::styleInterval(0, c('#b94a48', '#3c763d')))
})

# ─── Table Download Handlers ────────────────────────────────────────────────────

output$download_matrix <- downloadHandler(
    filename = function() "expression_matrix.csv",
    content  = function(file) {
        mat <- MatrixReactive(); req(mat)
        fwrite(mat, file, sep=",")
    }
)

output$download_metadata <- downloadHandler(
    filename = function() "sample_metadata.csv",
    content  = function(file) {
        meta <- MetadataReactive(); req(meta)
        fwrite(meta, file, sep=",")
    }
)

output$download_signatures <- downloadHandler(
    filename = function() "loaded_signatures_summary.csv",
    content  = function(file) {
        df <- GenesetPreviewDF(); req(df)
        fwrite(df, file, sep=",")
    }
)
