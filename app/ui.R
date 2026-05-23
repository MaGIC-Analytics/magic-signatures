library(shiny)
require(shinyjs)
library(shinythemes)
require(shinycssloaders)
library(shinyWidgets)

library(DT)
library(tidyverse)
library(data.table)
library(colourpicker)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)
library(ggpubr)

# Scoring engines
library(GSVA)
library(AUCell)
library(singscore)
library(GSEABase)

tagList(
    tags$head(
        includeHTML(("www/GA.html")),
        tags$style(type = 'text/css','.navbar-brand{display:none;}'),
        tags$style(HTML("
            .control-group-panel {
                border: 1px solid #ddd;
                border-radius: 6px;
                padding: 10px 12px;
                margin-bottom: 10px;
                background-color: #f9f9f9;
            }
            .control-group-title {
                font-weight: bold;
                font-size: 14px;
                color: #0F344C;
                margin-bottom: 8px;
            }
            #show_help_float {
                position: fixed;
                bottom: 28px;
                right: 28px;
                z-index: 9999;
                border-radius: 50%;
                width: 46px;
                height: 46px;
                font-size: 20px;
                padding: 0;
                box-shadow: 0 3px 8px rgba(0,0,0,0.25);
            }
        "))
    ),
    ## Global always-visible help button (fixed bottom-right)
    actionButton("show_help_float", label=NULL,
        icon=icon("circle-question"),
        title="Help & documentation",
        class="btn btn-info"
    ),
    fluidPage(theme = shinytheme('yeti'),
            windowTitle = "MaGIC Signature Scoring Tool",
            useShinyjs(),
            titlePanel(
                fluidRow(
                column(2, tags$a(href='http://www.bioinformagic.io/', tags$img(height=75, src="MaGIC_Icon_0f344c.svg")), align='center'),
                column(10, fluidRow(
                    column(10, h1(strong('MaGIC Signature Scoring Tool'), align='center', style="color:#0F344C;"))
                ))
                ),
                windowTitle = "MaGIC Signature Scoring Tool"),
                tags$style(type='text/css', '.navbar{font-size:20px;}'),
                tags$style(type='text/css', '.nav-tabs{padding-bottom:20px;}'),
                tags$style(type='text/css', '.navbar-default{background-color:#0F344C;}'),
                tags$style(type='text/css', HTML('.navbar { background-color: #0F344C;}
                          .tab-panel{ background-color: #0F344C;}
                          .navbar-default .navbar-nav > .active > a,
                           .navbar-default .navbar-nav > .active > a:focus,
                           .navbar-default .navbar-nav > .active > a:hover {
                                color: white;
                                background-color: #008cba;
                            }')
                          ),
                tags$head(tags$style(".modal-dialog{ width:1300px}")),

        navbarPage(title="", id='NAVTABS',

        ## Intro Page
##########################################################################################################################################################
            tabPanel('Introduction',
                fluidRow(
                    column(2),
                    column(8,
                        column(12, align="center", style="margin-bottom:25px;",
                            h3(markdown("Welcome to the Signature Scoring Tool by the
                            [Molecular and Genomics Informatics Core (MaGIC)](http://www.bioinformagic.io)."))),
                        hr(),
                        p("Score every sample against any number of gene-set signatures, producing a
                           signature × sample matrix that flows directly into the MaGIC Survival
                           and Correlation tools."),
                        h4("How to Use This Tool", style="color:#0F344C;"),
                        tags$ol(
                            tags$li(strong("Open the Data & Scoring tab."),
                                " Use the built-in demo data, or switch on 'Upload custom data' to load a normalized expression matrix (genes × samples, log-scale preferred) and sample metadata."),
                            tags$li(strong("Choose a signature source."),
                                " Built-in demo signatures, a GMT file, MSigDB collections, or pasted custom signatures."),
                            tags$li(strong("Pick a scoring method."),
                                " GSVA, ssGSEA, AUCell, or singscore, with method-specific parameters."),
                            tags$li(strong("Click Run Scoring."),
                                " The Loaded Signatures and Score Matrix sub-tabs populate (each table has its own download button), and the Score Heatmap and Distributions tabs appear."),
                            tags$li(strong("Explore & download."),
                                " View the annotated signature × sample heatmap and per-signature distributions; download the score matrix in the orientation Survival/Correlation expect, plus your plots.")
                        ),
                        hr(),
                        h4("Which scoring method should I use?", style="color:#0F344C;"),
                        fluidRow(
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("GSVA"), style="color:#0F344C;"),
                                    p("Non-parametric, sample-rank based enrichment. Scores are computed relative to the rest of the cohort, making it ideal for ", strong("stable comparative analysis"), " across samples. Use the Gaussian kcdf for continuous data (log-CPM/TPM/VST) and Poisson for integer counts.")
                                ),
                                div(class="control-group-panel",
                                    h5(strong("AUCell"), style="color:#0F344C;"),
                                    p("Rank-based area-under-the-curve scoring. Robust to ", strong("dropouts and varying detection depth"), " because it only considers whether signature genes fall in the top fraction of each sample's ranking.")
                                )
                            ),
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("ssGSEA"), style="color:#0F344C;"),
                                    p("Single-sample GSEA. Produces a ", strong("direct per-sample enrichment score"), " for each signature, computed independently per sample. Optional normalization rescales scores across the dataset.")
                                ),
                                div(class="control-group-panel",
                                    h5(strong("singscore"), style="color:#0F344C;"),
                                    p("Rank-based ", strong("directional single-sample"), " scoring. Choose whether each signature is interpreted as up-regulated, down-regulated, or undirected (bidirectional).")
                                )
                            )
                        ),
                        hr(),
                        h4("Required Input Data Formats", style="color:#0F344C;"),
                        fluidRow(
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("Expression Matrix"), style="color:#0F344C;"),
                                    tags$ul(
                                        tags$li("File format: CSV or TSV"),
                                        tags$li("Rows: Genes (one gene per row)"),
                                        tags$li("Columns: Samples (one sample per column)"),
                                        tags$li("First column: Gene identifiers (symbols, matched to gene sets)"),
                                        tags$li("Values: normalized, log-scale preferred (e.g. log2-CPM, VST) — the same schema as the MaGIC QC tool's output")
                                    ),
                                    tags$pre("Gene,   Sample1, Sample2\nBRCA1,  6.5,     7.1\nTP53,   8.2,     7.9")
                                )
                            ),
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("Sample Metadata"), style="color:#0F344C;"),
                                    tags$ul(
                                        tags$li("File format: CSV or TSV"),
                                        tags$li("Rows: Samples (one sample per row)"),
                                        tags$li("First column: Sample names — must match matrix column names exactly"),
                                        tags$li("Additional columns: Categorical or numeric metadata (e.g. Group, Sex, Batch)")
                                    ),
                                    tags$pre("Sample,  Group,   Sex\nSample1, Control, Male\nSample2, Treated, Female")
                                )
                            )
                        ),
                        hr()
                    ),
                    column(2)
                )
            ),


        ## Data & Scoring Page
##########################################################################################################################################################
            tabPanel('Data & Scoring',
                fluidRow(
                    column(3,
                        wellPanel(

                            ## ── 1. Data ──
                            h4('1. Data', style="color:#0F344C; margin-top:2px;"),
                            materialSwitch("DemoData", label="Upload custom data", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.DemoData",
                                fileInput('matrix_file', 'Expression Matrix (CSV/TSV)',
                                    accept=c('text/csv', 'text/comma-separated-values, text/plain', '.csv',
                                             'text/tsv', 'text/tab-separated-values, text/plain', '.tsv'),
                                    multiple=FALSE
                                ),
                                fileInput('metadata_file', 'Sample Metadata (CSV/TSV)',
                                    accept=c('text/csv', 'text/comma-separated-values, text/plain', '.csv',
                                             'text/tsv', 'text/tab-separated-values, text/plain', '.tsv'),
                                    multiple=FALSE
                                ),
                                uiOutput('gene_col_selector')
                            ),
                            conditionalPanel("input.DemoData==false",
                                p(em("Pre-loaded demo: 3 expression modules + background genes across 30 samples in 3 treatment groups, with matching built-in signatures."),
                                  style="font-size:12px; color:#777;")
                            ),
                            hr(),

                            ## ── 2. Gene set signatures ──
                            h4('2. Gene set signatures', style="color:#0F344C;"),
                            uiOutput('gs_mode_ui'),

                            conditionalPanel("input.gs_mode == 'demo'",
                                helpText("Scoring the 3 built-in demo signatures (DEMO_CELL_CYCLE, DEMO_INFLAMMATORY, DEMO_METABOLISM).")
                            ),
                            conditionalPanel("input.gs_mode == 'gmt'",
                                fileInput('gmt_file', 'Upload GMT File (Broad format)',
                                    accept=c('.gmt', 'text/plain'), multiple=FALSE),
                                helpText("Each line: signature name, description, then tab-separated gene symbols.")
                            ),
                            conditionalPanel("input.gs_mode == 'msigdb'",
                                selectInput("msigdb_species", "Organism:",
                                    choices=c("Human (Homo sapiens)"="Homo sapiens",
                                              "Mouse (Mus musculus)"="Mus musculus"),
                                    selected="Homo sapiens"),
                                selectInput("msigdb_collection", "Collection:",
                                    choices=c(
                                        "Hallmark (H)"="H",
                                        "C2: KEGG (legacy) pathways"="C2_KEGG",
                                        "C2: Reactome pathways"="C2_REACTOME",
                                        "C5: GO Biological Process"="C5_BP",
                                        "C7: Immunologic signatures"="C7"
                                    ),
                                    selected="H"),
                                selectizeInput("msigdb_sets", "Signatures (type to search):",
                                    choices=NULL, multiple=TRUE,
                                    options=list(placeholder='Select one or more gene sets...', maxOptions=2000)),
                                helpText("For Hallmark, leaving this empty loads all 50 sets. For larger collections, select the sets you want.")
                            ),
                            conditionalPanel("input.gs_mode == 'paste'",
                                textInput("paste_name", "Signature name:", placeholder="e.g. MY_SIGNATURE"),
                                textAreaInput("paste_genes", "Gene symbols (one per line or comma-separated):",
                                    rows=4, placeholder="TP53\nBRCA1\nMYC\n..."),
                                actionButton("add_paste", "Add signature", class="btn btn-default btn-block", icon=icon("plus")),
                                actionButton("clear_paste", "Clear all", class="btn btn-link btn-block", style="color:#b94a48;"),
                                strong("Staged signatures:"),
                                verbatimTextOutput("paste_staged")
                            ),
                            hr(),

                            ## ── 3. Scoring method ──
                            h4('3. Scoring method', style="color:#0F344C;"),
                            selectInput("score_method", label=NULL,
                                choices=c("GSVA"="gsva", "ssGSEA"="ssgsea", "AUCell"="aucell", "singscore"="singscore"),
                                selected="gsva"),
                            conditionalPanel("input.score_method == 'gsva'",
                                selectInput("gsva_kcdf", "Kernel (kcdf):",
                                    choices=c("Gaussian (continuous / log-scale)"="Gaussian",
                                              "Poisson (integer counts)"="Poisson"),
                                    selected="Gaussian")
                            ),
                            conditionalPanel("input.score_method == 'ssgsea'",
                                materialSwitch("ssgsea_norm", label="Normalize scores",
                                    value=TRUE, right=TRUE, status='info')
                            ),
                            conditionalPanel("input.score_method == 'aucell'",
                                sliderInput("aucell_thresh", "AUC threshold (top fraction of ranking):",
                                    min=0.01, max=0.50, step=0.01, value=0.05)
                            ),
                            conditionalPanel("input.score_method == 'singscore'",
                                radioButtons("singscore_dir", "Signature direction:", inline=TRUE,
                                    choices=c("Up"="up", "Down"="down", "Both"="both"), selected="up"),
                                helpText("Up: genes expected high. Down: genes expected low. Both: undirected (genes coordinately extreme).")
                            ),
                            hr(),

                            ## ── Run ──
                            actionButton("run_score", "Run Scoring",
                                class="btn btn-info btn-block", icon=icon("play")),
                            uiOutput("score_status")
                        )
                    ),
                    column(9,
                        tabsetPanel(id='ResultTabs',
                            tabPanel(title='Expression Matrix', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('matrix_table')
                                ),
                                div(style="margin-top:16px; text-align:center; padding-bottom:20px;",
                                    downloadButton('download_matrix', 'Download (CSV)')
                                )
                            ),
                            tabPanel(title='Sample Metadata', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('metadata_table')
                                ),
                                div(style="margin-top:16px; text-align:center; padding-bottom:20px;",
                                    downloadButton('download_metadata', 'Download (CSV)')
                                )
                            ),
                            tabPanel(title='Loaded Signatures', hr(),
                                p("Genes are matched against the genes in your expression matrix; signatures with fewer than 2 matched genes are skipped at scoring time.",
                                  style="font-size:12px; color:#777;"),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('geneset_preview')
                                ),
                                div(style="margin-top:16px; text-align:center; padding-bottom:20px;",
                                    downloadButton('download_signatures', 'Download (CSV)')
                                )
                            ),
                            tabPanel(title='Score Matrix', hr(),
                                p(em("Signatures (rows) × samples (columns) — the matrix consumed by downstream tools.")),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('score_matrix_table')
                                ),
                                div(style="margin-top:20px; padding-bottom:30px;",
                                    fluidRow(
                                        column(5, offset=1,
                                            radioButtons("dl_orientation", "Orientation:",
                                                choices=c("Signatures × samples (downstream-ready)"="sig_by_sample",
                                                          "Samples × signatures"="sample_by_sig"),
                                                selected="sig_by_sample")
                                        ),
                                        column(4,
                                            selectInput("dl_format", "File format:",
                                                choices=c("CSV"="csv", "TSV"="tsv"), selected="csv")
                                        )
                                    ),
                                    div(align="center",
                                        helpText("Signatures × samples (signature name in column 1) matches the expression-matrix schema read by the MaGIC Survival and Correlation tools."),
                                        downloadButton('download_scores', 'Download Score Matrix', class="btn btn-info")
                                    )
                                )
                            )
                        )
                    )
                )
            ),


        ## Score Heatmap Page (hidden until scored)
##########################################################################################################################################################
            tabPanel('Score Heatmap',
                fluidRow(
                    column(3,
                        wellPanel(

                            ## Scaling
                            h5(strong("Scaling"), style="color:#0F344C; margin-top:4px;"),
                            hr(),
                            radioButtons("hm_scale", label=NULL,
                                choices=c("Row z-score (per signature)"="row",
                                          "Column z-score (per sample)"="col",
                                          "Raw scores"="none"),
                                selected="row"),

                            ## Color
                            materialSwitch("show_color", label="Color Options", value=TRUE, right=TRUE, status='info'),
                            conditionalPanel("input.show_color",
                                hr(),
                                selectInput("hm_palette", "Color palette:",
                                    choices=c("RdBu (diverging)"="RdBu",
                                              "viridis"="viridis",
                                              "magma"="magma"),
                                    selected="RdBu"),
                                materialSwitch("hm_palette_reverse", label="Reverse palette",
                                    value=FALSE, right=TRUE, status='warning')
                            ),

                            ## Annotation
                            materialSwitch("show_anno", label="Annotation Tracks", value=TRUE, right=TRUE, status='info'),
                            conditionalPanel("input.show_anno",
                                hr(),
                                uiOutput('hm_anno_col_selector'),
                                sliderInput("anno_bar_size", "Bar size (mm):", min=1, max=30, step=1, value=5),
                                sliderInput("anno_font_size", "Label font size (pt):", min=4, max=20, step=1, value=10),
                                uiOutput('anno_color_ui')
                            ),

                            ## Clustering / Dendrograms
                            materialSwitch("show_clustering", label="Clustering & Dendrograms", value=TRUE, right=TRUE, status='info'),
                            conditionalPanel("input.show_clustering",
                                hr(),
                                materialSwitch("cluster_rows", label="Cluster rows (signatures)",
                                    value=TRUE, right=TRUE, status='info'),
                                materialSwitch("cluster_cols", label="Cluster columns (samples)",
                                    value=TRUE, right=TRUE, status='info'),
                                materialSwitch("show_dend", label="Show dendrograms",
                                    value=TRUE, right=TRUE, status='info')
                            ),

                            ## Fonts & Labels
                            materialSwitch("show_fonts", label="Fonts & Labels", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_fonts",
                                hr(),
                                materialSwitch("show_row_names", label="Show row names (signatures)",
                                    value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.show_row_names",
                                    sliderInput("row_font_size", "Row label size:", min=4, max=20, step=1, value=9),
                                    materialSwitch("wrap_row_names", label="Wrap long row labels",
                                        value=TRUE, right=TRUE, status='info'),
                                    conditionalPanel("input.wrap_row_names",
                                        sliderInput("wrap_row_width", "Row wrap width (chars):", min=8, max=60, step=2, value=24)
                                    )
                                ),
                                materialSwitch("show_col_names", label="Show column names (samples)",
                                    value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.show_col_names",
                                    sliderInput("col_font_size", "Column label size:", min=4, max=20, step=1, value=8),
                                    sliderInput("col_font_angle", "Column label angle:", min=0, max=360, step=5, value=45),
                                    materialSwitch("wrap_col_names", label="Wrap long column labels",
                                        value=TRUE, right=TRUE, status='info'),
                                    conditionalPanel("input.wrap_col_names",
                                        sliderInput("wrap_col_width", "Column wrap width (chars):", min=8, max=60, step=2, value=16)
                                    )
                                )
                            ),

                            ## Legend
                            materialSwitch("show_legend", label="Legend", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_legend",
                                hr(),
                                radioButtons("hm_legend_pos", "Position:", inline=TRUE,
                                    choices=c("Right"="right", "Left"="left", "Top"="top", "Bottom"="bottom"),
                                    selected="right"),
                                sliderInput("hm_legend_title_size", "Legend title size:", min=6, max=24, step=1, value=12),
                                sliderInput("hm_legend_label_size", "Legend label size:", min=6, max=24, step=1, value=10)
                            ),

                            ## Resize
                            materialSwitch("show_resize", label="Resize Plot", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_resize",
                                hr(),
                                sliderInput("hm_height", "Plot height (px):", min=200, max=2000, step=50, value=600),
                                sliderInput("hm_width",  "Plot width (px):",  min=200, max=2000, step=50, value=900)
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='HeatmapTabs',
                            tabPanel(title='Heatmap', hr(),
                                fluidRow(style="margin: 0 8px 4px 0;",
                                    column(12, align="right",
                                        actionButton("show_hm_code", label=NULL,
                                            icon=icon("file-code"),
                                            title="View R code to reproduce this heatmap",
                                            class="btn btn-default btn-sm",
                                            style="border-radius:6px; font-size:16px; padding:4px 8px;"
                                        )
                                    )
                                ),
                                hr(),
                                div(style="overflow-x:auto; width:100%;",
                                    withSpinner(type=6, color='#5bc0de',
                                        plotOutput("score_heatmap_out", height='100%')
                                    )
                                ),
                                div(style="margin-top:30px; text-align:center; padding-bottom:50px;",
                                    div(style="display:inline-block; width:250px; margin-bottom:10px;",
                                        selectInput("hm_download_format", "Download format:",
                                            choices=c('png','pdf','svg','tiff','jpeg','eps'))
                                    ),
                                    br(),
                                    downloadButton('download_heatmap', 'Download Heatmap')
                                )
                            )
                        )
                    )
                )
            ),


        ## Per-signature Distributions Page (hidden until scored)
##########################################################################################################################################################
            tabPanel('Distributions',
                fluidRow(
                    column(3,
                        wellPanel(
                            h4('Distribution Options', align='center', style="color:#0F344C;"),
                            hr(),
                            selectInput("dist_signature", "Signature:", choices=NULL),
                            uiOutput("dist_group_selector"),
                            radioButtons("dist_plot_type", "Plot type:", inline=TRUE,
                                choices=c("Boxplot"="boxplot", "Violin"="violin", "Violin + box"="violin_box"),
                                selected="boxplot"),
                            materialSwitch("dist_jitter", label="Show data points", value=TRUE, right=TRUE, status='info'),

                            materialSwitch("dist_show_stats", label="Statistics", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.dist_show_stats",
                                hr(),
                                materialSwitch("dist_pairwise", label="Pairwise comparisons", value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.dist_pairwise",
                                    selectInput("dist_pairwise_test", "Pairwise test:",
                                        choices=c("Wilcoxon"="wilcox.test", "t-test"="t.test"),
                                        selected="wilcox.test")
                                ),
                                materialSwitch("dist_global", label="Global test (Kruskal/ANOVA)", value=FALSE, right=TRUE, status='info'),
                                conditionalPanel("input.dist_global",
                                    selectInput("dist_global_test", "Global test:",
                                        choices=c("Kruskal-Wallis"="kruskal.test", "ANOVA"="anova"),
                                        selected="kruskal.test")
                                )
                            ),

                            materialSwitch("dist_show_colors", label="Group Colors", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.dist_show_colors",
                                hr(),
                                selectInput("dist_palette", "Palette:",
                                    choices=c("Default"="default", "Set1"="Set1", "Set2"="Set2",
                                              "viridis"="viridis", "Custom"="custom"),
                                    selected="default"),
                                conditionalPanel("input.dist_palette == 'custom'",
                                    uiOutput("dist_color_ui")
                                )
                            ),

                            materialSwitch("dist_show_fonts", label="Fonts", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.dist_show_fonts",
                                hr(),
                                sliderInput("dist_title_font", "Title size:", min=6, max=28, step=1, value=14),
                                sliderInput("dist_axis_font", "Axis label size:", min=6, max=24, step=1, value=12),
                                sliderInput("dist_tick_font", "Tick label size:", min=6, max=20, step=1, value=10),
                                radioButtons("dist_legend_pos", "Legend position:", inline=TRUE,
                                    choices=c("Right"="right", "Left"="left", "Top"="top", "Bottom"="bottom"),
                                    selected="right"),
                                sliderInput("dist_legend_title_size", "Legend title size:", min=6, max=24, step=1, value=12),
                                sliderInput("dist_legend_label_size", "Legend label size:", min=6, max=24, step=1, value=10),
                                sliderInput("dist_x_angle", "X label angle:", min=0, max=90, step=5, value=45)
                            ),

                            materialSwitch("dist_show_resize", label="Resize Plot", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.dist_show_resize",
                                hr(),
                                sliderInput("dist_height", "Plot height (px):", min=200, max=1600, step=50, value=550),
                                sliderInput("dist_width",  "Plot width (px):",  min=200, max=1600, step=50, value=700)
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='DistTabs',
                            tabPanel(title='Distribution', hr(),
                                fluidRow(style="margin: 0 8px 4px 0;",
                                    column(12, align="right",
                                        actionButton("show_dist_code", label=NULL,
                                            icon=icon("file-code"),
                                            title="View R code to reproduce this plot",
                                            class="btn btn-default btn-sm",
                                            style="border-radius:6px; font-size:16px; padding:4px 8px;"
                                        )
                                    )
                                ),
                                hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    plotOutput("dist_out", height='100%')
                                ),
                                div(style="margin-top:30px; text-align:center; padding-bottom:50px;",
                                    div(style="display:inline-block; width:250px; margin-bottom:10px;",
                                        selectInput("dist_download_format", "Download format:",
                                            choices=c('png','pdf','svg','tiff','jpeg','eps'))
                                    ),
                                    br(),
                                    downloadButton('download_dist', 'Download Plot')
                                )
                            )
                        )
                    )
                )
            ),


        ## Footer
##########################################################################################################################################################
            tags$footer(
                wellPanel(
                    fluidRow(
                        column(4, align='center',
                        tags$a(href="https://github.com/MaGIC-Analytics/magic-signature", icon("github", "fa-3x")),
                        tags$h4('GitHub to submit issues/requests')
                        ),
                        column(4, align='center',
                        tags$a(href="http://www.bioinformagic.io/", icon("magic", "fa-3x")),
                        tags$h4('MaGIC Home Page')
                        ),
                        column(4, align='center',
                        tags$a(href="https://github.com/MaGIC-Analytics", icon("address-card", "fa-3x")),
                        tags$h4("Developer's Page")
                        )
                    ),
                    fluidRow(
                        column(12, align='center',
                            HTML('<a href="https://www.youtube.com/watch?v=dQw4w9WgXcQ">
                            <p>&copy;
                                <script language="javascript" type="text/javascript">
                                var today = new Date()
                                var year = today.getFullYear()
                                document.write(year)
                                </script>
                            </p>
                            </a>
                            ')
                        )
                    )
                )
            )
        )# Ends navbarPage
    )# Ends fluidPage
)# Ends tagList
