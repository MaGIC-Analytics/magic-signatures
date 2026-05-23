# ─── Tab Visibility Management ──────────────────────────────────────────────────
# The Score Heatmap and Distributions tabs only appear once a score matrix has
# been computed on the Data & Scoring tab.

# Hide score-dependent tabs on initial load (runs once — no reactive deps).
observe({
    hideTab(inputId="NAVTABS", target="Score Heatmap")
    hideTab(inputId="NAVTABS", target="Distributions")
})

# Reveal them (and jump to the heatmap) once scoring succeeds.
.scores_revealed <- reactiveVal(FALSE)

observeEvent(ScoreMatrix(), {
    sm <- ScoreMatrix()
    if (!is.null(sm) && nrow(sm) > 0) {
        showTab(inputId="NAVTABS", target="Score Heatmap")
        showTab(inputId="NAVTABS", target="Distributions")
        if (!isTRUE(.scores_revealed())) {
            .scores_revealed(TRUE)
            updateTabsetPanel(session, inputId="NAVTABS", selected="Score Heatmap")
            shinyjs::delay(300, shinyjs::runjs("$(window).trigger('resize');"))
        }
    }
}, ignoreNULL=TRUE)
