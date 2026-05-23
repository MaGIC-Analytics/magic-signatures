FROM rocker/shiny-verse:4.5.3
LABEL authors="Alex Lemenze" \
    description="Docker image for MaGIC Signature Scoring Tool"

# ── System dependencies ────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    sudo \
    libhdf5-dev \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libv8-dev \
    libsodium-dev \
    libglpk40 \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libmagick++-dev \
    cmake && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ── CRAN packages (ggpubr excluded — installed after Bioconductor) ──────────
RUN R -e "install.packages(c( \
    'BiocManager', \
    'shiny', \
    'shinyjs', \
    'shinythemes', \
    'shinycssloaders', \
    'shinyWidgets', \
    'DT', \
    'tidyverse', \
    'data.table', \
    'RColorBrewer', \
    'colourpicker', \
    'circlize', \
    'msigdbr', \
    'viridis' \
    ), repos='https://cran.rstudio.com/', dependencies=TRUE)"

# ── Bioconductor: heatmap + scoring engines ─────────────────────────────────
RUN R -e "BiocManager::install(c( \
    'ComplexHeatmap', \
    'BiocGenerics', \
    'S4Vectors', \
    'IRanges', \
    'GSEABase', \
    'GSVA', \
    'AUCell', \
    'singscore' \
    ), ask=FALSE, update=FALSE)"

# ── ggpubr after Bioconductor (depends on rstatix which needs BiocManager) ──
RUN R -e "install.packages('ggpubr', repos='https://cran.rstudio.com/', dependencies=TRUE)"

# ── Copy application files ───────────────────────────────────────────────────
COPY ./app /srv/shiny-server/
COPY shiny-customized.config /etc/shiny-server/shiny-server.conf

# ── Permissions ──────────────────────────────────────────────────────────────
RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 8080
USER shiny
CMD ["/usr/bin/shiny-server"]
