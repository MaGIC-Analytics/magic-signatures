# magic-signature

![GitHub last commit](https://img.shields.io/github/last-commit/MaGIC-Analytics/magic-signature)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![made with Shiny](https://img.shields.io/badge/R-Shiny-blue)](https://shiny.rstudio.com/)

The **MaGIC Signature Scoring Tool** computes per-sample gene-set signature scores
(GSVA, ssGSEA, AUCell, or singscore) from a normalized expression matrix, producing a
signature × sample score matrix that feeds directly into the MaGIC Survival and
Correlation tools.

## Running the App
This Shiny app is packaged in a Docker container for easy deployment. Build the image
yourself (and customize ports as needed):
```
docker build -t signature_service .
docker run -d --rm -p 8080:8080 signature_service
```
It will be hosted at http://localhost:8080.

## What it does
- **Inputs:** a normalized expression matrix (genes × samples, CSV/TSV, log-scale preferred —
  the same schema as the MaGIC QC tool's output) plus sample metadata for grouping.
- **Gene sets**, from any of three sources:
  - **GMT upload** (standard Broad format).
  - **MSigDB** collections (Hallmark, C2:CP KEGG/Reactome, C5:GO BP, C7) for human or mouse via `msigdbr`.
  - **Paste-in** custom signatures (name a signature, paste a gene list, repeat).
- **Scoring methods:**
  - **GSVA** — cohort-relative, stable comparative analysis (kcdf Gaussian/Poisson).
  - **ssGSEA** — single-sample enrichment (optional normalization).
  - **AUCell** — rank-based, robust to dropouts (AUC threshold).
  - **singscore** — directional single-sample scoring (up / down / both).
- **Outputs:** a signature × sample score matrix (CSV/TSV) in the orientation the Survival and
  Correlation tools read directly, plus an annotated score heatmap and per-signature
  distribution plots — all downloadable.

## Output schema (pipeline contract)
The canonical export is **signatures × samples**: the signature name occupies the first
column and every other column is a sample. This matches the feature × sample expression-matrix
schema read by `survival_service` and `correlation_service`, so each signature can be used
downstream exactly like a gene. A samples × signatures orientation is also offered for joining
with metadata tables.

## Building the App
Built from `magic-modules-template`. To extend it:
- Modify the Dockerfile to include the correct library installations.
- Modify the UI to load the proper libraries and set the tool name.
- Get Google Analytics for the tool (if desired) and add it to `app/www/GA.html`.
- Push to Cloud Run following the deployment instructions, or set up Cloud Run to pull from the repo on commit.
