# Spatial profiling identifies a distinct and topographically-defined tumor microenvironment that emerges during multiple myeloma evolution — Analysis and Figure Code

This repository contains the R code used for the analysis and figures in
"Spatial profiling identifies a distinct and topographically-defined tumor microenvironment that emerges during multiple myeloma evolution"
(BCD, 2026). It covers processing of spatial
transcriptomics (Xenium) and imaging mass cytometry (IMC) data, through to
the plots and statistics reported in the paper.

## Repository contents

- `Analysis_Script_Paper.R` — main analysis pipeline: reading in raw
  Xenium/IMC output, QC, merging, Harmony integration, clustering,
  annotation, spatial region calling (hexagonal grid / tumor architecture),
  cell–cell interaction analysis, and pseudobulk DEG testing. Produces the
  processed `SpatialExperiment` / Seurat objects used by the figure script.
- `Figure_configuration_Paper.R` — generates all main and supplementary
  figures (UMAPs, spatial plots, boxplots, survival analysis, etc.) from
  the processed objects produced by `Analysis_Script_Paper.R`.


Intended order: run `Analysis_Script_Paper.R` first to go from raw data to
processed objects, then `Figure_configuration_Paper.R` to produce the
figures from those objects.


## Data availability

`Analysis_Script_Paper.R` expects raw Xenium/IMC instrument output as
input; `Figure_configuration_Paper.R` expects the processed `.rds` objects
produced by the analysis script. 

- Processed Xenium ST and IMC files can be downloaded from the Mendeley portal including 
  relevant metadata (https://data.mendeley.com/preview/8bcwmrxp6w) while unprocessed Xenium 
  and IMC files will be made available upon request to the corresponding authors. 
- File and folder paths in `Analysis_Script_Paper.R` have been
  genericized (e.g. `sample_dir_1`) — they do not reflect the actual raw
  data directory structure used internally.

## Requirements

- R version >= 4.X
- Key packages: `Seurat`, `imcRtools`, `SpatialExperiment`, `scran`,
  `scater`, `cytomapper`, `harmony`, `sf`, `data.table`, `arrow`,
  `survival`, `survminer`, and others listed in the `pkgs` vector at the
  top of each script.

To reproduce the computational environment:

```r
# Option 1: manual install
install.packages(c("Seurat", "ggplot2", "dplyr", ...))
BiocManager::install(c("imcRtools", "SpatialExperiment", "scran", "scater", "cytomapper"))


```

## Usage

1. **Analysis** (`Analysis_Script_Paper.R`): update the data directory
   variables at the top of the script to point to your local copy of the
   raw Xenium/IMC output, then run section by section (organized by
   processing step, e.g. `# QC ...`, `# Integration ...`,
   `# FindAllMarkers ...`) to produce the processed `.rds` objects.
2. **Figures** (`Figure_configuration_Paper.R`): update the `data_dir` and
   `output_dir` variables at the top to point to the processed objects
   from step 1 and your desired output location, then run section by
   section — organized by figure (e.g. `# Figure 1 ...`,
   `# Sup Figure 3 ...`).
3. Output figures are written as PDFs to `output_dir`; intermediate
   tables/objects from the analysis script are written as `.csv`/`.rds`
   files alongside it.

## Note on reproducibility

These scripts are provided for transparency of the analysis and plotting
code. Because patient-level data cannot be shared publicly, neither script
is fully executable "as-is" without access to the underlying raw or
processed data objects. Data can be requested as described above.
`Analysis_Script_Paper.R` in particular reflects an iterative analysis
process — some blocks represent exploratory steps rather than the final
pipeline; the code that feeds into the final figures is the version
retained/used in `Figure_configuration_Paper.R`.

## License

MIT License — Copyright (c) 2026 Marnix Koops. Free to use, copy, modify,
and redistribute, with attribution and no warranty. See
opensource.org/license/mit for the
full text.

## Citation

If you use this code, please cite:

> Koops M, Bertamini L, Papazian N, Korst C, Budai MJ, Hoogenboezem RM,
Sanders MA, van der Holt B, van Duin M, Broijl A, Balogh P, Zweegman S,
Raaijmakers MHGP, van de Donk NWCJ, Kellermayer Z, Sonneveld P, Cupedo T.
Spatial profiling identifies a distinct and topographically-defined
tumor microenvironment that emerges during multiple myeloma evolution.
2026. [Blood Cancer Discovery, in preparation/under review — DOI will be added upon
publication.]
