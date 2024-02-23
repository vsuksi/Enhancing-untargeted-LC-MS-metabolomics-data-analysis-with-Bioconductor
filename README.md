# Enhancing untargeted LC-MS metabolomics data analysis with Bioconductor

## Instructions
1. Download the repository

2. Make sure you have Quarto and TinyTex installed (do "quarto install tinytex", does not add to PATH). You may also need to have RTools installed to compile packages for which binary distributions are not available anymore.

3. Navigate to /Thesis, launch R to bootstrap renv and then do "renv::restore()" to create the environment

- To render the thesis, navigate to /Thesis and do "quarto render --to pdf"
- To render the standalone example analysis, navigate to /Thesis and do "quarto render example_analysis.qmd"
- To test visualizations or other functions, navigate to /Thesis, launch R and do "source("Input/Code/x")" to load the functions. Then do "readRDS("Input/Code/x")" to load tse instances for testing.


## Specs
- R version 4.3.2 (Bioconductor 3.18)
- Quarto version 1.3.450
- Pandoc version 3.1.1
- Dart Sass version 1.55
- knitr version 1.45
- rmarkdown version 2.25
