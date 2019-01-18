#!/bin/bash
#Run the R markdown script to generate a pdf with the plots resulting from the analysis
R -e "rmarkdown::render('ase_analysis_north_america.Rmd', output_file = 'results/ase_analysis_north_america.pdf')"
R -e "rmarkdown::render('ase_analysis_south_america.Rmd', output_file = 'results/ase_analysis_south_america.pdf')"

