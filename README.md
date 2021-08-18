![Github CI Status](https://github.com/grbot/TADA/actions/workflows/ci.yml/badge.svg)
[![Travis-CI Build Status](https://travis-ci.com/h3abionet/TADA.svg?branch=master)](https://travis-ci.com/h3abionet/TADA)
[![DOI](https://zenodo.org/badge/218786496.svg)](https://zenodo.org/badge/latestdoi/218786496)

<p>
<img align="left" src="./assets/cbio_logo.png" width="240" hspace="50"/>
<img align="left" src="./assets/HPCBio.png" width="350" hspace="50"/>
</br></br></br>
</p>

# Piperline - an automated DADA2 based workflow for pest and pathogen metabarcoding

A dada2-based workflow using the Nextflow workflow manager.  The basic pipeline is currently implemented, including some basic read-tracking. 

This pipeline is adapted from https://github.com/HPCBio/dada2-Nextflow

## Prerequisites

Nextflow (>=20.11.0), dada2 (>= 1.8), R (>= 3.2.0), Rcpp (>= 0.11.2), methods (>= 3.2.0), DECIPHER, phangorn, biomformat


Before starting, i recommend learning about the individual componenets of the workflow:
* [The indexing system we use](https://alexpiper.github.io/iMapPESTS/indexing.html)
* [The DADA2 algorithm](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)
* [The benefits of amplicon sequence variants over OTU's](https://www.nature.com/articles/ismej2017119)
* [Taxonomic assignment using IDTAXA](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0521-5)


## Credits

This pipeline was modified from the DADA2 pipeline Nextflow workflow developed by Chris Fields (https://github.com/HPCBio/dada2-Nextflow) fields

