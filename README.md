# Piperline - an automated DADA2 & NextFlow based pipeline for pest and pathogen metabarcoding

This pipeline is adapted from https://github.com/HPCBio/dada2-Nextflow

## Prerequisites

Software:
Nextflow (>=20.11.0)
R (>= 3.2.0)

R packages:
dada2 (>= 1.8), Rcpp (>= 0.11.2), methods (>= 3.2.0), DECIPHER, phangorn


## Run on AgVic BASC Cluster:
module load Nextflow/20.10.0

nextflow pull alexpiper/piperline
nextflow run alexpiper/piperline --reads '*_{1,2}.fastq.gz' --trimFor 24 --trimRev 25 --lengthvar false \
--fwdprimer 'ACCTGCGGARGGATCA' --revprimer 'GAGATCCRTTGYTRAAAGTT' --reference 'Fungal_LSU_v11_March2018.RData' \
-profile basc -r main


## Resources
* [The indexing system we use](https://alexpiper.github.io/iMapPESTS/indexing.html)
* [The DADA2 algorithm](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)
* [The benefits of amplicon sequence variants over OTU's](https://www.nature.com/articles/ismej2017119)
* [Taxonomic assignment using IDTAXA](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0521-5)


## Credits

Alexander Piper is the main developer of this pipeline

This pipeline was modified from the 16s DADA2 Nextflow workflow developed by Chris Fields (https://github.com/HPCBio/dada2-Nextflow).

## Licence

This project is licensed under the MIT License - see the LICENSE.md file for details

