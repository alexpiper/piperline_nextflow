#!/bin/bash
#SBATCH --job-name=nf_dada
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB 
#SBATCH --time=640:00:00
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au
#SBATCH --mail-type=ALL
#SBATCH --account=pathogens
#SBATCH --export=none
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

# SETUP 
# install the latest nextflow within a virtual enviornemnt
#module purge
#module load Python/3.8.2-GCCcore-9.3.0
#virtualenv ~/nf
#source ~/nf/bin/activate
#Install latest release of nextflow
#wget -O ~/nf/nextflow --no-check-certificate --content-disposition https://github.com/nextflow-io/nextflow/releases/download/v21.04.0-edge/nextflow

module purge
module load Nextflow/20.10.0

# Run COI data

# Data from /group/pathogens/Alexp/Metabarcoding/imappests/data/JDYG3
cp /group/pathogens/Alexp/Metabarcoding/imappests/data/JDYG3/* /group/pathogens/Alexp/Metabarcoding/piperline_testing/imap_tests/COI

for d in ./*/ ; do
(cd "$d" && mv *.fastq.gz ../. );
done

cd /group/pathogens/Alexp/Metabarcoding/piperline_testing/imap_tests/COI
nextflow pull alexpiper/piperline
nextflow run alexpiper/piperline --reads '*_R{1,2}_001.fastq.gz' --lengthvar false \
--fwdprimer 'GGDACWGGWTGAACWGTWTAYCCHCC' --fwdprimer_name 'fwhF2' --revprimer 'GTRATWGCHCCDGCTARWACWGG' --revprimer_name 'fwhR2n' \
--reference 'idtaxa_bftrimmed.rds' --samplesheet 'SampleSheet_JDYG3.csv' --runparams 'runParameters.xml' --phmm 'folmer_fullength_model.rds' --blast_db 'merged_final.fa.gz' --species_db 'merged_arthropoda_rdp_species.fa.gz' \
-profile basc --subsample 10000 -resume

nextflow run alexpiper/piperline --reads '*_R{1,2}_001.fastq.gz' --lengthvar false \
--fwdprimer 'GGDACWGGWTGAACWGTWTAYCCHCC' --fwdprimer_name 'fwhF2' --revprimer 'GTRATWGCHCCDGCTARWACWGG' --revprimer_name 'fwhR2n' \
--reference 'merged_final.fa.gz' --samplesheet 'SampleSheet_JDYG3.csv' --runparams 'runParameters.xml' --phmm 'folmer_fullength_model.rds' --blast_db 'merged_final.fa.gz' --species_db 'merged_arthropoda_rdp_species.fa.gz' \
-profile basc --subsample 10000 --taxassignment 'rdp'



taxassignment
# -bg runs nextflow in background to ensure conneciton issues dont stop background


# Test with multiplexed ITS data
cd /group/pathogens/Alexp/Metabarcoding/piperline_testing/imap_tests/ITS
nextflow pull alexpiper/piperline
nextflow run alexpiper/piperline --reads '*_R{1,2}_001.fastq.gz' --lengthvar true \
--fwdprimer 'GGTCATTTAGAGGAAGTAA;ARTCTTTGAACGCACATTG;GTGAAATTGYTRAAAGGGAA;CARGAYATGATYAACGAGG;GTTCGAGGCTGGTATCTCC' --fwdprimer_name 'ITS1Fngs;ITS7ngsOphio;LF402F;CALngsF1;EF1087Fngs' \
--revprimer 'TTYRCKRCGTTCTTCATCG;CCTSCSCTTANTDATATGC;CGATCGATTTGCACGTCAG;GRATCATCTCRTCVACCTC;GTCTCRATACGGCCRAC' --revprimer_name 'ITS2ngs;ITS4ngsUni;LR5Fungngs;CALngsR1;EF646Rngs' \
--reference '/group/pathogens/Alexp/Metabarcoding/piperline_testing/imap_tests/ITS/reference/CTedited_UNITE_RefS_Aug2021.rds' --samplesheet 'SampleSheet.csv' --runparams 'runParameters.xml' --coding false \
-profile basc --subsample 10000 -resume



