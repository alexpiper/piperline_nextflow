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


module load Nextflow/20.10.0
module load charliecloud/0.23-GCCcore-9.3.0
module load Java/12.0.1

# Run a test set using the dada2 test samples
# Get example samples
#wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip -O miseqsopdata.zip 
#unzip miseqsopdata.zip

# Get reference sequences
#wget https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1 -O gg_13_8_train_set_97.fa.gz



# Run COI data

# Data from /group/pathogens/Alexp/Metabarcoding/imappests/data/JDYG3
cp /group/pathogens/Alexp/Metabarcoding/imappests/data/JDYG3/* /group/pathogens/Alexp/Metabarcoding/test/COI

cd /group/pathogens/Alexp/Metabarcoding/test/COI
nextflow pull alexpiper/piperline
nextflow run alexpiper/piperline --reads '*_R{1,2}_001.fastq.gz' --lengthvar false \
--fwdprimer 'GGDACWGGWTGAACWGTWTAYCCHCC' --fwdprimer_name 'fwhF2' --revprimer 'GTRATWGCHCCDGCTARWACWGG' --revprimer_name 'fwhR2n' \
--reference 'idtaxa_bftrimmed.rds' \
-profile basc --subsample true

# Test with lenthvar false

for d in ./*/ ; do
(cd "$d" && mv *.fastq.gz ../. );
done

# Run dev version

cd /group/pathogens/Alexp/Metabarcoding/test/COI
nextflow run alexpiper/piperline --reads '*_R{1,2}_001.fastq.gz' --lengthvar false \
--fwdprimer 'GGDACWGGWTGAACWGTWTAYCCHCC' --fwdprimer_name 'fwhF2' --revprimer 'GTRATWGCHCCDGCTARWACWGG' --revprimer_name 'fwhR2n' \
--reference 'idtaxa_bftrimmed.rds' \
-profile basc --subsample true -resume