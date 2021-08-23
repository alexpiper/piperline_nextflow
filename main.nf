#!/usr/bin/env nextflow

/*
========================================================================================
                    P I P E R L I N E
========================================================================================
 DADA2 NEXTFLOW PIPELINE FOR BIOSECURITY METABARCODING

----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    ===================================
     ${workflow.repository}/piperline  ~  version ${params.version}
    ===================================
    Usage:

    This pipeline can be run specifying parameters in a config file or with command line flags.
    The typical example for running the pipeline with command line flags is as follows:
    nextflow run alexpiper/piperline --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex

    The typical command for running the pipeline with your own config (instead of command line flags) is as follows:
    nextflow run alexpiper/piperline -c dada2_user_input.config -profile uct_hex
    where:
    dada2_user_input.config is the configuration file (see example 'dada2_user_input.config')
    NB: -profile uct_hex still needs to be specified from the command line

    To override existing values from the command line, please type these parameters:

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. Currently profile available for UCT's HPC 'uct_hex' - create your own if necessary
                                    NB -profile should always be specified on the command line, not in the config file
      --trimFor                     integer. headcrop of read1 (set 0 if no trimming is needed)
      --trimRev                     integer. headcrop of read2 (set 0 if no trimming is needed)
      --reference                   Path to taxonomic database to be used for annotation (e.g. gg_13_8_train_set_97.fa.gz)

    All available read preparation parameters:
      --trimFor                     integer. headcrop of read1
      --trimRev                     integer. headcrop of read2
      --truncFor                    integer. truncate read1 here (i.e. if you want to trim 10bp off the end of a 250bp R1, truncFor should be set to 240). enforced before trimFor/trimRev
      --truncRev                    integer. truncate read2 here ((i.e. if you want to trim 10bp off the end of a 250bp R2, truncRev should be set to 240). enforced before trimFor/trimRev
      --maxEEFor                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
      --maxEERev                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
      --truncQ                      integer. Truncate reads at the first instance of a quality score less than or equal to truncQ; default=2
      --maxN                        integer. Discard reads with more than maxN number of Ns in read; default=0
      --maxLen                      integer. maximum length of trimmed sequence; maxLen is enforced before trimming and truncation; default=Inf (no maximum)
      --minLen                      integer. minLen is enforced after trimming and truncation; default=50
      --rmPhiX                      {"T","F"}. remove PhiX from read
      --minOverlap                  integer. minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 integer. The maximum mismatches allowed in the overlap region; default=0
      --trimOverhang                {"T","F"}. If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off.
                                    "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction                                               primer region. Default="F" (false)

    Other arguments:
      --dadaOpt.XXX                 Set as e.g. --dadaOpt.HOMOPOLYMER_GAP_PENALTY=-1 Global defaults for the dada function, see ?setDadaOpt in R for available options and their defaults
      --pool                        Should sample pooling be used to aid identification of low-abundance ASVs? Options are
                                    pseudo pooling: "pseudo", true: "T", false: "F"
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail adadaRsess to get a summary e-mail with details of the run
                                    sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --idType                      The ASV IDs are renamed to simplify downstream analysis, in particular with downstream tools.  The
                                    default is "ASV", which simply renames the sequences in sequencial order.  Alternatively, this can be
                                    set to "md5" which will run MD5 on the sequence and generate a QIIME2-like unique hash.
      --subsample                   subsample a random 10,000 reads to speed up process of testing

    Help:
      --help                        Will print out summary above when executing nextflow run alexpiper/piperline

    Merging arguments (optional):
      --minOverlap                  The minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 The maximum mismatches allowed in the overlap region; default=0.
      --trimOverhang                If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off. "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen
                                    when reads are longer than the amplicon and read into the other-direction primer region. Default="F" (false)
      --minMergedLen                Minimum length of fragment *after* merging
      --maxMergedLen                Maximum length of fragment *after* merging

    Taxonomic arguments (optional):
      --species                     Specify path to fasta file. See dada2 addSpecies() for more detail.
    """.stripIndent()
}

// TODO: add checks on options

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//Validate inputs
if ( params.trimFor == false && params.lengthvar == false) {
    exit 1, "Must set length of R1 (--trimFor) that needs to be trimmed (set 0 if no trimming is needed)"
}

if ( params.trimRev == false && params.lengthvar == false) {
    exit 1, "Must set length of R2 (--trimRev) that needs to be trimmed (set 0 if no trimming is needed)"
}

// if ( params.reference == false ) {
//     exit 1, "Must set reference database using --reference"
// }

if (params.fwdprimer == false && params.lengthvar == true){
    exit 1, "Must set forward primer using --fwdprimer"
}

if (params.revprimer == false && params.lengthvar == true){
    exit 1, "Must set reverse primer using --revprimer"
}

if (!(['simple','md5'].contains(params.idType))) {
    exit 1, "--idType can only be set to 'simple' or 'md5', got ${params.idType}"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


// Create main channels
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { samples_ch }

// Header log info
log.info "==================================="
log.info " ${params.base}/piperline  ~  version ${params.version}"
log.info "==================================="
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Forward primer'] = params.fwdprimer
summary['Reverse primer'] = params.revprimer
summary['Length variable']= params.lengthvar
summary['trimFor']        = params.trimFor
summary['trimRev']        = params.trimRev
summary['truncFor']       = params.truncFor
summary['truncRev']       = params.truncRev
summary['truncQ']         = params.truncQ
summary['maxEEFor']       = params.maxEEFor
summary['maxEERev']       = params.maxEERev
summary['maxN']           = params.maxN
summary['maxLen']         = params.maxLen
summary['minLen']         = params.minLen
summary['rmPhiX']         = params.rmPhiX
summary['minOverlap']     = params.minOverlap
summary['maxMismatch']    = params.maxMismatch
summary['trimOverhang']   = params.trimOverhang
summary['species']        = params.species
summary['dadaOpt']        = params.dadaOpt
summary['pool']           = params.pool
summary['qualityBinning'] = params.qualityBinning
summary['Reference']      = params.reference
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) {
    summary['E-mail AdadaRsess'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 *
 * Step 0: Copy read files and optionally subsample for testing
 *
 */

if (params.subsample == true) {
    process setup_subsample {
        tag { "subsample.${sample_id}" }

        input:
        tuple sample_id, file(reads) from samples_ch

        output:
        tuple sample_id, file("data/${sample_id}_R1_001.fastq.gz"), file("data/${sample_id}_R2_001.fastq.gz") into samples_toqual_ch
        tuple sample_id, file("data/*R[12]_001.fastq.gz") into samples_tofilt_ch

        """
        #!/bin/bash
        mkdir data
        seqtk sample -s100 ${reads[0]} 10000 | pigz -p ${task.cpus} > data/${sample_id}_R1_001.fastq.gz
        seqtk sample -s100 ${reads[1]} 10000 | pigz -p ${task.cpus} > data/${sample_id}_R2_001.fastq.gz
        """
    }
} else if (params.subsample == false){
    process setup {
        tag { "setup.${sample_id}" }

        input:
        tuple sample_id, file(reads) from samples_ch

        output:
        tuple sample_id, file("data/${sample_id}_R1_001.fastq.gz"), file("data/${sample_id}_R2_001.fastq.gz") into samples_toqual_ch
        tuple sample_id, file("data/*R[12]_001.fastq.gz") into samples_tofilt_ch

        """
        #!/bin/bash
        mkdir data
        cp ${reads[0]} data/${sample_id}_R1_001.fastq.gz
        cp ${reads[1]} data/${sample_id}_R2_001.fastq.gz
        """
    }
}

/*
 *
 * Step 1: Pre-filter Quality control
 *
 */

process runFastQC {
    tag { "rFQC.${sample_id}" }
    publishDir "${params.outdir}/qc/FASTQC-prefilter", mode: "copy", overwrite: true

    input:
    tuple sample_id, file(For), file(Rev) from samples_toqual_ch

    output:
    file '*_fastqc.{zip,html}' into fastqc_files_ch, fastqc_files2_ch

    """
    fastqc --nogroup -q ${For} ${Rev}
    """
}


// TODO: combine MultiQC reports and split by directory (no need for two)
process runMultiQC {
    tag { "runMultiQC" }
    publishDir "${params.outdir}/qc/MultiQC-prefilter", mode: 'copy', overwrite: true

    input:
    file('./raw-seq/*') from fastqc_files_ch.collect()

    output:
    file "*_report.html" into multiqc_report
    file "*_data"

    script:
    interactivePlots = params.interactiveMultiQC == true ? "-ip" : ""
    """
    multiqc ${interactivePlots} .
    """
}

// Summarise indexes used for each sample
//process summarise_index {
//    tag { "summarise_index_${sample_id}" }
//
//    input:
//    set sample_id, file(reads) from samples_tofilt_ch
//
//    output:
//    file '*_indexes.txt' into index_files
//
//    when:
//    params.precheck == false
//
//    script:
//    """
//  #!/bin/bash
//  zcat "${reads[0]}" | grep '^@M' | rev | cut -d':' -f 1 | rev | sort | uniq -c | sort -nr  | sed 's/+/ /' | sed 's/^ *//g' > ${sample_id}_indexes.txt
//  done
//    """
//}

// Calculate switch rate
//process index_calc {
//    tag { "index_calc_${sample_id}" }
//
//    input:
//    file('*_indexes.txt') from index_files.collect()
//
//    output:
//    file "index_switch_calc.txt" into index_switch
//
//    when:
//    params.precheck == false
//
//    script:
//    """
//  #!/bin/bash
//  ls | grep '_indexes.txt' | sort > files
//  grep -v 'Undetermined' files | xargs cat > determined_counts.txt
//
//  # Get all potential switched combinations of used indexes
//  index1=$(cat determined_counts.txt | cut -d' ' -f 2)
//  index2=$(cat determined_counts.txt | cut -d' ' -f 3)
//
//  [ -e all_combinations.txt ] && rm all_combinations.txt
//  touch all_combinations.txt
//  for i in ${index1}
//  do
//    for j in ${index2}
//    do
//      if [ "$i" \< "$j" ]
//      then
//       echo $i $j >> all_combinations.txt
//      fi
//    done
//  done
//
//  # Count number of undetermined reads
//  grep 'Undetermined' files | xargs cat > undetermined_counts.txt
//  cat undetermined_counts.txt | cut -d' ' -f 2,3 > undetermined_index.txt
//
//  # Count number of correctly demultiplexed reads
//  correct_counts=$(cat determined_counts.txt | cut -d' ' -f 1 | awk '{ SUM += $1} END { print SUM }')
//
//  # Count number of switched reads
//  comm -12 <(sort all_combinations.txt) <(sort undetermined_index.txt) > switched_indexes.txt
//  switched_counts=$(grep -f "switched_indexes.txt" "undetermined_counts.txt" | cut -d' ' -f 1 | awk '{ SUM += $1} END { print SUM }')
//
//  # Count number of other reads (these can be sequencing errors, PhiX and other junk)
//  other_counts=$(grep -v -f "switched_indexes.txt" "undetermined_counts.txt" | cut -d' ' -f 1 | awk '{ SUM += $1} END { print SUM }')
//
//  # Calculate switch rate (in percentage)
//  calc(){ awk "BEGIN { print "$*" }"; }
//  switch_rate=$(calc $switched_counts/$correct_counts)
//  switch_rate_perc=$(calc $switched_counts/$correct_counts*100)
//
//  # Print results to file
//  touch index_switch_calc.txt
//  echo "Correctly demultiplexed reads: ${correct_counts}" >> index_switch_calc.txt
//  echo "Switched reads: ${switched_counts}" >> index_switch_calc.txt
//  echo "Other undetermined reads: ${other_counts}" >> index_switch_calc.txt
//  echo "Index switching rate: ${switch_rate} (${switch_rate_perc}%)" >> index_switch_calc.txt
//    """
//}


/*
 *
 * Step 2: Filter N bases
 *
 */
process Nfilter {
    tag { "nfilter_${sample_id}" }

    input:
    set sample_id, reads from samples_tofilt_ch

    output:
    set val(sample_id), "${sample_id}.R[12].noN.fastq.gz" optional true into filt_step2
    set val(sample_id), "${sample_id}.out.RDS" into filt_step3Trimming // needed for join() later
    file "forwardP.fa" into forprimers
    file "reverseP.fa" into revprimers
    file "forwardP_rc.fa" into rcfor
    file "reverseP_rc.fa" into rcrev

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")
    library(ShortRead); packageVersion("ShortRead")
    library(Biostrings); packageVersion("Biostrings")
    library(stringr); packageVersion("stringr")
    
    #Filter out reads with N's
    out1 <- filterAndTrim(fwd = "${reads[0]}",
                        filt = paste0("${sample_id}", ".R1.noN.fastq.gz"),
                        rev = "${reads[1]}",
                        filt.rev = paste0("${sample_id}", ".R2.noN.fastq.gz"),
                        maxN = 0,
                        matchIDs = as.logical(${params.matchIDs}),
                        multithread = ${task.cpus})
                        
    # Write out fasta of primers - Handles multiple primers
    Fprimer_name <- unlist(stringr::str_split("${params.fwdprimer_name}", ";"))
    Rprimer_name <- unlist(stringr::str_split("${params.revprimer_name}", ";"))
    
    Fprimers <- unlist(stringr::str_split("${params.fwdprimer}", ";"))
    names(Fprimers) <- Fprimer_name
    Rprimers <- unlist(stringr::str_split("${params.revprimer}", ";"))
    names(Rprimers) <- Rprimer_name
    
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(Fprimers), "forwardP.fa")
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(Rprimers), "reverseP.fa")
    
    # Write out fasta of reverse complement primers
    # Used for checking for read-through into the other end of molecule for variable length markers
    fwd_rc <- sapply(Fprimers, dada2:::rc)
    names(fwd_rc) <- Fprimer_name

    rev_rc <- sapply(Rprimers, dada2:::rc)
    names(rev_rc) <- Rprimer_name
    
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(fwd_rc), "forwardP_rc.fa")
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(rev_rc), "reverseP_rc.fa")
        
    saveRDS(out1, "${sample_id}.out.RDS")
    """
}

/*
 *
 * Step 3: Remove primers with cutadapt
 *
 */
 
// TODO: Add lengthvar demulti because cores isnt supported when demultiplexing-cores ${task.cpus} \\
if (params.lengthvar == false) {
    process cutadapt {
        tag { "filt_step2_${sample_id}" }

        input:
        set sample_id, reads from filt_step2
        file("forwardP.fa") from forprimers
        file("reverseP.fa") from revprimers
        
        output:
        set val(sample_id), "${sample_id}*.R[12].cutadapt.fastq.gz" optional true into filt_step3
        file "*.cutadapt.out" into cutadaptToMultiQC

        when:
        params.precheck == false

        script:
        """
        #!/bin/bash
        
        if [ ! -z "${params.demux}" ];
        then
        echo "More than one primer detected, demultiplexing (single core)";
        cutadapt \\
            -g file:forwardP.fa \\
            -G file:reverseP.fa \\
            -n 2 \\
            --no-indels \\
            -o "${sample_id}.{name}.R1.cutadapt.fastq.gz" \\
            -p "${sample_id}.{name}.R2.cutadapt.fastq.gz" \\
            "${reads[0]}" "${reads[1]}" > "${sample_id}.cutadapt.out"           

        else
        echo "Single primer detected (multi-core)";
        cutadapt \\
            -g "${params.fwdprimer}" \\
            -G "${params.revprimer}" \\
            --cores ${task.cpus} \\
            -n 2 \\
            -o "${sample_id}.R1.cutadapt.fastq.gz" \\
            -p "${sample_id}.R2.cutadapt.fastq.gz" \\
            "${reads[0]}" "${reads[1]}" > "${sample_id}.cutadapt.out"
        fi;
        """
    }
}
/* Length variable amplicon filtering - Trim both sides*/
else if (params.lengthvar == true) {
    //TODO: calculate the -e parameter in order to allow 1 mismatch (1/max primer length)
    process cutadapt_var {
        tag { "varfilt_step2_${sample_id}" }

        input:
        set sample_id, reads from filt_step2
        file("forwardP.fa") from forprimers
        file("reverseP.fa") from revprimers
        file("forwardP_rc.fa") from rcfor
        file("reverseP_rc.fa") from rcrev
        
        output:
        set val(sample_id), "${sample_id}*.R[12].cutadapt.fastq.gz" optional true into filt_step3
        file "*.cutadapt.out" into cutadaptToMultiQC

        when:
        params.precheck == false

        script:
        """
        #!/bin/bash
        
        if [ ! -z "${params.demux}" ];
        then
        echo "More than one primer detected, demultiplexing (single-core)";
        cutadapt \\
            -g file:forwardP.fa -a file:reverseP_rc.fa \\
            -G file:reverseP.fa -a file:forwardP_rc.fa\\
            --cores ${task.cpus} \\
            -n 2 \\
            --no-indels \\
            -o "${sample_id}.{name}.R1.cutadapt.fastq.gz" \\
            -p "${sample_id}.{name}.R2.cutadapt.fastq.gz" \\
            "${reads[0]}" "${reads[1]}" > "${sample_id}.cutadapt.out"       

        else
        echo "Single primer detected (multi-core)";
        fwd_rc=\$(cat forwardP_rc.fa | tail -1)
        rev_rc=\$(cat reverseP_rc.fa | tail -1)
        
        cutadapt \\
            -g "${params.fwdprimer}" -a ${rev_rc}\\
            -G "${params.revprimer}" -A ${fwd_rc}\\
            --cores ${task.cpus} \\
            -n 2 \\
            -o "${sample_id}.R1.cutadapt.fastq.gz" \\
            -p "${sample_id}.R2.cutadapt.fastq.gz" \\
            "${reads[0]}" "${reads[1]}" > "${sample_id}.cutadapt.out"
        fi;
        """
    }    
} else {
    // We need to shut this down!
    cutadaptToMultiQC = Channel.empty()
    filteredReads = Channel.empty()
    filteredReadsforQC = Channel.empty()
}


/*
 *
 * Step 4: Filter reads
 *
 */
 
process FilterAndTrim {
    tag { "filt_step3_${sample_id}" }

    input:
    set sample_id, file(reads), file(trimming) from filt_step3.join(filt_step3Trimming)

    output:
    set val(sample_id), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" optional true into filteredReadsforQC, filteredReads
    file "*.R1.filtered.fastq.gz" optional true into forReads
    file "*.R2.filtered.fastq.gz" optional true into revReads
    file "*.trimmed.txt" into trimTracking

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")
    library(ShortRead); packageVersion("ShortRead")
    library(Biostrings); packageVersion("Biostrings")
    library(stringr); packageVersion("stringr")
    
    fastqFs <- sort(list.files(pattern="*.R1.cutadapt.fastq.gz"))
    fastqFs <- fastqFs[!stringr::str_detect(fastqFs, "unknown")]
    fastqRs <- sort(list.files(pattern="*.R2.cutadapt.fastq.gz"))
    fastqRs <- fastqRs[!stringr::str_detect(fastqRs, "unknown")]

    out1 <- readRDS("${trimming}")
    out2 <- filterAndTrim(fwd = fastqFs,
                        filt = stringr::str_replace(fastqFs, ".cutadapt.fastq.gz", ".filtered.fastq.gz"),
                        rev = fastqRs,
                        filt.rev = stringr::str_replace(fastqRs, ".cutadapt.fastq.gz", ".filtered.fastq.gz"),
                        maxEE = c(${params.maxEEFor},${params.maxEERev}),
                        truncLen = c(${params.truncFor},${params.truncRev}),
                        truncQ = ${params.truncQ},
                        maxN = ${params.maxN},
                        rm.phix = as.logical(${params.rmPhiX}),
                        maxLen = ${params.maxLen},
                        minLen = ${params.minLen},
                        compress = TRUE,
                        verbose = TRUE,
                        matchIDs = as.logical(${params.matchIDs}),
                        multithread = ${task.cpus})
    #Change input read counts to actual raw read counts
    out3 <- cbind(out1, out2)
    colnames(out3) <- c('input', 'filterN', 'cutadapt', 'filtered')
    write.csv(out3, paste0("${sample_id}", ".trimmed.txt"))
    """
}


/*
 *
 * Step 5: Post-filter Quality control
 *
 */

process runFastQC_postfilterandtrim {
    tag { "rFQC_post_FT.${sample_id}" }
    publishDir "${params.outdir}/qc/FastQC-postfilter", mode: "copy", overwrite: true

    input:
    set val(sample_id), file(filtFor), file(filtRev) from filteredReadsforQC

    output:
    file '*_fastqc.{zip,html}' into fastqc_files_post

    when:
    params.precheck == false

    """
    fastqc --nogroup -q ${filtFor} ${filtRev}
    """
}

process runMultiQC_postfilterandtrim {
    tag { "runMultiQC_postfilterandtrim" }
    publishDir "${params.outdir}/qc/MultiQC-postfilter", mode: 'copy', overwrite: true

    input:
    file('./raw-seq/*') from fastqc_files2_ch.collect()
    file('./trimmed-seq/*') from fastqc_files_post.collect()
    file('./cutadapt/*') from cutadaptToMultiQC.collect()

    output:
    file "*_report.html" into multiqc_report_post
    file "*_data"

    when:
    params.precheck == false

    script:
    interactivePlots = params.interactiveMultiQC == true ? "-ip" : ""
    """
    multiqc ${interactivePlots} .
    """
}

// TODO: Make this work with demultiplexing
process mergeTrimmedTable {
    tag { "mergeTrimmedTable" }
    publishDir "${params.outdir}/csv", mode: "copy", overwrite: true

    input:
    file trimData from trimTracking.collect()

    output:
    file "all.trimmed.csv" into trimmedReadTracking

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    trimmedFiles <- list.files(path = '.', pattern = '*.trimmed.txt')
    sample.names <- sub('.trimmed.txt', '', trimmedFiles)
    trimmed <- do.call("rbind", lapply(trimmedFiles, function (x) as.data.frame(read.csv(x))))
    colnames(trimmed)[1] <- "Sequence"
    trimmed\$SampleID <- sample.names
    write.csv(trimmed, "all.trimmed.csv", row.names = FALSE)
    """
}

/*
 *
 * Step 6: Learn error rates (run on all samples)
 *
 */

process LearnErrors {
    tag { "LearnErrors" }
    publishDir "${params.outdir}/qc", mode: "copy", overwrite: true

    input:
    file fReads from forReads.collect()
    file rReads from revReads.collect()

    output:
    file "errorsF.RDS" into errorsFor
    file "errorsR.RDS" into errorsRev
    file "*.pdf"
    
    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")    
    setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})

    # File parsing
    filtFs <- list.files('.', pattern="R1.filtered.fastq.gz", full.names = TRUE)
    filtRs <- list.files('.', pattern="R2.filtered.fastq.gz", full.names = TRUE)
    
    sample.namesF <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
    sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
    set.seed(100)

    # Learn error rates
    errF <- learnErrors(filtFs,         
        nbases=${params.errbases},
        errorEstimationFunction=${params.errfun},
        multithread=${task.cpus},
        randomize=FALSE
        )
    errR <- learnErrors(filtRs,
        nbases=${params.errbases},
        errorEstimationFunction=${params.errfun},
        multithread=${task.cpus},
        randomize=FALSE
        )

    # optional NovaSeq binning error correction
    if (as.logical('${params.qualityBinning}') == TRUE ) {
        print("Running binning correction")
        errs <- t(apply(getErrors(errF), 1, function(x) { x[x < x[40]] = x[40]; return(x)} ))
        errF\$err_out <- errs
        errs <- t(apply(getErrors(errR), 1, function(x) { x[x < x[40]] = x[40]; return(x)} ))
        errR\$err_out <- errs
    }

    pdf("err.pdf")
    plotErrors(errF, nominalQ=TRUE)
    plotErrors(errR, nominalQ=TRUE)
    dev.off()

    saveRDS(errF, "errorsF.RDS")
    saveRDS(errR, "errorsR.RDS")
    """
}

/*
 *
 * Step 7: Dereplication, ASV Inference, Merge Pairs
 *
 */

if (params.pool == "T" || params.pool == 'pseudo') {

    process PoolSamplesInferDerepAndMerge {
        tag { "PoolSamplesInferDerepAndMerge" }
        publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

        // TODO: filteredReads channel has ID and two files, should fix this
        // with a closure, something like  { it[1:2] }, or correct the channel
        // as the ID can't be used anyway

        input:
        file filts from filteredReads.collect( )
        file errFor from errorsFor
        file errRev from errorsRev

        output:
        file "seqtab.RDS" into seqTable,rawSeqTableToRename
        file "all.mergers.RDS" into mergerTracking
        file "all.dadaFs.RDS" into dadaForReadTracking
        file "all.dadaRs.RDS" into dadaRevReadTracking
        file "seqtab.*"

        when:
        params.precheck == false

        script:
        """
        #!/usr/bin/env Rscript
        library(dada2); packageVersion("dada2")
        setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})
        filtFs <- list.files('.', pattern="R1.filtered.fastq.gz", full.names = TRUE)
        filtRs <- list.files('.', pattern="R2.filtered.fastq.gz", full.names = TRUE)

        errF <- readRDS("${errFor}")
        errR <- readRDS("${errRev}")
        cat("Processing all samples\n")

        #Variable selection from CLI input flag --pool
        pool <- "${params.pool}"
        if(pool == "T" || pool == "TRUE"){
          pool <- as.logical(pool)
        }
        dadaFss <- dada(filtFs, err=errF, multithread=${task.cpus}, pool=pool)
        dadaRss <- dada(filtRs, err=errR, multithread=${task.cpus}, pool=pool)

        mergers <- mergePairs(dadaFss, filtFs, dadaRss, filtRs,
            returnRejects = TRUE,
            minOverlap = ${params.minOverlap},
            maxMismatch = ${params.maxMismatch},
            trimOverhang = as.logical("${params.trimOverhang}"),
            justConcatenate = as.logical("${params.justConcatenate}")
            )

        # TODO: make this a single item list with ID as the name, this is lost
        # further on
        saveRDS(mergers, "all.mergers.RDS")

        saveRDS(dadaFss, "all.dadaFs.RDS")
        saveRDS(dadaRss, "all.dadaRs.RDS")

        # go ahead and make seqtable
        seqtab <- makeSequenceTable(mergers)
        saveRDS(seqtab, "seqtab.RDS")
        
        # Track reads
        getN <- function(x) sum(getUniques(x))
        dada_out <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN)) %>%
            magrittr::set_colnames(c("dadaFs", "dadaRs", "merged")) %>%
            as.data.frame() %>%
            rownames_to_column("sample_id") %>%
            mutate(sample_id = str_replace(basename(sample_id), pattern="_S.*$", replacement=""))
        
        # TODO: change this to flow cell id
        write.csv(dada_out, "dada_out.csv") 
        """
        }
} else {
    // pool = F, process per sample
    process PerSampleInferDerepAndMerge {
        tag { "PerSampleInferDerepAndMerge" }
        publishDir "${params.outdir}/qc", mode: "copy", overwrite: true

        input:
        set val(sample_id), file(filtFor), file(filtRev) from filteredReads
        file errFor from errorsFor
        file errRev from errorsRev

        output:
        file "seqtab.RDS" into seqTable
        file "all.mergers.RDS" into mergerTracking
        file "all.dadaFs.RDS" into dadaForReadTracking
        file "all.dadaRs.RDS" into dadaRevReadTracking
        file "seqtab.*"

        when:
        params.precheck == false

        script:
        """
        #!/usr/bin/env Rscript
        library(dada2); packageVersion("dada2")        
        setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})

        errF <- readRDS("${errFor}")
        errR <- readRDS("${errRev}")
        cat("Processing:", "${sample_id}", "\\n")
        
        filtFs <- "${filtFor}"
        filtRs <- "${filtRev}"

        dadaFs <- dada(filtFs, err=errF, multithread=${task.cpus}, pool=as.logical("${params.pool}"))
        dadaRs <- dada(filtRs, err=errR, multithread=${task.cpus}, pool=as.logical("${params.pool}"))

        merger <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
            returnRejects = TRUE,
            minOverlap = ${params.minOverlap},
            maxMismatch = ${params.maxMismatch},
            trimOverhang = as.logical("${params.trimOverhang}"),
            justConcatenate=as.logical("${params.justConcatenate}")
            )

        saveRDS(merger, paste("${sample_id}", "merged", "RDS", sep="."))

        saveRDS(dadaFss, "all.dadaFs.RDS")
        saveRDS(dadaRss, "all.dadaRs.RDS")
        """
    }

    process mergeDadaRDS {
        tag { "mergeDadaRDS" }
        publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

        input:
        file dadaFss from dadaFor.collect()
        file dadaRss from dadaRev.collect()

        output:
        file "all.dadaFs.RDS" into dadaForReadTracking
        file "all.dadaRs.RDS" into dadaRevReadTracking

        when:
        params.precheck == false

        script:
        '''
        #!/usr/bin/env Rscript
        library(dada2)
        packageVersion("dada2")

        dadaFs <- lapply(list.files(path = '.', pattern = '.dadaFs.RDS$'), function (x) readRDS(x))
        names(dadaFs) <- sub('.dadaFs.RDS', '', list.files('.', pattern = '.dadaFs.RDS'))
        dadaRs <- lapply(list.files(path = '.', pattern = '.dadaRs.RDS$'), function (x) readRDS(x))
        names(dadaRs) <- sub('.dadaRs.RDS', '', list.files('.', pattern = '.dadaRs.RDS'))
        saveRDS(dadaFs, "all.dadaFs.RDS")
        saveRDS(dadaRs, "all.dadaRs.RDS")
        '''
    }

    process SequenceTable {
        tag { "SequenceTable" }
        publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

        input:
        file mr from mergedReads.collect()

        output:
        file "seqtab.RDS" into seqTable,rawSeqTableToRename
        file "all.mergers.RDS" into mergerTracking

        when:
        params.precheck == false

        script:
        '''
        #!/usr/bin/env Rscript
        library(dada2); packageVersion("dada2")        
        
        mergerFiles <- list.files(path = '.', pattern = '.*.RDS$')
        pairIds <- sub('.merged.RDS', '', mergerFiles)
        mergers <- lapply(mergerFiles, function (x) readRDS(x))
        names(mergers) <- pairIds
        seqtab <- makeSequenceTable(mergers)
        seqtab <- seqtab[,nchar(colnames(seqtab)) >= ${params.minLen}]

        saveRDS(seqtab, "seqtab.RDS")
        saveRDS(mergers, "all.mergers.RDS")
        '''
    }
}

// TODO: Add seqtab merge here for multiple runs


/*
 *
 * Step 8: ASV filtering
 *
 */

// TODO: add length filter
// If coding, then do codon checks
// if PHMM provided, align to phmm, output this instead of later alignment below
if (params.coding) {
    process coding_asv_filter {
        tag { "coding_asv_filter" }
        publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

        input:
        file st from seqTable

        output:
        file "seqtab_final.RDS" into seqTableToTax,seqTableToRename

        when:
        params.precheck == false

        script:
        chimOpts = params.removeBimeraDenovoOptions != false ? ", ${params.removeBimeraDenovoOptions}" : ''
        """
        #!/usr/bin/env Rscript
        library(dada2); packageVersion("dada2")        
        library(tidyverse); packageVersion("tidyverse")
        library(biostrings); packageVersion("Biostrings")
        source("functions.R") # Needed for codon_filter & PHMM - taken from taxreturn
        
        st.all <- readRDS("${st}")

        # Remove chimeras
        seqtab_nochim <- removeBimeraDenovo(
            st.all, 
            method="consensus", 
            multithread=${task.cpus}, 
            verbose=TRUE ${chimOpts} 
            )

        #cut to expected size
        seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) %in% as.numeric(${params.min_asv_len}):as.numeric(${params.max_asv_len})]
        
        # TODO: ADD PHMM HERE       
        
        #Filter sequences containing stop codons
        seqs <- Biostrings::DNAStringSet(getSequences(seqtab_cut))
        codon_filt <- codon_filter(seqs, genetic_code = "${params.genetic_code}") # Internal function
        seqtab_final <- seqtab_cut[,colnames(seqtab_cut) %in% codon_filt]                
        
        saveRDS(seqtab_final, "seqtab_final.RDS")
        
        # summarise cleanup
        cleanup <- st.all %>%
          as.data.frame() %>%
          tidyr::pivot_longer( everything(),
            names_to = "OTU",
            values_to = "Abundance") %>%
          dplyr::group_by(OTU) %>%
          dplyr::summarise(Abundance = sum(Abundance)) %>%
          dplyr::mutate(length  = nchar(OTU)) %>%
          dplyr::mutate(type = case_when(
            !OTU %in% getSequences(seqtab_nochim) ~ "Chimera",
            !OTU %in% getSequences(seqtab_cut) ~ "Incorrect size",
            !OTU %in% getSequences(seqtab_final) ~ "Stop codons",
            TRUE ~ "Real"
          )) 
        write_csv(cleanup, "ASV_cleanup_summary.csv")
        
        # Output length distribution plots
        gg.abundance <- ggplot(cleanup, aes(x=length, y=Abundance, fill=type))+
            geom_bar(stat="identity") + 
            labs(title = "Abundance of sequences")

        gg.unique <- ggplot(cleanup, aes(x=length, fill=type))+
            geom_histogram() + 
            labs(title = "Number of unique sequences")

        pdf(paste0("seqtab_length_dist.pdf"), width = 11, height = 8 , paper="a4r")
          plot(gg.abundance / gg.unique)
        try(dev.off(), silent=TRUE)
        """
    }
} else {
    process noncoding_asv_filter {
        tag { "noncoding_asv_filter" }
        publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

        input:
        file st from seqTable

        output:
        file "seqtab_final.RDS" into seqTableToTax,seqTableToRename

        when:
        params.precheck == false

        script:
        chimOpts = params.removeBimeraDenovoOptions != false ? ", ${params.removeBimeraDenovoOptions}" : ''
        """
        #!/usr/bin/env Rscript
        library(dada2); packageVersion("dada2")        
        library(tidyverse); packageVersion("tidyverse")
        library(biostrings); packageVersion("Biostrings")
        
        st.all <- readRDS("${st}")

        # Remove chimeras
        seqtab_nochim <- removeBimeraDenovo(
            st.all, 
            method="consensus", 
            multithread=${task.cpus}, 
            verbose=TRUE ${chimOpts} 
            )

        #cut to expected size
        seqtab_final <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) %in% as.numeric(${params.min_asv_len}):as.numeric(${params.max_asv_len})]
        
        saveRDS(seqtab_final, "seqtab_final.RDS")        
        
        # summarise cleanup
        cleanup <- st.all %>%
          as.data.frame() %>%
          tidyr::pivot_longer( everything(),
            names_to = "OTU",
            values_to = "Abundance") %>%
          dplyr::group_by(OTU) %>%
          dplyr::summarise(Abundance = sum(Abundance)) %>%
          dplyr::mutate(length  = nchar(OTU)) %>%
          dplyr::mutate(type = case_when(
            !OTU %in% getSequences(seqtab_nochim) ~ "Chimera",
            !OTU %in% getSequences(seqtab_final) ~ "Incorrect size",
            TRUE ~ "Real"
          )) 
        readr::write_csv(cleanup, "ASV_cleanup_summary.csv")
        
        # Output length distribution plots
        gg.abundance <- ggplot(cleanup, aes(x=length, y=Abundance, fill=type))+
            geom_bar(stat="identity") + 
            labs(title = "Abundance of sequences")

        gg.unique <- ggplot(cleanup, aes(x=length, fill=type))+
            geom_histogram() + 
            labs(title = "Number of unique sequences")

        pdf(paste0("seqtab_length_dist.pdf"), width = 11, height = 8 , paper="a4r")
          plot(gg.abundance / gg.unique)
        try(dev.off(), silent=TRUE)
        """
    }
}

/*
 *
 * Step 9: Taxonomic assignment
 *
 */

// TODO: Add IDTAXA + Blast
if (params.reference) {
    refFile = file(params.reference)
    if (params.taxassignment == 'rdp') {
        // TODO: we could combine these into the same script
        

        if (params.species) {

            speciesFile = file(params.species)

            process AssignTaxSpeciesRDP {
                tag { "AssignTaxSpeciesRDP" }
                publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

                input:
                file st from seqTableToTax
                file ref from refFile
                file sp from speciesFile

                output:
                file "tax_final.RDS" into taxFinal,taxTableToTable
                file "bootstrap_final.RDS" into bootstrapFinal

                when:
                params.precheck == false

                script:
                """
                #!/usr/bin/env Rscript
                library(dada2); packageVersion("dada2")                

                seqtab <- readRDS("${st}")

                # Assign taxonomy
                tax <- assignTaxonomy(seqtab, "${ref}",
                                        multithread=${task.cpus},
                                        tryRC = TRUE,
                                        outputBootstraps = TRUE,
                                        minBoot = ${params.minBoot},
                                        verbose = TRUE)
                boots <- tax\$boot

                tax <- addSpecies(tax\$tax, "${sp}",
                                 tryRC = TRUE,
                                 verbose = TRUE)

                rownames(tax) <- colnames(seqtab)

                # Write original data
                saveRDS(tax, "tax_final.RDS")
                saveRDS(boots, "bootstrap_final.RDS")
                """
            }

        } else {

            process AssignTaxonomyRDP {
                tag { "TaxonomyRDP" }
                publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

                input:
                file st from seqTableToTax
                file ref from refFile

                output:
                file "tax_final.RDS" into taxFinal,taxTableToTable
                file "bootstrap_final.RDS" into bootstrapFinal

                when:
                params.precheck == false

                script:
                taxLevels = params.taxLevels ? "c( ${params.taxLevels} )," : ''
                """
                #!/usr/bin/env Rscript
                library(dada2); packageVersion("dada2")
                
                seqtab <- readRDS("${st}")

                # Assign taxonomy
                tax <- assignTaxonomy(seqtab, "${ref}",
                                      multithread=${task.cpus},
                                      minBoot = ${params.minBoot},
                                      tryRC = TRUE,
                                      outputBootstraps = TRUE, ${taxLevels}
                                      verbose = TRUE 
                                      )

                # Write to disk
                saveRDS(tax\$tax, "tax_final.RDS")
                saveRDS(tax\$boot, "bootstrap_final.RDS")
                """
            }
        }
    } else if (params.taxassignment == 'idtaxa') {
        process TaxonomyIDTAXA {
            tag { "TaxonomyIDTAXA" }
            publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

            input:
            file st from seqTableToTax
            file ref from refFile // this needs to be a database from the IDTAXA site

            output:
            file "tax_final.RDS" into taxFinal,taxTableToTable
            file "bootstrap_final.RDS" into bootstrapFinal
            file "raw_idtaxa.RDS"

            when:
            params.precheck == false

            script:
            """
            #!/usr/bin/env Rscript
            library(dada2); packageVersion("dada2")
            library(DECIPHER); packageVersion("DECIPHER")
            library(stringr); packageVersion("stringr")

            seqtab <- readRDS("${st}")

            # Create a DNAStringSet from the ASVs
            dna <- DNAStringSet(getSequences(seqtab))

            # load database; this should be a RData file
            if (stringr::str_detect("${refFile}", ".RData")){
                load("${refFile}")
            } else if(stringr::str_detect("${refFile}", ".rds")){
                trainingSet <- readRDS("${refFile}")
            }

            ids <- IdTaxa(dna, trainingSet,
                strand="both",
                processors=${task.cpus},
                verbose=TRUE)
            # ranks of interest
            ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
            saveRDS(ids, 'raw_idtaxa.RDS')

            # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
            taxid <- t(sapply(ids, function(x) {
                    m <- match(ranks, x\$rank)
                    taxa <- x\$taxon[m]
                    taxa[startsWith(taxa, "unclassified_")] <- NA
                    taxa
            }))
            colnames(taxid) <- ranks
            rownames(taxid) <- getSequences(seqtab)

            boots <- t(sapply(ids, function(x) {
                    m <- match(ranks, x\$rank)
                    bs <- x\$confidence[m]
                    bs
            }))
            colnames(boots) <- ranks
            rownames(boots) <- getSequences(seqtab)

            # Write to disk
            saveRDS(taxid, "tax_final.RDS")
            saveRDS(boots, "bootstrap_final.RDS")
            """
        }

    } else if (params.taxassignment) {
        exit 1, "Unknown taxonomic assignment method set: ${params.taxassignment}"
    } else {
        exit 1, "No taxonomic assignment method set, but reference passed"
    }
} else {
    // set tax channels to 'false', do NOT assign taxonomy
    taxFinal = Channel.empty()
    taxTableToTable = Channel.empty()
    bootstrapFinal = Channel.empty()
}

// Note: this is currently a text dump.  We've found the primary issue with
// downstream analysis is getting the data in a form that can be useful as
// input, and there isn't much consistency with this as of yet.  So for now
// we're using the spaghetti approach (see what sticks).  Also, we  are running
// into issues with longer sequences (e.g. concatenated ones) used as IDs with
// tools like Fasttree (it doesn't seem to like that).

// Safest way may be to save the simpleID -> seqs as a mapping file, use that in
// any downstream steps (e.g. alignment/tree), then munge the seq names back
// from the mapping table

/*
 *
 * Step 10: Rename ASVs & Generate seuqence tables
 *
 * A number of downstream programs have issues with sequences as IDs, here we
 * (optionally) rename these
 *
 */

process RenameASVs {
    tag { "RenameASVs" }
    publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

    input:
    file st from seqTableToRename
    file rawst from rawSeqTableToRename

    output:
    file "seqtab_final.simple.RDS" into seqTableFinalToBiom,seqTableFinalToTax,seqTableFinalTree,seqTableFinalTracking,seqTableToTable,seqtabToPhyloseq,seqtabToTaxTable
    file "asvs.${params.idType}.nochim.fna" into seqsToAln, seqsToQIIME2
    file "readmap.RDS" into readsToRenameTaxIDs // needed for remapping tax IDs
    file "asvs.${params.idType}.raw.fna"

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")
    library(ShortRead); packageVersion("ShortRead")
    library(digest); packageVersion("digest")

    # read RDS w/ data
    st <- readRDS("${st}")
    st.raw <- readRDS("${rawst}")

    # get sequences
    seqs <- colnames(st)
    seqs.raw <- colnames(st.raw)

    # get IDs based on idType
    ids_study <- switch("${params.idType}", simple=paste("ASV", 1:ncol(st), sep = ""),
                                md5=sapply(colnames(st), digest, algo="md5"))
    ids_study.raw <- switch("${params.idType}", simple=paste("ASV", 1:ncol(st.raw), sep = ""),
                                md5=sapply(colnames(st.raw), digest, algo="md5"))
    
    # sub IDs
    colnames(st) <- ids_study
    colnames(st.raw) <- ids_study.raw

    # generate FASTA
    seqs.dna <- ShortRead(sread = DNAStringSet(seqs), id = BStringSet(ids_study))
    # Write out fasta file.
    writeFasta(seqs.dna, file = 'asvs.${params.idType}.nochim.fna')

    seqs.dna.raw <- ShortRead(sread = DNAStringSet(seqs.raw), id = BStringSet(ids_study.raw))
    writeFasta(seqs.dna.raw, file = 'asvs.${params.idType}.raw.fna')

    # Write modified data (note we only keep the no-chimera reads for the next stage)
    saveRDS(st, "seqtab_final.simple.RDS")
    saveRDS(data.frame(id = ids_study, seq = seqs), "readmap.RDS")
    """
}

process GenerateSeqTables {
    tag { "GenerateSeqTables" }
    publishDir "${params.outdir}/rds", mode: "link", overwrite: true

    input:
    file st from seqTableToTable

    output:
    file "seqtab_final.simple.qiime2.txt" into featuretableToQIIME2
    file "*.txt"

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")
    library(ShortRead); packageVersion("ShortRead")

    seqtab <- readRDS("${st}")

    if (as.logical('${params.sampleRegex}' != FALSE )) {
        rownames(seqtab) <- gsub('${params.sampleRegex}', "\\\\1", rownames(seqtab), perl = TRUE)
    }

    # Generate table output
    write.table(data.frame('SampleID' = row.names(seqtab), seqtab),
        file = 'seqtab_final.txt',
        row.names = FALSE,
        col.names=c('#SampleID', colnames(seqtab)), sep = "\t")

    ######################################################################
    # Convert to simple table + FASTA, from
    # https://github.com/LangilleLab/microbiome_helper/blob/master/convert_dada2_out.R#L69
    ######################################################################

    # Generate OTU table output (rows = samples, cols = ASV)
    write.table(data.frame('SampleID' = row.names(seqtab), seqtab),
        file = 'seqtab_final.simple.txt',
        row.names = FALSE,
        col.names=c('#SampleID', colnames(seqtab)),
        sep = "\t")

    # Generate OTU table for QIIME2 import (rows = ASVs, cols = samples)
    write.table(
        data.frame('Taxa' = colnames(seqtab), t(seqtab)),
        file = 'seqtab_final.simple.qiime2.txt',
        row.names = FALSE,
        quote=FALSE,
        sep = "\t")

    # Write modified data
    saveRDS(seqtab, "seqtab_final.simple.RDS")
    """
}

process GenerateTaxTables {
    tag { "GenerateTaxTables" }
    publishDir "${params.outdir}/rds", mode: "link", overwrite: true

    input:
    file tax from taxTableToTable
    file bt from bootstrapFinal
    file map from readsToRenameTaxIDs

    output:
    file "tax_final.simple.RDS" into taxtabToPhyloseq
    file "tax_final.simple.txt" into taxtableToQIIME2
    file "*.txt"

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")
    library(ShortRead); packageVersion("ShortRead")

    tax <- readRDS("${tax}")
    map <- readRDS("${map}")

    # Note that we use the old ASV ID for output here
    write.table(data.frame('ASVID' = row.names(tax), tax),
        file = 'tax_final.txt',
        row.names = FALSE,
        col.names=c('#OTU ID', colnames(tax)), sep = "\t")

    # Tax table
    if(!identical(rownames(tax), as.character(map\$seq))){
        stop("sequences in taxa and sequence table are not ordered the same.")
    }

    tax[is.na(tax)] <- "Unclassified"
    rownames(tax) <- map\$id
    taxa_combined <- apply(tax, 1, function(x) paste(x, collapse=";"))
    taxa_out <- data.frame(names(taxa_combined), taxa_combined)
    colnames(taxa_out) <- c("#OTU ID", "taxonomy")

    write.table(data.frame('ASVID' = row.names(tax), tax),
        file = 'tax_final.simple.full.txt',
        row.names = FALSE,
        col.names=c('#OTU ID', colnames(tax)), sep = "\t")

    write.table(taxa_out,
        file = 'tax_final.simple.txt',
        row.names = FALSE,
        sep = "\t")

    if (file.exists('bootstrap_final.RDS')) {
        boots <- readRDS("${bt}")
        if(!identical(rownames(boots), as.character(map\$seq))){
            stop("sequences in bootstrap and sequence table are not ordered the same.")
        }
        rownames(boots) <- map\$id
        write.table(data.frame('ASVID' = row.names(boots), boots),
            file = 'tax_final.bootstraps.simple.full.txt',
            row.names = FALSE,
            col.names=c('#OTU ID', colnames(boots)), sep = "\t")
    }

    # Write modified data
    saveRDS(tax, "tax_final.simple.RDS")
    """
}

/*
 *
 * Step 11: Alignment & Phylogenetic tree
 *
 */

// NOTE: 'when' directive doesn't work if channels have the same name in
// two processes
// TODO: if PHMM provided, subset the previous alignment to just the retained sequences and output that

if (!params.precheck && params.runTree && params.lengthvar == false) {

    if (params.aligner == 'DECIPHER') {

        process AlignReadsDECIPHER {
            tag { "AlignReadsDECIPHER" }
            publishDir "${params.outdir}/fasta", mode: "copy", overwrite: true
            errorStrategy 'ignore'

            input:
            file seqs from seqsToAln

            output:
            file "aligned_seqs.fasta" optional true into alnFile,alnToQIIME2
            
            script:
            """
            #!/usr/bin/env Rscript
            library(dada2); packageVersion("dada2")
            library(DECIPHER); packageVersion("DECIPHER")

            seqs <- readDNAStringSet("${seqs}")
            alignment <- AlignSeqs(seqs,
                                   anchor=NA,
                                   processors = ${task.cpus})
            writeXStringSet(alignment, "aligned_seqs.fasta")
            """
        }
    } else {
        exit 1, "Unknown aligner option: ${params.aligner}"
    }

    /*
     *
     * Step 10b: Construct phylogenetic tree
     *
     */

    if (params.runTree == 'phangorn') {

        process GenerateTreePhangorn {
            tag { "GenerateTreePhangorn" }
            publishDir "${params.outdir}/trees", mode: "copy", overwrite: true

            input:
            file aln from alnFile

            output:
            file "phangorn.tree.RDS" into treeRDS
            file "tree.newick" into treeFile
            file "tree.GTR.newick" into treeGTRFile

            script:
            """
            #!/usr/bin/env Rscript
            library(phangorn); packageVersion("phangorn")

            phang.align <- read.phyDat("aligned_seqs.fasta",
                                        format = "fasta",
                                        type = "DNA")

            dm <- dist.ml(phang.align)
            treeNJ <- NJ(dm) # Note, tip order != sequence order
            fit = pml(treeNJ, data=phang.align)
            write.tree(fit\$tree, file = "tree.newick")

            ## negative edges length changed to 0!
            fitGTR <- update(fit, k=4, inv=0.2)
            fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                                  rearrangement = "stochastic", control = pml.control(trace = 0))
            saveRDS(fitGTR, "phangorn.tree.RDS")
            write.tree(fitGTR\$tree, file = "tree.GTR.newick")
            """
        }
    } else if (params.runTree == 'fasttree') {

        process GenerateTreeFasttree {
            tag { "GenerateTreeFasttree" }
            publishDir "${params.outdir}/trees", mode: "copy", overwrite: true

            input:
            file aln from alnFile

            output:
            file "fasttree.tree" into treeGTRFile, treeToQIIME2
            // need to deadend the other channels, they're hanging here

            script:
            """
            OMP_NUM_THREADS=${task.cpus} FastTree -nt \\
                -gtr -gamma -spr 4 -mlacc 2 -slownni \\
                -out fasttree.tree \\
                aligned_seqs.fasta
            """
        }

    } else {
        // dead-end channels generated above
    }

    process RootTree {
        tag { "RootTree" }
        publishDir "${params.outdir}/trees", mode: "link"

        input:
        file tree from treeGTRFile

        output:
        file "rooted.newick" into rootedTreeFile, rootedToQIIME2
        // need to deadend the other channels, they're hanging here

        script:
        """
        #!/usr/bin/env Rscript
        library(phangorn); packageVersion("phangorn")
        library(ape); packageVersion("ape")

        tree <- read.tree(file = "${tree}")

        midtree <- midpoint(tree)

        write.tree(midtree, file = "rooted.newick")
        """
    }
} else {
    // Note these are caught downstream
    alnToQIIME2 = false
    treeToQIIME2 = false
    rootedToQIIME2 = false
}

/*
 *
 * Step 12: Track reads
 *
 */

// Broken?: needs a left-join on the initial table

process ReadTracking {
    tag { "ReadTracking" }
    publishDir "${params.outdir}/qc", mode: "copy", overwrite: true

    input:
    file trimmedTable from trimmedReadTracking
    file sTable from seqTableFinalTracking
    file mergers from mergerTracking
    file dadaFss from dadaForReadTracking
    file dadaRss from dadaRevReadTracking

    output:
    file "all.readtracking.txt"

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")
    library(dplyr); packageVersion("dplyr")

    getN <- function(x) sum(getUniques(x))

    # the gsub here might be a bit brittle...
    dadaFs <- as.data.frame(sapply(readRDS("${dadaFss}"), getN))
    rownames(dadaFs) <- gsub('.R1.filtered.fastq.gz', '',rownames(dadaFs))
    colnames(dadaFs) <- c("denoisedF")
    dadaFs\$SampleID <- rownames(dadaFs)

    dadaRs <- as.data.frame(sapply(readRDS("${dadaRss}"), getN))
    rownames(dadaRs) <- gsub('.R2.filtered.fastq.gz', '',rownames(dadaRs))
    colnames(dadaRs) <- c("denoisedR")
    dadaRs\$SampleID <- rownames(dadaRs)

    all.mergers <- readRDS("${mergers}")
    mergers <- as.data.frame(sapply(all.mergers, function(x) sum(getUniques(x %>% filter(accept)))))
    rownames(mergers) <- gsub('.R1.filtered.fastq.gz', '',rownames(mergers))
    colnames(mergers) <- c("merged")
    mergers\$SampleID <- rownames(mergers)

    seqtab.nochim <- as.data.frame(rowSums(readRDS("${sTable}")))
    rownames(seqtab.nochim) <- gsub('.R1.filtered.fastq.gz', '',rownames(seqtab.nochim))
    colnames(seqtab.nochim) <- c("seqtab.nochim")
    seqtab.nochim\$SampleID <- rownames(seqtab.nochim)

    trimmed <- read.csv("${trimmedTable}")

    track <- Reduce(function(...) merge(..., by = "SampleID",  all.x=TRUE),  list(trimmed, dadaFs, dadaRs, mergers, seqtab.nochim))
    # dropped data in later steps gets converted to NA on the join
    # these are effectively 0
    track[is.na(track)] <- 0
    
    write.table(track, "all.readtracking.txt", sep = "\t", row.names = FALSE)
    """
}

// TODO: Add phyloseq object

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    def subject = "[${params.base}/piperline] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[${params.base}/piperline] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[${params.base}/piperline] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[${params.base}/piperline] Sent summary e-mail to $params.email (mail)"
        }
    }
}
