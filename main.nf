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
    nextflow run alexpiper/piperline --reads '*_R{1,2}_001.fastq.gz' --lengthvar false \
    --fwdprimer 'GGDACWGGWTGAACWGTWTAYCCHCC' --fwdprimer_name 'fwhF2' --revprimer 'GTRATWGCHCCDGCTARWACWGG' --revprimer_name 'fwhR2n' \
    --reference 'idtaxa_bftrimmed.rds' -profile basc 

    The typical command for running the pipeline with your own config (instead of command line flags) is as follows:
    nextflow run alexpiper/piperline -c dada2_user_input.config -profile uct_hex
    where:
    dada2_user_input.config is the configuration file (see example 'dada2_user_input.config')
    NB: -profile uct_hex still needs to be specified from the command line

    To override existing values from the command line, please type these parameters:

    Mandatory inputs:
      --reads                       Path to input data (must be surrounded with quotes)
      --read_format                 Format of input read filenames. default: "fcid_sampleid_ext_pcr_S_R_L.fastq.gz" (must be surrounded with quotes)
      --samplesheet                 SampleSheet.csv file used for the sequencing run (must be surrounded with quotes)
      --runparams                   RunParameters.xml file from the sequencing run (must be surrounded with quotes)
      --reference                   Path to taxonomic database to be used for annotation (e.g. gg_13_8_train_set_97.fa.gz) (must be surrounded with quotes)
      --fwdprimer                   Sequence of the forward primer (must be surrounded with quotes)
      --revprimer                   Sequence of the reverse primer (must be surrounded with quotes)
      --fwdprimer_name              name of the forward primer (must be surrounded with quotes)
      --revprimer_name              name of the reverse primer (must be surrounded with quotes)
      -profile                      Hardware config to use. Currently profile available for BASC 'basc' - create your own if necessary
                                    NB -profile should always be specified on the command line, not in the config file

    Optional inputs:
      --interop                     Path to InterOp directory from the sequencing run (must be surrounded with quotes)
      --calc_switchrate             Option to calculate index switching rate using undetermined reads (default true)

    Read filtering parameters:
      --trimFor                     integer. The number of nucleotides to remove from the start of read 1
      --trimRev                     integer. The number of nucleotides to remove from the start of read 2
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
                                    
    Merging arguments (optional):
      --minOverlap                  The minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 The maximum mismatches allowed in the overlap region; default=0.
      --trimOverhang                If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off. "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen
                                    when reads are longer than the amplicon and read into the other-direction primer region. Default="F" (false)
      --minMergedLen                Minimum length of fragment *after* merging
      --maxMergedLen                Maximum length of fragment *after* merging                              

    Taxonomic arguments (optional):
      --species_db                   Specify path to fasta file. See dada2 addSpecies() for more detail.  
      --blast_db                     Specify path to fasta file.
      
    Sample and taxon filtering parameters (optional):
      --min_sample_reads            integer. The minimum number of reads a sample must have to be retained
      --phylum                      semicolon delimited string indicating the phyla to retain
      --order                       semicolon delimited string indicating the orders to retain
      --family                      semicolon delimited string indicating the families to retain
      --genus                       semicolon delimited string indicating the genera to retain
      
    Outputs
      --check_ala                   Whether to check taxon detections against the Atlas of Living Australia.
      --check_afd                   Whether to check taxon detections against the Australian Faunal Directory.

    Other arguments (optional):
      --dadaOpt.XXX                 Set as e.g. --dadaOpt.HOMOPOLYMER_GAP_PENALTY=-1 Global defaults for the dada function, see ?setDadaOpt in R for available options and their defaults
      --subsample                   subsample a random number of sequencing reads from each fastq.
      --subsample_seed              Seed for random subsampling
      --pool                        Should sample pooling be used to aid identification of low-abundance ASVs? Options are
                                    pseudo pooling: "pseudo", true: "T", false: "F"
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run
                                    sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    Help:
      --help                        Will print out summary above when executing nextflow run alexpiper/piperline




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

if ( params.reference == false ) {
    exit 1, "Must set reference database using --reference"
}

if ( params.samplesheet == false ) {
    exit 1, "Must set samplesheet using --samplesheet"
}

if ( params.runparams == false ) {
    exit 1, "Must set runparams using --runparams"
}

if (params.fwdprimer == false){
    exit 1, "Must set forward primer using --fwdprimer"
}

if (params.revprimer == false){
    exit 1, "Must set reverse primer using --revprimer"
}

if (params.fwdprimer_name == false){
    exit 1, "Must set name of the reverse primer using--fwdprimer_name"
}

if (params.revprimer_name == false){
    exit 1, "Must set name of the reverse primer using --revprimer_name"
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
        
Channel
    .fromPath( params.samplesheet )
    .ifEmpty { error "Cannot find any sample sheet matching: ${params.samplesheet}" }
    .set { samplesheet_ch }
    
Channel
    .fromPath( params.runparams )
    .ifEmpty { error "Cannot find any runparameters file matching: ${params.runparams}" }
    .set { runparams_ch }

if(params.phmm != false){
    Channel
        .fromPath( params.phmm )
        //.ifEmpty { "" }
        .set { phmm_ch }
}
else {
    phmm_ch = Channel.empty()
}

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
summary['species_db']     = params.species_db
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
    summary['E-mail address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 *
 * Step 0: Copy read files and validate read files
 *
 */

// TODO: Add FCID to start if cant find it - Get from runparameters
// TODO: Parse tuple rather than files to samples_to_validate

process validate_reads {
    tag { "validate_reads.${fastq_id}" }

    input:
    tuple fastq_id, file(reads) from samples_ch

    output:
    tuple fastq_id, file("data/${fastq_id}_R1_001.fastq.gz"), file("data/${fastq_id}_R2_001.fastq.gz") into samples_to_qual
    tuple fastq_id, env(fcid), env(sampleid), env(ext_id), env(pcr_id), file("data/*R[12]_001.fastq.gz") into samples_to_filt
    file ('*fastq_id.txt') into samples_to_validate
    tuple fastq_id, file("data/${fastq_id}_R1_001.fastq.gz") into samples_to_index
    
    script:
    """
    #!/bin/bash
    mkdir data
    
    # Subsample if a value has been provided to params.subsample
    if [ -z "${params.subsample}" ]
    then
        cp ${reads[0]} data/${fastq_id}_R1_001.fastq.gz
        cp ${reads[1]} data/${fastq_id}_R2_001.fastq.gz
    else
        seqtk sample -s"${params.subsample_seed}" ${reads[0]} "${params.subsample}" | pigz -p ${task.cpus} > data/${fastq_id}_R1_001.fastq.gz
        seqtk sample -s"${params.subsample_seed}" ${reads[1]} "${params.subsample}" | pigz -p ${task.cpus} > data/${fastq_id}_R2_001.fastq.gz
    fi  
    
    # Get expected positions of elements from read_format
    fcid_pos="\$(echo "${params.read_format}" | tr '_' '\n' | grep -n 'fcid' | cut -d : -f 1 )"
    sampleid_pos="\$(echo "${params.read_format}" | tr '_' '\n' | grep -Fxn 'sampleid' | cut -d : -f 1 )"
    ext_pos="\$(echo "${params.read_format}" | tr '_' '\n' | grep -Fxn 'ext' | cut -d : -f 1 )"
    pcr_pos="\$(echo "${params.read_format}" | tr '_' '\n' | grep -Fxn 'pcr' | cut -d : -f 1 )"
    
    # Check filepaths match input_format  
    if [[ "${fastq_id}" =~ .*"Undetermined".* ]]; then
        echo "Undetermined reads file."
        # Extract elements                      
        fcid="\$(echo "${fastq_id}" | cut -d'_' -f\${fcid_pos})"
        sampleid="\$(echo "${fastq_id}" | cut -d'_' -f\${sampleid_pos})"
        ext_id="ext1"
        pcr_id="pcr1"
    else
        # Check number of underscores match the read format
        read_format_start="\$(echo "${params.read_format}" | sed 's/_R.*\$//g')"
        fmt_split="\$(awk -F"_" '{print NF-1}' <<< "\${read_format_start}")"
        name_split="\$(awk -F"_" '{print NF-1}' <<< "${fastq_id}")"
        
        # if different - exit
        if [ "\${fmt_split}" -ne "\${name_split}" ]; then
            echo 'Number of underscores in filename "${fastq_id}" do not match expected read format: "${params.read_format}"'
            exit 1
        fi
        
        # Extract elements
        fcid="\$(echo "${fastq_id}" | cut -d'_' -f\${fcid_pos})"
        sampleid="\$(echo "${fastq_id}" | cut -d'_' -f\${sampleid_pos})"
        ext_id="\$(echo "${fastq_id}" | cut -d'_' -f\${ext_pos})"
        pcr_id="\$(echo "${fastq_id}" | cut -d'_' -f\${pcr_pos})"         
    fi
    
    # write out info for samplesheet creation
    echo "${fastq_id} \${fcid} \${sampleid} \${ext_id} \${pcr_id}" > "${fastq_id}_fastq_id.txt"
    """
}


/*
 *
 * Step 1: Pre-filter Quality control
 *
 */

process fastqc {
    tag { "rFQC.${fastq_id}" }
    publishDir "${params.outdir}/qc/FASTQC-prefilter", mode: "copy", overwrite: true

    input:
    tuple fastq_id, file(For), file(Rev) from samples_to_qual

    output:
    file '*_fastqc.{zip,html}' into fastqc_files_ch, fastqc_files2_ch

    """
    fastqc --nogroup -q ${For} ${Rev}
    """
}


// TODO: combine MultiQC reports and split by directory (no need for two)
process multiqc {
    tag { "multiqc" }
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

// TODO add a filter on _indexes.txt before collectFile instead of outputting empty undetermined_counts.txt files
if (params.calc_switchrate == true) {
    // Summarise indexes used for each sample
    process summarise_index {
        tag { "summarise_index_${fastq_id}" }

        input:
        set fastq_id, file(For) from samples_to_index

        output:
        file "*_indexes.txt" into index_count
        file "*_undetermined.txt" into undetermined_counts

        script:
        """
        #!/bin/bash        
        if [[ "${fastq_id}" =~ .*"Undetermined".* ]]; then
            echo "Undetermined reads file."
            zcat "${For}" | grep '^@M' | rev | cut -d':' -f 1 | rev | sort | uniq -c | sort -nr  | sed 's/+/ /' | sed 's/^ *//g' > ${fastq_id}_undetermined.txt
            touch ${fastq_id}_indexes.txt
        else
            zcat "${For}" | grep '^@M' | rev | cut -d':' -f 1 | rev | sort | uniq -c | sort -nr  | sed 's/+/ /' | sed 's/^ *//g' > ${fastq_id}_indexes.txt
            touch ${fastq_id}_undetermined.txt
        fi   
        """
    }
    // Calculate switch rate
    process index_calc {
        tag { "index_calc" }
        publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true
    
        input:
        file(indexes) from index_count.collectFile(name: 'determined_counts.txt', newLine: true)
        file(undetermined) from undetermined_counts.collectFile(name: 'undetermined_counts.txt', newLine: true)
        
        output:
        file "index_switch_calc.txt"
    
        script:
        """
        #!/bin/bash
        
        # Remove empty lines of collated files
        sed -i '/^\$/d' determined_counts.txt
        sed -i '/^\$/d' undetermined_counts.txt
        
        # Get all potential switched combinations of used indexes
        index1="\$(cat determined_counts.txt | cut -d' ' -f 2)"
        index2="\$(cat determined_counts.txt | cut -d' ' -f 3)"
        
        # get all possible combinations of determined indexes
        touch all_combinations.txt
        for i in \${index1} ; do
          for j in \${index2} ; do
            if [ "\${i}" \\< "\${j}" ]
            then
             echo \${i} \${j} >> all_combinations.txt
            fi
          done
        done
            
        # Count number of undetermined reads
        cat undetermined_counts.txt | cut -d' ' -f 2,3 > undetermined_index.txt
    
        # Count number of correctly demultiplexed reads
        correct_counts="\$(cat determined_counts.txt | cut -d' ' -f 1 | awk '{ SUM += \$1} END { print SUM }')"
    
        # Count number of switched reads
        comm -12 <(sort all_combinations.txt) <(sort undetermined_index.txt) > switched_indexes.txt
        switched_counts="\$(grep -f "switched_indexes.txt" "undetermined_counts.txt" | cut -d' ' -f 1 | awk '{ SUM += \$1} END { print SUM }')"
    
        # Count number of other reads (these can be sequencing errors, PhiX and other junk)
        other_counts="\$(grep -v -f "switched_indexes.txt" "undetermined_counts.txt" | cut -d' ' -f 1 | awk '{ SUM += \$1} END { print SUM }')"
    
        # Calculate switch rate (in percentage)
        calc(){ awk "BEGIN { print "\$*" }"; }
        switch_rate="\$(calc "\${switched_counts}"/"\${correct_counts}")"
        switch_rate_perc="\$(calc "\${switched_counts}"/"\${correct_counts}"*100)"
    
        # Print results to file
        touch index_switch_calc.txt
        echo Correctly demultiplexed reads: "\${correct_counts}" >> index_switch_calc.txt
        echo Switched reads: "\${switched_counts}" >> index_switch_calc.txt
        echo Other undetermined reads: "\${other_counts}" >> index_switch_calc.txt
        echo Index switching rate: "\${switch_rate}" \\("\${switch_rate_perc}"%\\) >> index_switch_calc.txt
        """
    }
}
else if (params.calc_switchrate == false) {
    // Set channels to empty
    index_count = Channel.empty()
    undetermined_counts = Channel.empty()
}


/*
 *
 * Step 2: Filter N bases
 *
 */
process nfilter {
    tag { "nfilter_${fastq_id}" }

    input:
    tuple fastq_id, fcid, sampleid, ext_id, pcr_id, reads from samples_to_filt.filter { !it[0].contains('Undetermined') }
    
    output:
    set val(fastq_id), val(fcid), val(sampleid), val(ext_id), val(pcr_id), "${fastq_id}.R[12].noN.fastq.gz" optional true into filt_step2
    set val(fastq_id), "${fastq_id}.out.RDS" into filt_step3Trimming // needed for join() later
    file "forwardP.fa" into forprimers
    file "reverseP.fa" into revprimers
    file "forwardP_rc.fa" into rcfor
    file "reverseP_rc.fa" into rcrev

    script:
    """
    #!/usr/bin/env Rscript
    require(dada2); packageVersion("dada2")
    require(ShortRead); packageVersion("ShortRead")
    require(Biostrings); packageVersion("Biostrings")
    require(stringr); packageVersion("stringr")
    
    #Filter out reads with N's
    out1 <- filterAndTrim(fwd = "${reads[0]}",
                        filt = paste0("${fastq_id}", ".R1.noN.fastq.gz"),
                        rev = "${reads[1]}",
                        filt.rev = paste0("${fastq_id}", ".R2.noN.fastq.gz"),
                        maxN = 0,
                        matchIDs = as.logical(${params.matchIDs}),
                        multithread = ${task.cpus})
                        
    # Write out fasta of primers - Handles multiple primers
    Fprimer_name <- unlist(stringr::str_split("${params.fwdprimer_name}", ";"))
    Rprimer_name <- unlist(stringr::str_split("${params.revprimer_name}", ";"))
    
    Fprimers <- unlist(stringr::str_split("${params.fwdprimer}", ";"))
    names(Fprimers) <- paste0(Fprimer_name,"-", Rprimer_name)
    Rprimers <- unlist(stringr::str_split("${params.revprimer}", ";"))
    names(Rprimers) <- paste0(Fprimer_name,"-", Rprimer_name)
    
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(Fprimers), "forwardP.fa")
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(Rprimers), "reverseP.fa")
    
    # Write out fasta of reverse complement primers
    # Used for checking for read-through into the other end of molecule for variable length markers
    fwd_rc <- sapply(Fprimers, dada2:::rc)
    names(fwd_rc) <- paste0(Fprimer_name,"-", Rprimer_name)

    rev_rc <- sapply(Rprimers, dada2:::rc)
    names(rev_rc) <- paste0(Fprimer_name,"-", Rprimer_name)
    
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(fwd_rc), "forwardP_rc.fa")
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(rev_rc), "reverseP_rc.fa")
        
    saveRDS(out1, "${fastq_id}.out.RDS")
    """
}

/*
 *
 * Step 3: Remove primers with cutadapt
 *
 */
 
// TODO:Drop single primer and multi primer workflow now that both can work with multicore
// TODO: Add a ^ to the start of the primers to make it restrict left?

process cutadapt {
    tag { "filt_step2_${fastq_id}" }

    input:
    set fastq_id, fcid, sampleid, ext_id, pcr_id, reads from filt_step2
    file("forwardP.fa") from forprimers
    file("reverseP.fa") from revprimers
    file("forwardP_rc.fa") from rcfor
    file("reverseP_rc.fa") from rcrev
    
    output:
    set val(fastq_id), val(fcid), val(sampleid), val(ext_id), val(pcr_id), "${fastq_id}*.R[12].cutadapt.fastq.gz" optional true into filt_step3
    file "*.cutadapt.out" into cutadaptToMultiQC
    file ('*demux.txt') into demux_to_validate

    script:
    if (params.lengthvar == false) 
        """
        #!/bin/bash
        cutadapt \\
            -g file:forwardP.fa \\
            -G file:reverseP.fa \\
            --cores ${task.cpus} \\
            -n 2 -e 1 \\
            --no-indels \\
            -o "${fastq_id}.{name}.R1.cutadapt.fastq.gz" \\
            -p "${fastq_id}.{name}.R2.cutadapt.fastq.gz" \\
            "${reads[0]}" "${reads[1]}" > "${fastq_id}.cutadapt.out"           

        # write out demuxed fastq names for sample sheet
        echo "\$(ls | grep .fastq.gz)" > "${fastq_id}_demux.txt"
        """
    else if (params.lengthvar == true)
        """
        cutadapt \\
            -g file:forwardP.fa -a file:reverseP_rc.fa \\
            -G file:reverseP.fa -a file:forwardP_rc.fa\\
            --cores ${task.cpus} \\
            -n 2 -e 1 \\
            --no-indels \\
            -o "${fastq_id}.{name}.R1.cutadapt.fastq.gz" \\
            -p "${fastq_id}.{name}.R2.cutadapt.fastq.gz" \\
            "${reads[0]}" "${reads[1]}" > "${fastq_id}.cutadapt.out"
            
        # write out demuxed fastq names for sample sheet
        echo "\$(ls | grep .fastq.gz)" > "${fastq_id}_demux.txt"
        """
}

/*
 *
 * Step 4: Filter reads
 *
 */
 
process filter_and_trim {
    tag { "filt_step3_${fastq_id}" }

    input:
    set fastq_id, fcid, sampleid, ext_id, pcr_id, file(reads), file(trimming) from filt_step3.join(filt_step3Trimming)

    output:
    set val(fastq_id), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" optional true into filteredReadsforQC, filteredReads   
    file("*.filtered.fastq.gz") into filtReadsErr
    file("*.filtered.fastq.gz") into filtReads
    file "*.trimmed.txt" into trimTracking

    script:
    """
    #!/usr/bin/env Rscript
    require(dada2); packageVersion("dada2")
    require(ShortRead); packageVersion("ShortRead")
    require(Biostrings); packageVersion("Biostrings")
    require(stringr); packageVersion("stringr")
    require(dplyr); packageVersion("dplyr")
    
    fastqFs <- sort(list.files(pattern="*.R1.cutadapt.fastq.gz"))
    fastqFs <- fastqFs[!stringr::str_detect(fastqFs, "unknown")]
    fastqRs <- sort(list.files(pattern="*.R2.cutadapt.fastq.gz"))
    fastqRs <- fastqRs[!stringr::str_detect(fastqRs, "unknown")]
    
    out1 <- readRDS("${trimming}")
    out2 <- filterAndTrim(
                        fwd = fastqFs,
                        filt = stringr::str_replace(fastqFs, ".cutadapt.fastq.gz", ".filtered.fastq.gz"),
                        rev = fastqRs,
                        filt.rev = stringr::str_replace(fastqRs, ".cutadapt.fastq.gz", ".filtered.fastq.gz"),
                        maxEE = c(${params.maxEEFor},${params.maxEERev}),
                        trimLeft = c(${params.trimFor},${params.trimRev}),
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
    out1 <- out1 %>% 
        as.data.frame() %>%
        dplyr::slice(rep(dplyr::row_number(), nrow(out2)))
    out3 <- cbind(out1, out2)      
    rownames(out3) <- rownames(out2)
    colnames(out3) <- c('cutadapt', 'filtered', 'input', 'filterN')    
    write.csv(out3, paste0("${fastq_id}", ".trimmed.txt"))
    """
}


/*
 *
 * Step 5: Post-filter Quality control
 *
 */

process fastqc_filtered {
    tag { "fastqc_filtered.${fastq_id}" }
    publishDir "${params.outdir}/qc/FastQC-postfilter", mode: "copy", overwrite: true

    input:
    set val(fastq_id), file(filtFor), file(filtRev) from filteredReadsforQC

    output:
    file '*_fastqc.{zip,html}' into fastqc_files_post

    """
    fastqc --nogroup -q ${filtFor} ${filtRev}
    """
}

process multiqc_filtered {
    tag { "multiqc_filtered" }
    publishDir "${params.outdir}/qc/MultiQC-postfilter", mode: 'copy', overwrite: true

    input:
    file('./raw-seq/*') from fastqc_files2_ch.collect()
    file('./trimmed-seq/*') from fastqc_files_post.collect()
    file('./cutadapt/*') from cutadaptToMultiQC.collect()

    output:
    file "*_report.html" into multiqc_report_post
    file "*_data"

    script:
    interactivePlots = params.interactiveMultiQC == true ? "-ip" : ""
    """
    multiqc ${interactivePlots} .
    """
}

// TODO: Make this work with demultiplexing
process merge_trimmed_table {
    tag { "merge_trimmed_table" }
    publishDir "${params.outdir}/csv", mode: "copy", overwrite: true

    input:
    file trimData from trimTracking.collect()

    output:
    file "all.trimmed.csv" into trimmed_read_tracking

    script:
    """
    #!/usr/bin/env Rscript
    require(tidyverse); packageVersion("tidyverse")
    readr::read_csv(list.files(path = '.', pattern = '*.trimmed.txt')) %>%
        magrittr::set_colnames(c("sample_id", "cutadapt", "filtered", "input", "filterN")) %>%
        dplyr::mutate(sample_id = sample_id %>% stringr::str_remove(".R[1-2].cutadapt.fastq.gz")) %>%
        write_csv("all.trimmed.csv")
    """
}

/*
 *
 * Step 6: Learn error rates (run on all samples)
 *
 */

// Make new channels with fcid as key
filtReadsErrKey = filtReadsErr.collect().flatten().map {
        file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }.groupTuple()
     
filtReadsKey = filtReads.collect().flatten().map {
        file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }.groupTuple()
     
process dada2_learn_errors {
    tag { "dada2_learn_errors" }
    publishDir "${params.outdir}/qc", mode: "copy", overwrite: true

    input:
    set key, file(reads) from filtReadsErrKey.view()

    output:
    file "errorsF_${key}.RDS" into errorsFor
    file "errorsR_${key}.RDS" into errorsRev
    file "*.pdf"
    
    script:
    """
    #!/usr/bin/env Rscript
    require(dada2); packageVersion("dada2")    
    setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})

    # File parsing
    filtFs <- list.files('.', pattern="${key}.*R1.filtered.fastq.gz", full.names = TRUE)
    filtRs <- list.files('.', pattern="${key}.*R2.filtered.fastq.gz", full.names = TRUE)
    
    #TODO Change this to new param: dada2 seed?
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

    pdf("err_${key}.pdf")
    plotErrors(errF, nominalQ=TRUE)
    plotErrors(errR, nominalQ=TRUE)
    dev.off()

    saveRDS(errF, "errorsF_${key}.RDS")
    saveRDS(errR, "errorsR_${key}.RDS")
    """
}

/*
 *
 * Step 7: Dereplication, ASV Inference, Merge Pairs
 *
 */

if (params.pool == "T" || params.pool == 'pseudo') {

    process dada2_pooled {
        tag { "dada2_pooled" }
        publishDir "${params.outdir}/qc/dada2", mode: "copy", overwrite: true

        input:
        set key, file(reads) from filtReadsKey.view()
        file(errFor) from errorsFor
        file(errRev) from errorsRev

        output:
        file "seqtab.RDS" into seqTable,rawSeqTableToRename
        file "*.mergers.RDS" into merged_read_tracking
        file "*.dadaFs.RDS" into dada_for_read_tracking
        file "*.dadaRs.RDS" into dada_rev_read_tracking
        file "seqtab.*"

        script:
        """
        #!/usr/bin/env Rscript
        require(dada2); packageVersion("dada2")
        require(tidyverse); packageVersion("tidyverse")
        setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})
        filtFs <- list.files('.', pattern="${key}.*R1.filtered.fastq.gz", full.names = TRUE)
        filtRs <- list.files('.', pattern="${key}.*R2.filtered.fastq.gz", full.names = TRUE)

        errF <- readRDS("${errFor}")
        errR <- readRDS("${errRev}")
        cat("Processing all samples")

        #Variable selection from CLI input flag --pool
        pool <- "${params.pool}"
        if(pool == "T" || pool == "TRUE"){
          pool <- as.logical(pool)
        }
        dadaFs <- dada(filtFs, err=errF, multithread=${task.cpus}, pool=pool)
        dadaRs <- dada(filtRs, err=errR, multithread=${task.cpus}, pool=pool)

        mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
            returnRejects = TRUE,
            minOverlap = ${params.minOverlap},
            maxMismatch = ${params.maxMismatch},
            trimOverhang = as.logical("${params.trimOverhang}"),
            justConcatenate = as.logical("${params.justConcatenate}")
            )

        # TODO: make this a single item list with ID as the name, this is lost
        # further on
        saveRDS(mergers, "${key}.mergers.RDS")

        saveRDS(dadaFs, "${key}.dadaFs.RDS")
        saveRDS(dadaRs, "${key}.dadaRs.RDS")

        # go ahead and make seqtable
        seqtab <- makeSequenceTable(mergers)
        saveRDS(seqtab, "seqtab.RDS")
        
        # Track reads
        getN <- function(x){sum(getUniques(x))}
        dada_out <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN)) %>%
            magrittr::set_colnames(c("dadaFs", "dadaRs", "merged")) %>%
            as.data.frame()
            
        # TODO: change this to flow cell id
        write.csv(dada_out, "dada_${key}.csv") 
        """
        }
} else {
    // pool = F, process per sample
    process dada2_persample {
        tag { "dada2_persample" }
        publishDir "${params.outdir}/qc/dada2", mode: "copy", overwrite: true

        input:
        set val(fastq_id), file(filtFor), file(filtRev) from filteredReads
        file errFor from errorsFor
        file errRev from errorsRev

        output:
        file "seqtab.RDS" into seqTable
        file "all.mergers.RDS" into merged_read_tracking
        file "all.dadaFs.RDS" into dada_for_read_tracking
        file "all.dadaRs.RDS" into dada_rev_read_tracking
        file "seqtab.*"

        script:
        """
        #!/usr/bin/env Rscript
        require(dada2); packageVersion("dada2")        
        setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})

        errF <- readRDS("${errFor}")
        errR <- readRDS("${errRev}")
        cat("Processing:", "${fastq_id}", "\\n")
        
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

        saveRDS(merger, paste("${fastq_id}", "merged", "RDS", sep="."))

        saveRDS(dadaFs, "all.dadaFs.RDS")
        saveRDS(dadaRs, "all.dadaRs.RDS")
        """
    }

    process merge_dada_rds {
        tag { "merge_dada_rds" }
        publishDir "${params.outdir}/qc/dada2", mode: "copy", overwrite: true

        input:
        file dadaFs from dadaFor.collect()
        file dadaRs from dadaRev.collect()

        output:
        file "all.dadaFs.RDS" into dada_for_read_tracking
        file "all.dadaRs.RDS" into dada_rev_read_tracking

        script:
        '''
        #!/usr/bin/env Rscript
        require(dada2); packageVersion("dada2")

        dadaFs <- lapply(list.files(path = '.', pattern = '.dadaFs.RDS$'), function (x) readRDS(x))
        names(dadaFs) <- sub('.dadaFs.RDS', '', list.files('.', pattern = '.dadaFs.RDS'))
        dadaRs <- lapply(list.files(path = '.', pattern = '.dadaRs.RDS$'), function (x) readRDS(x))
        names(dadaRs) <- sub('.dadaRs.RDS', '', list.files('.', pattern = '.dadaRs.RDS'))
        saveRDS(dadaFs, "all.dadaFs.RDS")
        saveRDS(dadaRs, "all.dadaRs.RDS")
        '''
    }

    process make_seqtab {
        tag { "make_seqtab" }
        publishDir "${params.outdir}/qc/dada2", mode: "copy", overwrite: true

        input:
        file mr from mergedReads.collect()

        output:
        file "seqtab.RDS" into seqTable,rawSeqTableToRename
        file "all.mergers.RDS" into merged_read_tracking

        script:
        '''
        #!/usr/bin/env Rscript
        require(dada2); packageVersion("dada2")        
        
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


process asv_filter {
    tag { "asv_filter" }
    publishDir "${params.outdir}/rds", pattern: 'seqtab_final.rds', mode: "copy", overwrite: true
    publishDir "${params.outdir}/qc", pattern: 'seqtab_length_dist.pdf', mode: "copy", overwrite: true
    publishDir "${params.outdir}/qc", pattern: 'ASV_cleanup_summary.csv', mode: "copy", overwrite: true

    input:
    file st from seqTable
    file phmm from phmm_ch.ifEmpty( "" )

    output:
    file "seqtab_final.RDS" into seqtab_to_tax,seqtab_to_exact,seqtab_to_blast,seqtab_to_rename,seqtab_to_output,seqtab_read_tracking
    file "ASV_cleanup_summary.csv" 
    file "seqtab_length_dist.pdf"
    file "sample_cleanup_summary.csv" into asv_read_tracking

    script:
    chimOpts = params.removeBimeraDenovoOptions != false ? ", ${params.removeBimeraDenovoOptions}" : ''
    
    if (params.coding == true) 
        """
        #!/usr/bin/env Rscript
        require(dada2); packageVersion("dada2")        
        require(tidyverse); packageVersion("tidyverse")
        require(Biostrings); packageVersion("Biostrings")
        require(taxreturn); packageVersion("taxreturn")
        require(patchwork); packageVersion("patchwork")
        st.all <- readRDS("${st}")

        # Remove chimeras
        seqtab_nochim <- removeBimeraDenovo(
            st.all, 
            method="consensus", 
            multithread=${task.cpus}, 
            verbose=TRUE ${chimOpts} 
            )

        #cut to expected size
        min_asv_len <- as.numeric(${params.min_asv_len})
        max_asv_len <- as.numeric(${params.max_asv_len})
        
        if(min_asv_len > 0 && max_asv_len > 0){
            seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) %in% as.numeric(${params.min_asv_len}):as.numeric(${params.max_asv_len})]
        } else {
            seqtab_cut <- seqtab_nochim
        }
        
        seqs <- Biostrings::DNAStringSet(getSequences(seqtab_cut))
        names(seqs) <- getSequences(seqtab_cut)
        
        # Align against phmm if provided    
        if(file.info("${phmm}")\$size > 0){
            model <- readRDS("${phmm}")
            phmm_filt <- taxreturn::map_to_model(seqs, model = model, min_score = 100, min_length = 100, shave = FALSE, check_frame = TRUE, kmer_threshold = 0.5, k=5, extra = "fill")
            seqs <- Biostrings::DNAStringSet(names(phmm_filt))
            names(seqs) <- names(phmm_filt)
        }
        
        #Filter sequences containing stop codons        
        codon_filt <- seqateurs::codon_filter(seqs, genetic.code = "${params.genetic_code}")
        seqtab_final <- seqtab_cut[,colnames(seqtab_cut) %in% names(codon_filt)]                
        
        saveRDS(seqtab_final, "seqtab_final.RDS")
        
        # Summarise reads retained per sample
        asv_read_tracker <- tibble::enframe(rowSums(seqtab_nochim), name="sample_id", value="chimera_filt") %>%
            left_join(tibble::enframe(rowSums(seqtab_cut), name="sample_id", value="asv_size_filt")) %>%
            left_join(tibble::enframe(rowSums(seqtab_final), name="sample_id", value="asv_codon_filt"))

        if(file.info("${phmm}")\$size > 0){
            asv_read_tracker <- asv_read_tracker %>%
                left_join(tibble::enframe(rowSums(seqtab_cut[,colnames(seqtab_cut) %in% names(seqs)]), name="sample_id", value="asv_phmm_filt")) %>%
                dplyr::select(sample_id, chimera_filt, asv_size_filt, asv_phmm_filt, asv_codon_filt)
        }
        write_csv(asv_read_tracker, "sample_cleanup_summary.csv")
        
        # Track fate of ASV's
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
    else if (params.coding == false) 
        """
        #!/usr/bin/env Rscript
        require(dada2); packageVersion("dada2")        
        require(tidyverse); packageVersion("tidyverse")
        require(Biostrings); packageVersion("Biostrings")
        require(patchwork); packageVersion("patchwork")
        
        st.all <- readRDS("${st}")

        # Remove chimeras
        seqtab_nochim <- removeBimeraDenovo(
            st.all, 
            method="consensus", 
            multithread=${task.cpus}, 
            verbose=TRUE ${chimOpts} 
            )

        #cut to expected size
        min_asv_len <- as.numeric(${params.min_asv_len})
        max_asv_len <- as.numeric(${params.max_asv_len})
        
        if(min_asv_len > 0 && max_asv_len > 0){
            seqtab_final <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) %in% as.numeric(${params.min_asv_len}):as.numeric(${params.max_asv_len})]
        } else {
            seqtab_final <- seqtab_nochim
        }       
        
        saveRDS(seqtab_final, "seqtab_final.RDS")       

        # Summarise reads retained per sample
        asv_read_tracker <- tibble::enframe(rowSums(seqtab_nochim), name="sample_id", value="chimera_filt") %>%
            left_join(tibble::enframe(rowSums(seqtab_final), name="sample_id", value="asv_size_filt"))

        write_csv(asv_read_tracker, "sample_cleanup_summary.csv")
        
        # Track fate of ASV's
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


/*
 *
 * Step 9: Taxonomic assignment
 *
 */

if (params.reference) {
    ref_file = file(params.reference)

    if (params.taxassignment == 'rdp') {
        // TODO: we could combine these into the same script
        
        process assign_tax_rdp {
            tag { "assign_tax_rdp" }
            publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

            input:
            file st from seqtab_to_tax
            file ref from ref_file

            output:
            file "tax_rdp.RDS" into tax_to_combine
            file "bootstrap_final.RDS" into bootstrapFinal

            script:
            """
            #!/usr/bin/env Rscript
            require(dada2); packageVersion("dada2")             
            
            if (!stringr::str_detect("${ref_file}", ".fasta|.fa|.fa.gz|.fasta.gz|.fas")){
                stop("reference file for RDP must be a fasta file")
            }
            seqtab <- readRDS("${st}")

            # Assign taxonomy
            tax <- assignTaxonomy(seqtab, "${ref}",
                                    multithread=${task.cpus},
                                    tryRC = TRUE,
                                    outputBootstraps = TRUE,
                                    minBoot = ${params.minBoot},
                                    verbose = TRUE)
            boots <- tax\$boot
                         
            # Write original data
            saveRDS(tax, "tax_rdp.RDS")
            saveRDS(boots, "bootstrap_final.RDS")
            """
        }
        
    } else if (params.taxassignment == 'idtaxa') {
        process assign_tax_idtaxa {
            tag { "assign_tax_idtaxa" }
            publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

            input:
            file st from seqtab_to_tax
            file ref from ref_file

            output:
            file "tax_idtaxa.RDS" into tax_to_combine
            file "bootstrap_final.RDS" into bootstrapFinal
            file "raw_idtaxa.RDS"

            script:
            """
            #!/usr/bin/env Rscript
            require(dada2); packageVersion("dada2")
            require(DECIPHER); packageVersion("DECIPHER")
            require(tidyverse); packageVersion("tidyverse")
            require(seqateurs); packageVersion("seqateurs")
            require(stringr); packageVersion("stringr")
            require(stringi); packageVersion("stringi")
            
            seqtab <- readRDS("${st}")

            # Create a DNAStringSet from the ASVs
            dna <- DNAStringSet(getSequences(seqtab))

            # load database; this should be a RData file
            if (stringr::str_detect("${ref_file}", ".RData")){
                load("${ref_file}")
            } else if(stringr::str_detect("${ref_file}", ".rds")){
                trainingSet <- readRDS("${ref_file}")
            } else {
                stop("reference file for IDTAXA must be in .RData or .RDS format")
            }

            ids <- IdTaxa(dna, trainingSet,
                strand="top",
                processors=${task.cpus},
                verbose=TRUE)
                
            # ranks of interest
            ranks <-  c("root", "kingdom", "phylum","class", "order", "family", "genus","species") 
            saveRDS(ids, 'raw_idtaxa.RDS')

            #Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
            tax <- t(sapply(ids, function(x) {
              taxa <- paste0(x\$taxon,"_", x\$confidence)
              taxa[startsWith(taxa, "unclassified_")] <- NA
              taxa
            })) %>%
              purrr::map(unlist) %>%
              stri_list2matrix(byrow=TRUE, fill=NA) %>%
              magrittr::set_colnames(ranks) %>%
              as.data.frame() %>%
              magrittr::set_rownames(getSequences(seqtab)) %>%
              mutate_all(str_replace,pattern="(?:.(?!_))+\$", replacement="")

            boots <- t(sapply(ids, function(x) {
                    m <- match(ranks, x\$rank)
                    bs <- x\$confidence[m]
                    bs
            }))
            colnames(boots) <- ranks
            rownames(boots) <- getSequences(seqtab)

            # Write to disk
            saveRDS(tax, "tax_idtaxa.RDS")
            saveRDS(boots, "bootstrap_final.RDS")
            """
        }

    } else if (params.taxassignment) {
        exit 1, "Unknown taxonomic assignment method set: ${params.taxassignment}"
    } else {
        exit 1, "No taxonomic assignment method set, but reference passed"
    }
    
    if (params.species_db ){
      species_file = file(params.species_db)
      
        process assign_tax_exact {
            tag { "assign_tax_exact" }
            publishDir "${params.outdir}/rds", mode: "copy", overwrite: true
            
            input:
            file st from seqtab_to_exact
            file ref from species_file
            
            output:
            file "tax_exact.RDS" into exact_to_combine
            
            script:
            """
            #!/usr/bin/env Rscript
            require(dada2); packageVersion("dada2")
            require(tidyverse); packageVersion("tidyverse")
            
            seqtab <- readRDS("${st}")
            
            # Assign species using exact matching
            tax_exact <- assignSpecies(seqtab, "${ref}", allowMultiple = TRUE, tryRC = TRUE, verbose = FALSE) %>%
                as_tibble(rownames = "OTU") %>%
                filter(!is.na(Species)) %>%
                dplyr::mutate(binomial = paste0(Genus," ",Species)) %>%
                dplyr::rename(exact_genus = Genus, exact_species = Species)
            
            # Write to disk
            saveRDS(tax_exact, "tax_exact.RDS")
            """
            
           }
    } else {
        exact_to_combine = Channel.empty()
    }
    
    if (params.blast_db ){
      blast_file = file(params.blast_db)
      
        process assign_tax_blast {
            tag { "assign_tax_blast" }
            publishDir "${params.outdir}/rds", mode: "copy", overwrite: true
            
            input:
            file st from seqtab_to_blast
            file ref from blast_file
            
            output:
            file "tax_blast.RDS" into blast_to_combine
            
            script:
            """
            #!/usr/bin/env Rscript
            require(dada2); packageVersion("dada2")
            require(DECIPHER); packageVersion("DECIPHER")
            require(tidyverse); packageVersion("tidyverse")
            require(seqateurs); packageVersion("seqateurs")
            require(taxreturn); packageVersion("taxreturn")
            
            seqtab <- readRDS("${st}")
            
            # Add species using BLAST
            seqs <- taxreturn::char2DNAbin(colnames(seqtab))
            names(seqs) <- colnames(seqtab) 
            
            # TODO: allow editing of these BLAST parameters
            tax_blast <- taxreturn::blast_assign_species(query=seqs, db="${ref}", identity=97, coverage=95, evalue=1e06, max_target_seqs=5, max_hsp=5,
                ranks=c("Root", "Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"), delim=";") %>%
                dplyr::rename(blast_genus = Genus, blast_spp = Species) %>%
                dplyr::filter(!is.na(blast_spp))
            
            # Write to disk
            saveRDS(tax_blast, "tax_blast.RDS")            
            """
            
           }
    } else {
        blast_to_combine = Channel.empty()
    }
    
    process combine_tax {
            tag { "combine_tax" }
            publishDir "${params.outdir}/rds", mode: "copy", overwrite: true

            input:
            file(tax) from tax_to_combine
            file(exact) from exact_to_combine.ifEmpty( "" )
            file(blast) from blast_to_combine.ifEmpty( "" )
            
            output:
            file "tax_final.RDS" into taxtab_to_output

            script:
            """
            #!/usr/bin/env Rscript
            require(tidyverse); packageVersion("tidyverse")
            require(seqateurs); packageVersion("seqateurs")
            
            tax <- readRDS("${tax}") %>%
                  as_tibble(rownames = "OTU")
            
            #Join tax with blast
            if(file.info("${blast}")\$size > 0){
                blast_spp <- readRDS("${blast}")
                tax <- tax  %>%
                  left_join(blast_spp , by="OTU") %>%
                  dplyr::mutate(species = case_when(
                    is.na(species) & genus == blast_genus ~ blast_spp,
                    !is.na(species) ~ species
                  )) 
            }
                  
            #Join tax with exact matching
            if(file.info("${exact}")\$size > 0){
                exact <- readRDS("${exact}")
                tax <- tax %>%
                  left_join(exact, by="OTU") %>%
                  dplyr::mutate(species = case_when(
                    is.na(species) & genus == exact_genus ~ binomial,
                    !is.na(species) ~ species
                  )) 
            }
            
            # Make final tax object
            ranks <-  c("root", "kingdom", "phylum","class", "order", "family", "genus","species") 
            tax_final <- tax %>%
                dplyr::select(OTU, all_of(ranks)) %>%
                column_to_rownames("OTU") %>%
                seqateurs::na_to_unclassified() %>% #Propagate high order ranks to unassigned ASVs
                as.matrix()
             
            # Write to disk
            saveRDS(tax_final, "tax_final.RDS")
            """
    }
     
} else {
    // set tax channels to 'false', do NOT assign taxonomy
    taxtab_to_output = Channel.empty()
    bootstrapFinal = Channel.empty()
}

/*
 *
 * Step 10: Alignment & Phylogenetic tree
 *
 */
 
process output_asvs {
    tag { "output_asvs" }

    input:
    file st from seqtab_to_rename
    file rawst from rawSeqTableToRename

    output:
    file "asvs.${params.idType}.nochim.fna" into seqsToAln

    script:
    """
    #!/usr/bin/env Rscript
    require(dada2)
    require(Biostrings)    
    require(digest)
    
    # read RDS w/ data
    st <- readRDS("${st}")
    
    # get sequences
    seqs <- colnames(st)
   
    # generate FASTA
    seqs.dna <- Biostrings::DNAStringSet(seqs)
    names(seqs.dna) <- colnames(st)
    # Write out fasta file
    Biostrings::writeXStringSet(seqs.dna, filepath = 'asvs.${params.idType}.nochim.fna')
    """
}

// TODO: if PHMM provided, subset the previous PHMM alignment to just the retained sequences and output that?
if (params.runTree && params.lengthvar == false) {

    if (params.aligner == 'DECIPHER') {

        process align_seqs_decipher {
            tag { "align_seqs_decipher" }
            publishDir "${params.outdir}/fasta", mode: "copy", overwrite: true
            errorStrategy 'ignore'

            input:
            file seqs from seqsToAln

            output:
            file "aligned_seqs.fasta" optional true into alnFile,aln_to_output
            
            script:
            """
            #!/usr/bin/env Rscript
            require(dada2); packageVersion("dada2")
            require(DECIPHER); packageVersion("DECIPHER")

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

        process make_tree_phangorn {
            tag { "make_tree_phangorn" }
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
            require(phangorn); packageVersion("phangorn")

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

        process make_tree_fasttree {
            tag { "make_tree_fasttree" }
            publishDir "${params.outdir}/trees", mode: "copy", overwrite: true

            input:
            file aln from alnFile

            output:
            file "fasttree.tree" into treeGTRFile
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

    process root_tree {
        tag { "root_tree" }
        publishDir "${params.outdir}/trees", mode: "link"

        input:
        file tree from treeGTRFile

        output:
        file "rooted.newick" into tree_to_output

        script:
        """
        #!/usr/bin/env Rscript
        require(phangorn); packageVersion("phangorn")
        require(ape); packageVersion("ape")

        tree <- read.tree(file = "${tree}")

        midtree <- midpoint(tree)

        write.tree(midtree, file = "rooted.newick")
        """
    }
} else {
    // Note these are caught downstream
    aln_to_output = false
    tree_to_output = false
}


/*
 *
 * Step 12: Track reads
 *
 */

// TODO: add the results of the seqtable tracking
process track_reads {
    tag { "track_reads" }
    publishDir "${params.outdir}/qc", mode: "copy", overwrite: true

    input:
    file(trimmed) from trimmed_read_tracking
    file(mergers) from merged_read_tracking
    file(dadaFs) from dada_for_read_tracking
    file(dadaRs) from dada_rev_read_tracking
    file(st) from seqtab_read_tracking
    file(asv) from asv_read_tracking

    output:
    file "all.readtracking.tsv"

    script:
    """
    #!/usr/bin/env Rscript
    require(dada2); packageVersion("dada2")
    require(tidyverse); packageVersion("tidyverse")

    getN <- function(x) sum(getUniques(x))

    # the gsub here might be a bit brittle...
    dadaFs <- as.data.frame(sapply(readRDS("${dadaFs}"), getN))
    rownames(dadaFs) <- stringr::str_remove(rownames(dadaFs), '.R1.filtered.fastq.gz')
    colnames(dadaFs) <- c("denoisedF")
    dadaFs\$sample_id <- rownames(dadaFs)

    dadaRs <- as.data.frame(sapply(readRDS("${dadaRs}"), getN))
    rownames(dadaRs) <- stringr::str_remove(rownames(dadaRs), '.R2.filtered.fastq.gz')
    colnames(dadaRs) <- c("denoisedR")
    dadaRs\$sample_id <- rownames(dadaRs)

    all.mergers <- readRDS("${mergers}")
    mergers <- as.data.frame(sapply(all.mergers, function(x) sum(getUniques(x %>% dplyr::filter(accept)))))
    rownames(mergers) <- stringr::str_remove(rownames(mergers), '.R1.filtered.fastq.gz')
    colnames(mergers) <- c("merged")
    mergers\$sample_id <- rownames(mergers)

    trimmed <- read.csv("${trimmed}")
    
    asv_cleanup <- read.csv("${asv}") %>%
        dplyr::mutate(sample_id = stringr::str_remove(sample_id, '.R1.filtered.fastq.gz'))

    track <- Reduce(function(...) merge(..., by = "sample_id",  all.x=TRUE),  list(trimmed, dadaFs, dadaRs, mergers, asv_cleanup))
    # dropped data in later steps gets converted to NA on the join
    # these are effectively 0
    track[is.na(track)] <- 0

    write_tsv(track, "all.readtracking.tsv")
    """
}

/*
 *
 * Create and check samplesheet
 *
 */
 
// TODO: need to have input fastq list from demultiplexed files
process create_samdf {
    tag { "create_samdf" }
    publishDir "${params.outdir}/sample_info", mode: "copy", overwrite: true

    input:
    file(samplesheet) from samplesheet_ch
    file(runparams) from runparams_ch
    file(fastq_names) from samples_to_validate.collectFile(name: 'fastq_list.txt', newLine: true)
    file(fastq_names) from demux_to_validate.collectFile(name: 'demux_list.txt', newLine: true)
    
    output:
    file "*.csv" into samdf_to_output
    
    script:
    """
    #!/usr/bin/env Rscript
    require(seqateurs)
    require(tidyverse)
    SampleSheet <- normalizePath( "${samplesheet}")
    runParameters <- normalizePath( "${runparams}")
    
    # Create samplesheet containing samples and run parameters for all runs
    samdf <- dplyr::distinct(seqateurs::create_samplesheet(SampleSheet = SampleSheet, runParameters = runParameters, template = "V4"))

    # Add fcid if not present
    samdf <- samdf %>%
      dplyr::mutate(
      sample_id = case_when(
        !stringr::str_detect(sample_id, fcid) ~ paste0(fcid,"_", sample_id),
        TRUE ~ sample_id
      ))

    # Check if samples match samplesheet
    fastqFs <- read_delim("fastq_list.txt", delim=" ", col_names=c("sample_id", "fcid", "sample_name", "extraction_rep", "amp_rep")) %>%
        dplyr::filter(!stringr::str_detect(sample_id, "Undetermined")) %>%
        dplyr::mutate(sample_id = str_remove(sample_id, pattern = "(?:.(?!_S))+\$"))
    
    #Check missing in samplesheet
    if (length(setdiff(fastqFs\$sample_id, samdf\$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(fastqFs\$sample_id, samdf\$sample_id), " are not in the sample sheet") }

    #Check missing fastqs
    if (length(setdiff(samdf\$sample_id, fastqFs\$sample_id)) > 0) {
      samdf <- samdf %>%
        filter(!sample_id %in% setdiff(samdf\$sample_id, fastqFs\$sample_id))
    }

    # Add mising fields
    primer_names <- purrr::map2(unlist(str_split("${params.fwdprimer_name}", ";")), unlist(str_split("${params.revprimer_name}", ";")), ~{
            paste0(.x, "-", .y)
        }) %>%
        unlist() %>%
        paste0(collapse=";")
    
    samdf <- samdf %>%
      dplyr::mutate(
      for_primer_seq = "${params.fwdprimer}",
      rev_primer_seq = "${params.revprimer}",
      pcr_primers = primer_names,
      sample_name = NA_character_
      ) %>%
      seqateurs::coalesce_join(fastqFs, by="sample_id")

    #Update the sample sheet and logging sheet to deal with any newly demultiplexed files
    samdf <- samdf %>%
     group_by(sample_id) %>%
     group_split() %>%
     purrr::map(function(x){
       #if(any(str_detect(x\$pcr_primers, ";"))){
         primer_names <- unlist(str_split(unique(x\$pcr_primers), ";"))
         x %>%
           mutate(count = length(primer_names)) %>% #Replicate the samples
           uncount(count) %>%
           mutate(pcr_primers = unlist(str_split(unique(x\$pcr_primers), ";")),
                  for_primer_seq = unlist(str_split(unique(x\$for_primer_seq), ";")),
                  rev_primer_seq = unlist(str_split(unique(x\$rev_primer_seq), ";")),
                  sample_id = paste0(sample_id, "_", pcr_primers)
                )
       #} else (x)
     }) %>%
      bind_rows()

    # Check f these samples are present in the demultiplexed fastq list
    demux_fastqFs <- readLines("demux_list.txt")  %>%
        stringr::str_remove(".R[1-2].cutadapt.fastq.gz")%>%
        stringr::str_replace("_S[0-9].*\\\\.", "_")
    demux_fastqFs <- demux_fastqFs[!demux_fastqFs == ""]
    demux_fastqFs <- demux_fastqFs[!stringr::str_detect(demux_fastqFs, "_unknown\$")]

    #Check missing in samplesheet
    if (length(setdiff(demux_fastqFs, samdf\$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(demux_fastqFs, samdf\$sample_id), " are not in the sample sheet") }

    #Check missing fastqs
    if (length(setdiff(samdf\$sample_id, demux_fastqFs)) > 0) {
      samdf <- samdf %>%
        filter(!sample_id %in% setdiff(samdf\$sample_id, demux_fastqFs))
    }
    
    #Write out updated sample CSV for use
    write_csv(samdf, "Sample_info.csv")
    """
}
 
 
/*
 *
 * Step 12: Generate outputs
 *
 */
 
process make_phyloseq {
    tag { "make_phyloseq" }
    publishDir "${params.outdir}/results/unfiltered", mode: "link", overwrite: true

    input:
    file st from seqtab_to_output
    file samdf from samdf_to_output
    file tax from taxtab_to_output
    file bt from bootstrapFinal
    file tree from tree_to_output
    file aln from aln_to_output
    
    output:
    file "ps.rds" into ps_to_filter,ps_to_export
    
    script:
    if (tree_to_output != false)
        """
        #!/usr/bin/env Rscript
        require(phyloseq); packageVersion("phyloseq")
        require(tidyverse); packageVersion("tidyverse")
        require(ape); packageVersion("ape")
                
        # Read in files
        seqtab <- readRDS("${st}")
        
        #Extract start of sample names only
        rownames(seqtab) <- rownames(seqtab) %>%
                stringr::str_remove(".R[1-2]\\\\..*\\\\.fastq.gz") %>%
                stringr::str_replace("_S[0-9].*\\\\.", "_")
                
        tax <- readRDS("${tax}")
        colnames(tax) <- stringr::str_to_lower(colnames(tax))
        seqs <- Biostrings::readDNAStringSet("${aln}")
        tree <- read.tree(file = "${tree}")

        samdf <- read.csv("${samdf}", header=TRUE) %>%
          filter(!duplicated(sample_id)) %>%
          magrittr::set_rownames(.\$sample_id) 
        
        # Create phyloseq object
        ps <- phyloseq(tax_table(tax), 
                       sample_data(samdf),
                       otu_table(seqtab, taxa_are_rows = FALSE),
                       phy_tree(tree),
                       refseq(seqs))
                       
        saveRDS(ps, "ps.rds")
        """
    else
        """
        #!/usr/bin/env Rscript
        require(phyloseq); packageVersion("phyloseq")
        require(tidyverse); packageVersion("tidyverse")
        require(ape); packageVersion("ape")
        
        # Read in files
        seqtab <- readRDS("${st}")
        
        #Extract start of sample names only
        rownames(seqtab) <- rownames(seqtab) %>%
                stringr::str_remove(".R[1-2]\\\\..*\\\\.fastq.gz") %>%
                stringr::str_replace("_S[0-9].*\\\\.", "_")
                
        tax <- readRDS("${tax}")
        colnames(tax) <- stringr::str_to_lower(colnames(tax))
        
        seqs <- Biostrings::DNAStringSet(colnames(seqtab))
        names(seqs) <- colnames(seqtab)

        samdf <- read.csv("${samdf}", header=TRUE) %>%
          filter(!duplicated(sample_id)) %>%
          magrittr::set_rownames(.\$sample_id) 
        
        # Create phyloseq object
        ps <- phyloseq(tax_table(tax), 
                       sample_data(samdf),
                       otu_table(seqtab, taxa_are_rows = FALSE),
                       refseq(seqs))
        
        saveRDS(ps, "ps.rds")
        """
}

// TODO: how to output pdf to  figs?
process filter_phyloseq {
    tag { "filter_phyloseq" }
    publishDir "${params.outdir}/results/filtered", mode: "link", overwrite: true

    input:
    file ps from ps_to_filter
    
    output:
    file "ps_filtered.rds" into tax_check_ala,tax_check_afd,ps_to_export2
    file "rarefaction.pdf"

    script:
    """
    #!/usr/bin/env Rscript
    require(phyloseq); packageVersion("phyloseq")
    require(tidyverse); packageVersion("tidyverse")
    require(ape); packageVersion("ape")
    require(vegan); packageVersion("vegan")
    
    # Taxonomic filters
    ps0 <- readRDS("${ps}")
    
    # Phylum - Change all to lowercase at the taxtable step
    if(nchar(as.character("${params.phylum}")) > 0){
        ps0 <- phyloseq::subset_taxa(ps0, phylum %in% unlist(stringr::str_split(as.character("${params.phylum}"), ";",  n=Inf)))
    }
    # Order
    if(nchar(as.character("${params.order}")) > 0){
        ps0 <- phyloseq::subset_taxa(ps0, order %in% unlist(stringr::str_split(as.character("${params.order}"), ";",  n=Inf)))
    }
    # Family
    if(nchar(as.character("${params.family}")) > 0){
        ps0 <- phyloseq::subset_taxa(ps0, family %in% unlist(stringr::str_split(as.character("${params.family}"), ";",  n=Inf)))
    }
    # Family
    if(nchar(as.character("${params.family}")) > 0){
        ps0 <- phyloseq::subset_taxa(ps0, family %in% unlist(stringr::str_split(as.character("${params.family}"), ";",  n=Inf)))
    }
    # Genus
    if(nchar(as.character("${params.genus}")) > 0){
        ps0 <- phyloseq::subset_taxa(ps0, genus %in% unlist(stringr::str_split(as.character("${params.genus}"), ";",  n=Inf)))
    }
    # Species
    if(nchar(as.character("${params.species}")) > 0){
        ps0 <- phyloseq::subset_taxa(ps0, species %in% unlist(stringr::str_split(as.character("${params.species}"), ";",  n=Inf)))
    }
    
    # Drop missing taxa
    ps0 <- ps0 %>%
      phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE) %>%
      phyloseq::prune_samples(sample_sums(.) >0, .) 

    #Set a threshold for minimum reads per sample
    threshold <- as.numeric("${params.min_sample_reads}")
    
    #Create rarefaction curve
    rare <- phyloseq::otu_table(ps0) %>%
      as("matrix") %>%
      vegan::rarecurve(step=max(sample_sums(ps0))/100) %>%
      purrr::map(function(x){
        b <- as.data.frame(x)
        b <- data.frame(OTU = b[,1], count = rownames(b))
        b\$count <- as.numeric(gsub("N", "",  b\$count))
        return(b)
      }) %>%
      purrr::set_names(sample_names(ps0)) %>%
      dplyr::bind_rows(.id="sample_id")

    # Write out rarefaction curve
    pdf(file="rarefaction.pdf", width = 11, height = 8 , paper="a4r")
    ggplot(data = rare)+
      geom_line(aes(x = count, y = OTU, group=sample_id), alpha=0.5)+
      geom_point(data = rare %>% 
                   group_by(sample_id) %>% 
                   top_n(1, count),
                 aes(x = count, y = OTU, colour=(count > threshold))) +
      geom_label(data = rare %>% 
                   group_by(sample_id) %>% 
                   top_n(1, count),
                 aes(x = count, y = OTU,label=sample_id, colour=(count > threshold)),
                 hjust=-0.05)+
      scale_x_continuous(labels =  scales::scientific_format()) +
      geom_vline(xintercept=threshold, linetype="dashed") +
      labs(colour = "Sample kept?") +
      xlab("Sequence reads") +
      ylab("Observed ASV's")
    try(dev.off(), silent=TRUE)

    #Remove all samples under the minimum read threshold 
    ps1 <- ps0 %>%
      phyloseq::prune_samples(sample_sums(.)>=threshold, .) %>% 
      phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
     
    # output filtered phyloseq object
    saveRDS(ps1, "ps_filtered.rds") 
    """
}

process output_unfiltered {
    tag { "output_unfiltered" }
    publishDir "${params.outdir}/results/unfiltered", mode: "link", overwrite: true

    input:
    file ps from ps_to_export
    
    output:
    file "*.csv"
    file "*.fasta"
    file "*.nwk"

    script:
    """
    #!/usr/bin/env Rscript
    require(phyloseq); packageVersion("phyloseq")
    require(seqateurs); packageVersion("seqateurs")
    require(tidyverse); packageVersion("tidyverse")
    require(ape); packageVersion("ape")
    
    ps <- readRDS("${ps}")
    
    #Export raw csv
    speedyseq::psmelt(ps) %>%
      dplyr::filter(Abundance > 0) %>%
      dplyr::select(-Sample) %>%
      write_csv("raw_combined.csv")
  
    #Export species level summary
    seqateurs::summarise_taxa(ps, "species", "sample_id") %>%
      tidyr::spread(key="sample_id", value="totalRA") %>%
      write.csv(file = "spp_sum_unfiltered.csv")
      
    #Export genus level summary
    seqateurs::summarise_taxa(ps, "genus", "sample_id") %>%
      tidyr::spread(key="sample_id", value="totalRA") %>%
      write.csv(file = "gen_sum_unfiltered.csv")

    #Output newick tree
    if(!is.null(phy_tree(ps, errorIfNULL=FALSE))){
        write.tree(phy_tree(ps), file="tree_filtered.nwk")
    } else {
        cat("", file = "tree_filtered.nwk")
    }

    #Output fasta of all ASV's
    seqateurs::ps_to_fasta(ps, out.file = "asvs_unfiltered.fasta", seqnames = "species")
    """
}

process output_filtered {
    tag { "output_filtered" }
    publishDir "${params.outdir}/results/filtered", mode: "link", overwrite: true

    input:
    file ps from ps_to_export2
    
    output:
    file "*.csv"
    file "*.fasta"
    file "*.nwk"

    script:
    """
    #!/usr/bin/env Rscript
    require(phyloseq); packageVersion("phyloseq")
    require(seqateurs); packageVersion("seqateurs")
    require(tidyverse); packageVersion("tidyverse")
    require(ape); packageVersion("ape")
    require(vegan); packageVersion("vegan")
    
    # Taxonomic filters
    ps1 <- readRDS("${ps}")
    
    # Export summary of filtered results
    seqateurs::summarise_taxa(ps1, "species", "sample_id") %>%
      tidyr::spread(key="sample_id", value="totalRA") %>%
      write.csv(file = "spp_sum_filtered.csv")

    seqateurs::summarise_taxa(ps1, "genus", "sample_id") %>%
      tidyr::spread(key="sample_id", value="totalRA") %>%
      write.csv(file = "gen_sum_filtered.csv")

    #Output fasta of all ASV's
    seqateurs::ps_to_fasta(ps1, "asvs_filtered.fasta", seqnames="species")

    #Output newick tree
    if(!is.null(phy_tree(ps1, errorIfNULL=FALSE))){
        ape::write.tree(phy_tree(ps1), file="tree_filtered.nwk")
    } else {
        cat("", file = "tree_filtered.nwk")
    }
    
    # output filtered phyloseq object
    saveRDS(ps1, "ps_filtered.rds") 
    """
}

/*
 *
 * Step 14: Taxon checks
 *
 */ 
 
if (params.check_ala) {
    process check_ala {
        tag { "check_ala" }
        publishDir "${params.outdir}/results/filtered", mode: "copy", overwrite: true

        input:
        file ps from tax_check_ala

        output:
        file "*.csv"

        script:
        """
        #!/usr/bin/env Rscript
        require(galah); packageVersion("galah")
        require(tidyverse); packageVersion("tidyverse")
        
        ps1 <- readRDS("${ps}")    
        
        # Check presence on ALA
        # First we need to set some data quality filters for ALA
        # To view available filters, run: find_field_values("basis_of_record")
        ala_quality_filter <- galah::select_filters(
              basisOfRecord = c("PreservedSpecimen", "LivingSpecimen",
                              "MaterialSample", "NomenclaturalChecklist"),
              profile = "ALA")

        spp_to_check <- ps1 %>%
          speedyseq::psmelt() %>%
          dplyr::group_by(family, genus, species) %>%
          dplyr::summarise(metabarcoding_reads = sum(Abundance)) %>%
          dplyr::filter(!str_detect(species, "__")) %>%
          dplyr::mutate(species = species %>% str_replace_all("_", " ")) 
          
        if(nrow(spp_to_check) > 0 ){
            ala_check <- spp_to_check %>%
              dplyr::mutate(
                species_present = purrr::map(species, function(x){
                # first check name
                query <- galah::select_taxa(x) %>% 
                  tibble::as_tibble()%>%
                  dplyr::filter(across(any_of("match_type"), ~!.x == "higherMatch"))
                # Then get occurance counts
                if(!is.null(query\$scientific_name)){
                  ala_occur <- galah::ala_counts(taxa=query, filters=ala_quality_filter)
                  return(data.frame(species_present = ifelse(ala_occur > 0, TRUE, FALSE), ALA_counts = ala_occur))
                } else {
                  return(data.frame(species_present = FALSE, ALA_counts = 0))
                }
                })) %>%
              tidyr::unnest(species_present) %>%
              dplyr::select(family, genus, species, species_present, ALA_counts, metabarcoding_reads)
            
            write_csv(ala_check, "ala_check.csv")    
        } else {
            # Write out empty csv so nextflow doesnt fail
            write_csv(as.data.frame(""), "ala_check.csv")
        }
        """
    }
}

if (params.check_afd) {
    process check_afd {
        tag { "check_afd" }
        publishDir "${params.outdir}/results/filtered", mode: "copy", overwrite: true

        input:
        file ps from tax_check_afd

        output:
        file "*.csv"

        script:
        """
        #!/usr/bin/env Rscript
        require(afdscraper); packageVersion("afdscraper")
        require(tidyverse); packageVersion("tidyverse")
        require(speedyseq); packageVersion("speedyseq")
        
        ps1 <- readRDS("${ps}")    
        
        # Check presence on AFD
        spp_to_check <- ps1 %>%
          speedyseq::psmelt() %>%
          dplyr::group_by(family, genus, species) %>%
          dplyr::summarise(metabarcoding_reads = sum(Abundance)) %>%
          dplyr::filter(!str_detect(species, "__")) %>%
          dplyr::mutate(species = species %>% str_replace_all("_", " ")) 
          
        if(nrow(spp_to_check) > 0 ){
              afd_check <- spp_to_check %>%
              dplyr::mutate(
                family_present = afdscraper::check_afd_presence(family),
                genus_present = afdscraper::check_afd_presence(genus),
                species_present = afdscraper::check_afd_presence(species)
              ) %>%
              dplyr::select(family, family_present, genus, genus_present, 
                            species, species_present, metabarcoding_reads)
                
            write_csv(afd_check, "afd_check.csv")    
        } else {
            # Write out empty csv so nextflow doesnt fail
            write_csv(as.data.frame(""), "afd_check.csv")
        }
        """
    }
}
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
