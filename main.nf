#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// *
// * It uses an Apptainer container specified by the 'params.container' parameter.
// */

// Process to create the bcl_dir samplesheets
process bclconvert_samplesheets {
    container params.container_r
    cpus 1
    memory "4G"
    time "10m"

    input:
    path(sampletable)
    path(single_template)
    path(dual_template)
    path(single_index)
    path(dual_index)

    output:
    path("gex_samplesheet.csv"), emit: gex_samplesheet
    path("legotek_samplesheet.csv"), emit: legotek_samplesheet


    script:
    """
    make_bclconvert_samplesheets.R ${sampletable} \
        ${single_template} ${dual_template} \
        ${single_index} ${dual_index} 
    """
}


// Process for demultiplexing the LeGOteknik fastqs
process demux_legotek {
    container params.container_bclconvert 
    cpus 24
    memory "64.G"
    time {"30m"}

    input:
    path(bcl_dir)
    path(bcl_samplesheet)

    output:
    path("fastq/teknik")

    script:
    """
    mkdir -p fastq

    bcl-convert \
        --bcl-input-directory ${bcl_dir} \
        --output-directory "fastq/teknik" \
        --sample-sheet ${bcl_samplesheet} \
        --bcl-num-conversion-threads 16 \
        --bcl-num-compression-threads 8
    """
}

// Process for demultiplexing the gex fastqs
process demux_gex {
    container params.container_bclconvert 
    cpus 24
    memory "64.G"
    time {"30m"}

    input:
    path(bcl_dir)
    path(bcl_samplesheet)

    output:
    path("fastq/gex")

    script:
    """
    mkdir -p fastq

    bcl-convert \
        --bcl-input-directory ${bcl_dir} \
        --output-directory "fastq/gex" \
        --sample-sheet ${bcl_samplesheet} \
        --bcl-num-conversion-threads 16 \
        --bcl-num-compression-threads 8
    """
}

// Define the workflow by chaining the processes
workflow {
    bcl_samplesheets_ch = bclconvert_samplesheets(
        file(params.sampletable),
        file(params.single_template),
        file(params.dual_template),
        file(params.single_index),
        file(params.dual_index)
    )

    // bcl_samplesheets_ch.gex_samplesheet.view()
    // bcl_samplesheets_ch.legotek_samplesheet.view()


    // generate fastq files from the bcl directory
    fastq_gex_ch = demux_gex(
        file(params.bcl_dir),
        bcl_samplesheets_ch.gex_samplesheet
    )

    fastq_legotek_ch = demux_legotek(
        file(params.bcl_dir),
        bcl_samplesheets_ch.legotek_samplesheet
    )
    
    fastq_gex_ch.view()
    fastq_legotek_ch.view()

    // 
}