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

    output:
    path("gex_samplesheet.csv"), emit: gex_samplesheet


    script:
    """
    echo \\${PATH}

    Rscript make_bclconvert_samplesheets.R ${sampletable}   
    """
}


// Process for demultiplexing the LeGOteknik fastqs
process demux_teknik {
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
    bclconvert \
        --bcl-input-directory ${bcl_dir} \
        --output-directory "fastq/teknik" \
        --samplesheet ${bcl_samplesheet} \
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
    bclconvert \
        --bcl-input-directory ${bcl_dir} \
        --output-directory "fastq/gex" \
        --samplesheet ${bcl_samplesheet} \
        --bcl-num-conversion-threads 16 \
        --bcl-num-compression-threads 8
    """
}

// Define the workflow by chaining the processes
workflow {
    bcl_samplesheets_ch = bclconvert_samplesheets(
        file(params.sampletable)
    ).view ()
    
}