#!/usr/bin/env nextflow
nextflow.enable.dsl=2
nextflow.preview.output = true

// --- params banner (runs at script load) ---
def _mask = { value, key ->
    def k = key?.toLowerCase() ?: ''
    ['password','pass','token','secret','apikey','api_key','key'].any { k.contains(it) } ? '********' : value
}
def _pp = { v ->
    if( v instanceof java.nio.file.Path ) return v.toString()
    if( v instanceof Collection )         return v.join(', ')
    return v
}

log.info "-- Launching: ${workflow.manifest?.name ?: 'pipeline'} v${workflow.manifest?.version ?: 'dev'}"
log.info "-- Run name : ${workflow.runName}"
log.info "-- Workdir  : ${workflow.workDir}"
log.info "----- PARAMS -----"
params.keySet().sort().each { k ->
    log.info String.format("  %-22s = %s", k, _mask(_pp(params[k]), k))
}
log.info "-------------------"

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

// Process to create the library csv file for cellranger
process cellranger_library {
    container params.container_r
    cpus 1
    memory "4G"
    time "10m"

    input:
    val(sample)
    path(library_template)

    output:
    tuple val(sample), path("lib.csv")



    script:
    """
    sed "s/{{SAMPLE}}/${sample}/g" ${library_template} > lib.csv
    """ 
}

// Process to run cellranger on the gex + LeGO fastq files
// - note, could think about doing a light/heavy version of this 
//   where --create-bam and --nosecondary are toggled.
process cellranger_count {
    container params.container_cellranger 
    cpus 32
    memory "128.G"
    time {"1.0h"}
    // publishDir "${params.out_dir}/cellranger", mode: 'copy'

    input:
    tuple val(sample), path(library_csv)
    path(fastq_gex)
    path(fastq_teknik)
    path(txome)
    path(feature_reference)

    output:
    tuple val(sample), path(out)

    script:
    out=sample
    """
    cellranger count \
        --id=${out} \
        --libraries=${library_csv} \
        --transcriptome=${txome} \
        --feature-ref=${feature_reference} \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()} \
        --create-bam=false \
        --nosecondary
    """    
}

process run_multiqc {
    cpus 1
    memory "2.G"
    time "30.m"
    container params.container_multiqc
    // publishDir "${params.out_dir}/multiqc", mode: 'copy' 
    

    input:
    tuple val(sample), path(cr_dir) 


    output:
    // path("${sample}/*")
    tuple val(sample), path(out) 


    script:
    out="${sample}/multiqc*"
    """
    mkdir -p ${sample}
    multiqc -f -o ${sample} $cr_dir
    """
}

// Define the workflow by chaining the processes
workflow {
    main:
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
    
    // fastq_gex_ch.view()
    // fastq_legotek_ch.view()

    // channel of sample names (from sampletable)
    Channel.fromPath( params.sampletable )
        .splitCsv( header: true)
        .map { it.sample_id }
        .set { sampletable_ch }

    // sampletable_ch.view()

    // create the library csv files
    library_csv_ch = cellranger_library(
        sampletable_ch,
        file(params.library_template)
    )

    // library_csv_ch.view()

    cellranger_ch = cellranger_count(
        library_csv_ch,
        fastq_gex_ch,
        fastq_legotek_ch,
        file(params.transcriptome),
        file(params.feature_reference)
    )

    // cellranger_ch.view()

    multiqc_ch = run_multiqc(
        cellranger_ch
    )

    publish:
    // cellranger_ch = cellranger_ch
    // multiqc_ch = multiqc_ch
    cellranger_ch >> "cellranger"
    multiqc_ch >> "multiqc"
}

output {
  directory params.outputDir
  mode 'copy'
}

// output {
// //     cellranger_ch {
// //     path { tup -> 
// //       def (sample, outs) = tup
// //       "cellranger/${sample}/"
// //     }
// //   }
// //     multiqc_ch {
// //     path { tup ->
// //       def (sample, report) = tup
// //       "multiqc/${sample}/"
// //       }
// //     }

// }
