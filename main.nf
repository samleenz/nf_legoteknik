#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// nextflow.preview.output = true

/////////////////////////////////////////
///////////// NF_LEGOTEKNIK /////////////
/////////////////////////////////////////
// Processing of LeGO-teknik scRNA-seq //
/////////////////////////////////////////

// Pipeline overview
    // 1. make the samplesheet for bcl-convert
    // 2. generate fastq files from the bcl directory
    //      2a. gene expression fastqs
    //      2b. teknkik barcode fastqs
    // 3. make the sample table for cellranger
    //      3a. channel of sample IDs
    //      3b. create the library csv files
    // 4. run cellranger per sample
    // 5. run multiqc on the bcl-convert and cellranger outputs



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
    tuple val(run_id), path(bcl_dir), path(bcl_samplesheet)

    output:
    tuple val(run_id), path("fastq/teknik")

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
    tuple val(run_id), path(bcl_dir), path(bcl_samplesheet)

    output:
    tuple val(run_id), path("fastq/gex")

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

// Process to merge fastqs from multiple runs per sample
// - creates per-sample directories:
//   <sample>/gex and <sample>/teknik
// - concatenates run-level fastqs into one file per read type:
//   <sample>_S1_[R1|R2|I1|I2]_001.fastq.gz
process merge_fastqs {
    cpus 1
    memory "4G"
    time "20m"

    input:
    tuple val(sample), val(gex_run_manifest), val(teknik_run_manifest)

    output:
    tuple val(sample), path(sample)

    script:
    """
    mkdir -p ${sample}/gex ${sample}/teknik

    cat > gex_runs.tsv <<'EOF'
${gex_run_manifest}
EOF

    cat > teknik_runs.tsv <<'EOF'
${teknik_run_manifest}
EOF

    merge_modality() {
      manifest="\$1"
      out_dir="\$2"

      found_any=0
      for read_type in R1 R2 I1 I2; do
        list_file="\${out_dir}_\${read_type}.list"
        : > "\$list_file"

        while IFS=\$'\\t' read -r run_id run_dir; do
          [ -z "\$run_dir" ] && continue
          for fq in "\$run_dir"/${sample}_*.fastq.gz; do
            [ -e "\$fq" ] || continue
            base=\$(basename "\$fq")
            case "\$base" in
              *_\${read_type}_*.fastq.gz)
                echo "\$fq" >> "\$list_file"
                ;;
            esac
          done
        done < "\$manifest"

        if [ -s "\$list_file" ]; then
          found_any=1
          sort -u "\$list_file" > "\${list_file}.sorted"
          out_fastq="\$out_dir/${sample}_S1_\${read_type}_001.fastq.gz"
          : > "\$out_fastq"
          while IFS= read -r fq; do
            cat "\$fq" >> "\$out_fastq"
          done < "\${list_file}.sorted"
        fi
      done

      # Shell convention: 0 = success (files found), 1 = no files found
      if [ "\$found_any" -eq 1 ]; then
        return 0
      else
        return 1
      fi
    }

    if ! merge_modality gex_runs.tsv "${sample}/gex"; then
      echo "No GEX fastqs found for sample ${sample} across runs" >&2
      exit 1
    fi
    if [ ! -s "${sample}/gex/${sample}_S1_R1_001.fastq.gz" ] || [ ! -s "${sample}/gex/${sample}_S1_R2_001.fastq.gz" ]; then
      echo "Merged GEX fastqs for ${sample} must include both R1 and R2" >&2
      exit 1
    fi

    # teknik is optional; ignore non-zero return when absent
    merge_modality teknik_runs.tsv "${sample}/teknik" || true
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
    tuple val(sample), path(merged_fastqs)
    path(txome)
    path(feature_reference)

    output:
    path(sample)

    script:
    """
    gex_r1_count=\$(find "${merged_fastqs}/gex" -maxdepth 1 -type f -name "${sample}_S1_R1_001.fastq.gz" | wc -l)
    gex_r2_count=\$(find "${merged_fastqs}/gex" -maxdepth 1 -type f -name "${sample}_S1_R2_001.fastq.gz" | wc -l)
    teknik_any_count=\$(find "${merged_fastqs}/teknik" -maxdepth 1 -type f -name "${sample}_S1_*.fastq.gz" | wc -l)

    if [ "\$gex_r1_count" -eq 0 ] || [ "\$gex_r2_count" -eq 0 ]; then
      echo "Missing required GEX R1/R2 fastqs for sample ${sample} in ${merged_fastqs}/gex" >&2
      exit 1
    fi

    cat > lib.csv <<EOF
fastqs,sample,library_type
${merged_fastqs}/gex,${sample},Gene Expression
EOF

    feature_ref_arg=""
    if [ "\$teknik_any_count" -gt 0 ]; then
      echo "${merged_fastqs}/teknik,${sample},Antibody Capture" >> lib.csv
      feature_ref_arg="--feature-ref=${feature_reference}"
    fi

    cellranger count \
        --id=${sample} \
        --libraries=lib.csv \
        --transcriptome=${txome} \
        \$feature_ref_arg \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()} \
        --create-bam=false \
        --nosecondary
    """    
}

// process run_multiqc {
//     cpus 1
//     memory "2.G"
//     time "30.m"
//     container params.container_multiqc
//     // publishDir "${params.out_dir}/multiqc", mode: 'copy' 
    

//     input:
//     tuple val(sample), path(cr_dir) 


//     output:
//     // path("${sample}/*")
//     tuple val(sample), path(out) 


//     script:
//     out="${sample}/multiqc*"
//     """
//     mkdir -p ${sample}
//     multiqc -f -o ${sample} $cr_dir
//     """
// }

process run_multiqc {
    cpus 1
    memory "2.G"
    time "30.m"
    container params.container_multiqc
    

    input:
    path(inputs) 
    val(name)


    output:
    path(output) 


    script:
    output="multiqc"
    title="${name} QC report"
    """
    mkdir -p ${output}
    multiqc -f -o ${output} -b "${title}" ${inputs}
    """
}

// entry workflow 
workflow {
    main:
    // 0. parse run manifest (or fallback to single run from params.bcl_dir)
    // expected run_manifest columns: run_id,bcl_dir,has_gex,has_teknik
    if( params.run_manifest ) {
        Channel.fromPath(params.run_manifest)
            .ifEmpty { error "run_manifest not found: ${params.run_manifest}" }
            .splitCsv(header: true)
            .map { row ->
                ['run_id','bcl_dir','has_gex','has_teknik'].each { required_col ->
                    if( !row.containsKey(required_col) ) {
                        error "run_manifest is missing required column: ${required_col}"
                    }
                }

                def normBool = { v ->
                    def s = v?.toString()?.trim()?.toLowerCase()
                    if( s in ['true','t','1','yes','y'] ) return true
                    if( s in ['false','f','0','no','n'] ) return false
                    error "Invalid boolean '${v}' in run_manifest. Use true/false."
                }

                tuple(
                    row.run_id.toString(),
                    file(row.bcl_dir.toString()),
                    normBool(row.has_gex),
                    normBool(row.has_teknik)
                )
            }
            .set { run_manifest_ch }
    } else {
        if( !params.bcl_dir ) {
            error "Set either params.run_manifest or params.bcl_dir"
        }

        Channel.of(
            tuple(
                file(params.bcl_dir).baseName,
                file(params.bcl_dir),
                true,
                true
            )
        ).set { run_manifest_ch }
    }

    // run_manifest_ch.view()

    // 1. make the samplesheet for bcl-convert
    bcl_samplesheets_ch = bclconvert_samplesheets(
        file(params.sampletable),
        file(params.single_template),
        file(params.dual_template),
        file(params.single_index),
        file(params.dual_index)
    )

    // bcl_samplesheets_ch.gex_samplesheet.view()
    // bcl_samplesheets_ch.legotek_samplesheet.view()


    // 2. generate fastq files from bcl directories listed in run manifest

    // 2a. gene expression fastqs
    run_gex_ch = run_manifest_ch
        .filter { run_id, bcl_dir, has_gex, has_teknik -> has_gex }
        .map { run_id, bcl_dir, has_gex, has_teknik -> tuple(run_id, bcl_dir) }

    demux_gex_input_ch = run_gex_ch
        .combine(bcl_samplesheets_ch.gex_samplesheet)
        .map { row -> tuple(row[0], row[1], row[2]) }

    fastq_gex_by_run_ch = demux_gex(demux_gex_input_ch)
    fastq_gex_ch = fastq_gex_by_run_ch.map { run_id, p -> p }

    // 2b. teknkik barcode fastqs
    run_legotek_ch = run_manifest_ch
        .filter { run_id, bcl_dir, has_gex, has_teknik -> has_teknik }
        .map { run_id, bcl_dir, has_gex, has_teknik -> tuple(run_id, bcl_dir) }

    demux_legotek_input_ch = run_legotek_ch
        .combine(bcl_samplesheets_ch.legotek_samplesheet)
        .map { row -> tuple(row[0], row[1], row[2]) }

    fastq_legotek_by_run_ch = demux_legotek(demux_legotek_input_ch)
    fastq_legotek_ch = fastq_legotek_by_run_ch.map { run_id, p -> p }
    
    // fastq_gex_ch.view()
    // fastq_legotek_ch.view()

    // 3. make channels for per-sample merge + cellranger

    // 3a. channel of sample IDs
    Channel.fromPath( params.sampletable )
        .splitCsv( header: true)
        .map { it.sample_id }
        .set { sampletable_ch }

    // sampletable_ch.view()

    // 3b. collect run-level fastq outputs so each sample can be merged across runs
    gex_run_manifest_ch = fastq_gex_by_run_ch
        .collect(flat: false)
        .map { runs ->
            runs.collect { run_id, run_dir -> "${run_id}\t${run_dir}" }.join('\n')
        }

    teknik_run_manifest_ch = fastq_legotek_by_run_ch
        .collect(flat: false)
        .map { runs ->
            runs.collect { run_id, run_dir -> "${run_id}\t${run_dir}" }.join('\n')
        }

    merge_fastqs_input_ch = sampletable_ch
        .combine(gex_run_manifest_ch)
        .combine(teknik_run_manifest_ch)
        .map { row -> tuple(row[0], row[1], row[2]) }

    merged_fastqs_ch = merge_fastqs(merge_fastqs_input_ch)

    // 4. run cellranger per sample
    cellranger_ch = cellranger_count(
        merged_fastqs_ch,
        file(params.transcriptome),
        file(params.feature_reference)
    )

    // cellranger_ch.view()


    // 5. run multiqc on the bcl-convert and cellranger outputs
    // keep only the 2nd element from each tuple
    // bcl_gex_path = fastq_gex_ch.map { sample, path -> path }
    // bcl_legotek_path = fastq_legotek_ch.map { sample, path -> path }
    // cr_path = cellranger_ch.map { sample, path -> path }

    //     cr_path
    //         .concat(bcl_gex_path)
    //         .concat(bcl_legotek_path)
    //         .collect()
    //         .view()

    multiqc_ch = run_multiqc(
        cellranger_ch.concat(fastq_gex_ch).concat(fastq_legotek_ch).collect(),
        file(params.bcl_dir).baseName
    )

    publish:
    cellranger_ch = cellranger_ch
    multiqc_ch    = multiqc_ch
}

output {
  cellranger_ch {
    path { p -> "cellranger" }   // p is the <sample>/ directory
  }
  multiqc_ch {
    path { p -> "multiqc" }                    // p is the multiqc directory
  }
}

// output {
//     cellranger_ch {
//     path { tup -> 
//       def (sample, outs) = tup
//       "cellranger/${sample}/"
//     }
//   }
//     multiqc_ch {
//     path { tup ->
//       def (sample, report) = tup
//       "multiqc/${sample}/"
//       }
//     }

// }
