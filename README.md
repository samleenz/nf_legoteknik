# nf_legoteknik

A nextflow DSL2 pipeline to process 10X scrnaseq using LeGO-teknik barcodes

> **Note:** This pipeline runs on nextflow 24.02.2 (using the preview workflow output DSL)

## Prerequisites

- [Nextflow](https://www.nextflow.io/) installed.
- [Apptainer](https://apptainer.org/) installed.

## Repository Structure

```plaintext
.
```

## Getting Started

1. **Clone the Repository**  
   Ensure you have access to the private repo:
   ```bash
   git clone git@github.com:samleenz/nf_legoteknik.git
   cd nf_legoteknik
   ```

2. **Prepare the files**  
   
> *_sampletable.csv

```plaintext
sample_id,gex_index,feature_index
S001,SI-TT-A8,D702
S002,SI-TT-A9,D703
```

> `{{NAME}}_config.config`

```plaintext
params {
    outputDir      = '/data/Analysis/Lab_TPH/scrnaseq_legotek/{{NAME}}'
    bcl_dir      = '/data/Analysis/Lab_TPH/datasets_raw/{{NAME}}'
    sampletable  = '{{NAME}}_sampletable.csv'
    transcriptome = '/data/Analysis/Lab_TPH/references/10X/refdata-gex-GRCh38_and_GRCm39-2024-A'
}
```

3. **Run the Pipeline**  
   Launch the pipeline with:
   ```bash
   nextflow run main.nf -resume -with-report {{NAME}}_report.html -c {{NAME}}_config.config
   ```





