# nf_legoteknik

A nextflow DSL2 pipeline to process 10X scrnaseq using LeGO-teknik barcodes

> **Note:** This pipeline runs on nextflow 25.10+ (workflow output DSL)

## Prerequisites

- [Nextflow](https://www.nextflow.io/) installed.
- [Apptainer](https://apptainer.org/) installed.
- `rsync` available on the server.
- The lab fileserver (CIFS) mounted and readable from the server.

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
   
> `*_sampletable.csv`

```plaintext
sample_id,gex_index,feature_index
S001,SI-TT-A8,D702
S002,SI-TT-A9,D703
```

> `{{NAME}}_config.config`

```plaintext
params {
    outputDir      = '/data/Analysis/Lab_TPH/scrnaseq_legotek/{{NAME}}'
    run_dir      = '/path/to/fileserver/sequencing runs/{{RUN_DIR}}'
    sampletable  = '{{NAME}}_sampletable.csv'
    feature_index_type = 'single'
    gex_index_type     = 'dual'
    transcriptome = '/data/Analysis/Lab_TPH/references/10X/refdata-gex-GRCh38_and_GRCm39-2024-A'
}
```

`feature_index_type` and `gex_index_type` can be set to either `single` or `dual`.

`run_dir` is the sequencing run directory **on the fileserver**
(e.g. `241009_A01234_0123_AHXXXXXX`). The pipeline copies it into the work
directory itself (`STAGE_BCL`, skipping `Thumbnail_Images/` and `Logs/`), so no
manual copy is needed. For multiple runs, use a
`run_manifest` CSV with columns `run_id,run_dir,has_gex,has_teknik`
(see `data/example_run_manifest.csv`).

3. **Authenticate to the fileserver**  
   CIFS credentials live in the login session keyring, so Nextflow must be
   launched from a shell (or tmux/screen session) that has already run
   `cifscreds`. Authenticate *inside* the tmux/screen session you will launch from:
   ```bash
   tmux new -s {{NAME}}
   cifscreds add <fileserver-host>
   ```
   If the mount isn't readable, the run fails at launch with a
   "file does not exist" error for `run_dir`.

4. **Run the Pipeline**  
   Launch the pipeline with:
   ```bash
   nextflow run main.nf -with-report {{NAME}}_report.html -c {{NAME}}_config.config
   ```
   `run_dir` can also be given on the command line. Parent directories on the
   fileserver contain spaces, so **always quote it**:
   ```bash
   nextflow run main.nf -c {{NAME}}_config.config \
       --run_dir '/path/to/fileserver/sequencing runs/241009_A01234_0123_AHXXXXXX'
   ```

### Cleanup and `-resume`

`cleanup = true` is set in `nextflow.config`, so after a **successful** run
the work directory (including the staged BCL copy and all intermediates) is
deleted. Published outputs are copied (not symlinked) to the output directory
and are kept.

- After a successful run, `-resume` will not reuse anything: the work dir is gone.
- A **failed** run keeps its work dir and can be resumed with `-resume` as usual.
- When a resumed run succeeds, cleanup only removes tasks that ran in that
  session. Tasks reused from the failed attempt (including the staged BCLs)
  are left behind. The pipeline prints a warning at the end of a resumed run
  with the exact command to remove them; once you've checked the results, run
  it from the launch directory. It cleans every run in that session (the
  failed attempts and the successful resume) and nothing else:
  ```bash
  nextflow log | awk -F'\t' '$6 ~ /<session_id>/ {print $3}' | xargs -n1 nextflow clean -f
  ```
