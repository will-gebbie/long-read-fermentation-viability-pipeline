# Viability CLI — run-all workflow

One Snakemake workflow that runs four analyses on a sample and writes CSV/TSV to `results/`:

| Tool                     | Input (from `input_data/`)              | Output                                           |
| ------------------------ | --------------------------------------- | ------------------------------------------------ |
| counter (flow cytometry) | `flow_cytometry/*.fcs`                  | `results/<sample>` → `counter.csv`               |
| microscopy (viability)   | `microscopy_images/`                    | `microscopy.csv`                                 |
| strain_composition       | `reads/*.fastq.gz` + `references/*.fna` | `strain_composition.csv`                         |
| contamination (kraken2)  | `reads/*.fastq.gz`                      | `<sample>.kreport`, `<sample>.contamination.tsv` |

## One-time setup

1. Install conda/mamba, then a snakemake env:
   ```
   conda create -n snakemake -c conda-forge -c bioconda snakemake pandas
   conda activate snakemake
   ```
2. Per-machine assets (too large to ship — set the paths in `config.yaml`):
   - **kraken2 DB** → `contamination.kraken_db`. Get one from https://benlangmead.github.io/aws-indexes/k2 (e.g. PlusPF-8 or PlusPF-16 depending on available RAM on machine).
   - **cellpose model** → `microscopy.model`

## Input layout

```
input_data/
  reads/            <sample>.fastq.gz          (strain_composition + contamination) — exactly one file
  flow_cytometry/   <sample>.fcs               (counter)                            — exactly one file
  microscopy_images/ *.tif                      (microscopy)                          — one or more images
  references/       *.fna  (one per strain)    (strain_composition)                 — one file per strain
```

Required per tool (a tool is skipped only if you turn it off in `config.yaml`):

- **counter** — one `.fcs` file in `flow_cytometry/`. More than one (or none) is an error.
- **microscopy** — a `microscopy_images/` directory holding at least one `.tif`. Skip by setting `run_micro: false`.
- **strain_composition** — one `.fastq.gz` in `reads/`, plus one `.fna` per strain in `references/`.
- **contamination** — the same one `.fastq.gz` in `reads/`. Skip by setting `run_contamination: false`.

Set `sample_id` and the tool knobs (`od`, `uL`, `gpu`, …) in `config.yaml`.

## Run

```
conda activate snakemake
snakemake --use-conda --conda-prefix ~/.snakemake-envs --cores 20
```

- `--use-conda` builds each tool's env from its `environment.yml` on first run, then reuses it.
- `--conda-prefix ~/.snakemake-envs` builds those envs **once** and shares them across runs/checkouts
- Add `-n` for a dry run (shows the plan without executing).

Results land in `results/`.
