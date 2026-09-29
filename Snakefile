import glob

configfile: "config.yaml"

IN = config["input_dir"]
OUT = config["outdir"]
SID = config["sample_id"]
TARGET_GENUS = config["target_genus"]


def one(pattern):
    """The single input file matching pattern (errors if zero or many)."""
    hits = sorted(glob.glob(pattern))
    if len(hits) != 1:
        raise WorkflowError(f"expected exactly one file matching {pattern}, found {hits}")
    return hits[0]


def gpu_flag(_):
    if config["microscopy"]["gpu"]:
        return "--gpu"
    return "--no-gpu"


def wanted(_):
    targets = [OUT + "counter.csv"]
    if config["run_micro"]:
        targets.append(OUT + "microscopy.csv")
    if config["strain_composition"]["refdir"]:
        targets.append(OUT + "strain_composition.csv")
    if config["run_contamination"]:
        targets.append(OUT + SID + ".contamination.tsv")
    return targets


rule all:
    input:
        wanted,


rule counter:
    input:
        fcs=lambda wc: one(IN + "flow_cytometry/*.fcs"),
    output:
        OUT + "counter.csv",
    conda:
        "counter/environment.yml"
    params:
        c=config["counter"],
    shell:
        "python counter/fcs_count.py '{input.fcs}' "
        "--od {params.c[od]} --uL {params.c[uL]} "
        "--sample-id {SID} --out {output}"


rule microscopy:
    input:
        imgdir=IN + "microscopy_images",
    output:
        OUT + "microscopy.csv",
    conda:
        "microscopy/environment.yml"
    params:
        model=config["microscopy"]["model"],
        gpu=gpu_flag,
    shell:
        "python microscopy/microscopy_analysis.py --imgdir {input.imgdir} "
        "--model {params.model} --sample-id {SID} {params.gpu} --out {output}"


rule strain_composition:
    input:
        reads=lambda wc: one(IN + "reads/*.fastq.gz"),
    output:
        OUT + "strain_composition.csv",
    conda:
        "strain_composition/environment.yml"
    threads: workflow.cores
    params:
        refdir=config["strain_composition"]["refdir"],
    shell:
        "python strain_composition/em_abundance.py --reads {input.reads} "
        "--refdir {params.refdir} --sample-id {SID} --threads {threads} --out {output}"


rule kraken2:
    input:
        reads=lambda wc: one(IN + "reads/*.fastq.gz"),
    output:
        kraken=temp(OUT + SID + ".kraken"),
        report=OUT + SID + ".kreport",
    conda:
        "envs/kraken2.yaml"
    threads: workflow.cores
    params:
        db=config["contamination"]["kraken_db"],
    shell:
        "kraken2 --db {params.db} --threads {threads} --gzip-compressed "
        "--output {output.kraken} --report {output.report} {input.reads}"


rule contamination_report:
    input:
        report=OUT + SID + ".kreport",
    output:
        tsv=OUT + SID + ".contamination.tsv",
    run:
        rank_of = {"G": "genus", "S": "species"}
        rows = []
        with open(input.report) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                rank = cols[3]
                if rank not in rank_of:
                    continue
                name = cols[5].strip()
                if TARGET_GENUS.capitalize() in name:
                    group = TARGET_GENUS.capitalize()
                else:
                    group = "other"
                rows.append((rank_of[rank], name, group, cols[1], cols[0]))
        with open(output.tsv, "w") as out:
            out.write("level\ttaxon\tgroup\treads_classified\trelative_abundance_pct\n")
            for level, name, group, reads, pct in rows:
                out.write(f"{level}\t{name}\t{group}\t{reads}\t{pct}\n")
