"""Resolve read assignments across high-ANI strains by expectation-maximization.

At ~99% ANI a long read maps near-equally to several genomes, so best-hit
assignment is arbitrary. The signal that separates strains is the *difference*
in mismatches at strain-specific variant sites. We give each (read, strain)
pair a likelihood weight = penalty ** mismatches, then run EM:

    E-step: responsibility r[read,s] proportional to theta[s] * weight[read,s]
    M-step: theta[s] = sum of responsibilities for s, normalised

Random ONT errors are shared across a read's alignments, so they cancel in the
E-step normalisation; only the strain-discriminating mismatches drive the split.
theta is the fraction of *reads* per strain. We then soft-assign each read's
*bases* to strains by those EM responsibilities, divide by genome size and
renormalise for relative abundance (cells map ~ bases / genome_size).
"""

import argparse
import csv
import glob
import io
import json
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

import pandas as pd

PENALTY = 0.25
VALID_BASES = set("ATGCRYSWKMBDHVNX")


def default_threads():
    return os.cpu_count()


def tag(strain, contig):
    """Strain-qualified contig id, so identical contig names across refs stay distinct."""
    return f"{'_'.join(strain.split())}_{contig}"


def parse_refs(ref_paths):
    """{contig: strain} and {strain: genome_size}, straight from the reference FASTAs.

    Raises ValueError on anything that isn't nucleotide FASTA.
    """
    contig2strain = {}
    genome_size = {}
    for strain, path in ref_paths.items():
        size = 0
        with open(path, errors="replace") as fh:
            for lineno, line in enumerate(fh, 1):
                seq = line.strip()
                if lineno == 1 and not seq.startswith(">"):
                    raise ValueError(
                        f"{strain}: not a FASTA file — the first line must start with '>'"
                    )
                if seq.startswith(">"):
                    contig = seq[1:].split()
                    if not contig:
                        raise ValueError(f"{strain}: line {lineno} has an empty header")
                    contig2strain[tag(strain, contig[0])] = strain
                    continue
                bad = set(seq.upper()) - VALID_BASES
                if bad:
                    raise ValueError(
                        f"{strain}: line {lineno} has invalid sequence "
                        f"characters: {''.join(sorted(bad))}"
                    )
                size += len(seq)
        if not size:
            raise ValueError(f"{strain}: FASTA contains no sequence")
        genome_size[strain] = size
    return contig2strain, genome_size


def write_combined_ref(ref_paths, combined):
    """Concatenate the per-strain refs into one target, rewriting headers via tag().

    A read's near-equal hits to every high-ANI strain only survive minimap2's
    secondary output when all strains share a single target.
    """
    with open(combined, "w") as out:
        for strain, path in ref_paths.items():
            with open(path, errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        out.write(f">{tag(strain, line[1:].lstrip())}")
                    else:
                        out.write(line)


def run_minimap2(combined, reads_path, paf_path, threads):
    """Align reads to the combined ref, keeping secondary hits for the EM pass."""
    with open(paf_path, "w") as pafout:
        subprocess.run(
            [
                "minimap2",
                "-t",
                str(threads),
                "-x",
                "map-ont",
                "-N",
                "10",
                "-p",
                "0.8",
                combined,
                reads_path,
            ],
            stdout=pafout,
            check=True,
        )


def strain_abundance(reads_path, ref_paths, threads, penalty=PENALTY):
    """Full pipeline: refs + reads -> {strain: relative abundance}. Raises ValueError on bad FASTA."""
    contig2strain, genome_size = parse_refs(ref_paths)
    with tempfile.TemporaryDirectory() as d:
        combined = os.path.join(d, "combined_ref.fasta")
        write_combined_ref(ref_paths, combined)
        paf = os.path.join(d, "aligned.paf")
        run_minimap2(combined, reads_path, paf, threads)
        read_hits, read_block = load_paf(paf, contig2strain)
    return abundance(read_hits, read_block, genome_size, penalty)


def run_em(read_hits, strains, penalty, max_iter=200, tol=1e-7):
    """read_hits: {read_id: {strain: fewest_mismatches}}. Returns {strain: read_fraction}."""
    # Per-read, per-strain likelihood weight.
    weights = {}
    for read, hits in read_hits.items():
        w = {}
        for strain, mm in hits.items():
            w[strain] = penalty**mm
        weights[read] = w

    # Seed theta from a hard best-hit assignment (+1 pseudocount avoids zero-lock).
    theta = {}
    for strain in strains:
        theta[strain] = 1.0
    for read, hits in read_hits.items():
        best = min(hits, key=hits.get)
        theta[best] += 1.0
    total = sum(theta.values())
    for strain in strains:
        theta[strain] /= total

    for _ in range(max_iter):
        counts = {}
        for strain in strains:
            counts[strain] = 0.0
        for read, w in weights.items():
            denom = 0.0
            for strain, wt in w.items():
                denom += theta[strain] * wt
            if denom == 0.0:
                continue
            for strain, wt in w.items():
                counts[strain] += theta[strain] * wt / denom
        n = sum(counts.values())
        new = {}
        for strain in strains:
            if n > 0:
                new[strain] = counts[strain] / n
            else:
                new[strain] = theta[strain]
        delta = 0.0
        for strain in strains:
            delta += abs(new[strain] - theta[strain])
        theta = new
        if delta < tol:
            break
    return theta


def load_paf(path, contig2strain):
    """Best (fewest-mismatch) alignment of each read to each strain, plus its block length."""
    # PAF 0-indexed: 0=read, 5=target contig, 9=residue matches, 10=block length
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 5, 9, 10],
        names=["ReadID", "Contig", "Matches", "Block"],
    )
    df["Strain"] = df["Contig"].map(contig2strain)
    df = df.dropna(subset=["Strain"])
    df["Mismatches"] = df["Block"] - df["Matches"]

    read_hits = defaultdict(dict)
    read_block = defaultdict(
        dict
    )  # aligned bases of the chosen alignment per (read, strain)
    for row in df.itertuples(index=False):
        cur = read_hits[row.ReadID].get(row.Strain)
        if cur is None or row.Mismatches < cur:
            read_hits[row.ReadID][row.Strain] = row.Mismatches
            read_block[row.ReadID][row.Strain] = row.Block
    return read_hits, read_block


def base_totals(read_hits, read_block, theta, penalty):
    """Soft-assign each read's aligned bases to strains by final EM responsibility."""
    bases = defaultdict(float)
    for read, hits in read_hits.items():
        w = {}
        denom = 0.0
        for strain, mm in hits.items():
            wt = theta[strain] * penalty**mm
            w[strain] = wt
            denom += wt
        if denom == 0.0:
            continue
        for strain, wt in w.items():
            bases[strain] += read_block[read][strain] * wt / denom
    return bases


def abundance(read_hits, read_block, genome_size, penalty):
    strains = sorted(genome_size)
    theta = run_em(read_hits, strains, penalty)
    bases = base_totals(read_hits, read_block, theta, penalty)
    cov = {}
    for strain in strains:
        cov[strain] = bases[strain] / genome_size[strain]
    tot = sum(cov.values())
    out = {}
    for strain in strains:
        if tot > 0:
            out[strain] = cov[strain] / tot
        else:
            out[strain] = 0.0
    return out


def write_csv(rows, fieldnames, out):
    """Write rows (list of dicts) as CSV to stdout, and to `out` too if given."""
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)
    text = buf.getvalue()
    if out:
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        with open(out, "w", newline="") as f:
            f.write(text)
    print(text, end="")


def find_refs(refdir):
    """{strain: path} for every <strain>.fna in refdir; strain name is the filename stem."""
    ref_paths = {}
    for path in sorted(glob.glob(os.path.join(refdir, "*.fna"))):
        strain = os.path.splitext(os.path.basename(path))[0]
        ref_paths[strain] = path
    return ref_paths


def main():
    p = argparse.ArgumentParser(
        description="Strain composition from Nanopore reads by EM over minimap2 hits."
    )
    p.add_argument("--reads", help="Nanopore raw reads (FASTQ)")
    p.add_argument(
        "--refdir",
        help="Directory of <strain>.fna reference FASTAs (filename stem = strain name)",
    )
    p.add_argument(
        "--sample-id", dest="sample_id", default="sample", help="Id for this sample"
    )
    p.add_argument("--out", help="Write result CSV here")
    p.add_argument("--threads", type=int, default=default_threads())
    p.add_argument("--penalty", type=float, default=PENALTY)
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args()

    if args.selftest:
        _selfcheck()
        return

    if not args.reads or not args.refdir:
        p.error("--reads and --refdir are required")
    if not os.path.isfile(args.reads):
        p.error(f"reads file not found: {args.reads}")
    ref_paths = find_refs(args.refdir)
    if not ref_paths:
        p.error(f"no <strain>.fna reference files in {args.refdir}")

    try:
        abund = strain_abundance(args.reads, ref_paths, args.threads, args.penalty)
    except ValueError as e:
        print(json.dumps({"error": str(e)}), file=sys.stderr)
        sys.exit(1)

    rows = []
    for strain in sorted(abund):
        rows.append(
            {"sample_id": args.sample_id, "strain": strain, "abundance": abund[strain]}
        )
    write_csv(rows, ["sample_id", "strain", "abundance"], args.out)


def _selfcheck():
    # 70 reads truly from A, 20 from B, 10 from C; true strain has 0 mismatches,
    # the others 5. EM should recover ~0.7 / 0.2 / 0.1 (equal genome sizes).
    read_hits = {}
    i = 0
    for strain, n in (("A", 70), ("B", 20), ("C", 10)):
        for _ in range(n):
            hits = {"A": 5, "B": 5, "C": 5}
            hits[strain] = 0
            read_hits[i] = hits
            i += 1
    theta = run_em(read_hits, ["A", "B", "C"], penalty=0.25)
    assert abs(theta["A"] - 0.7) < 0.02, theta
    assert abs(theta["B"] - 0.2) < 0.02, theta
    assert abs(theta["C"] - 0.1) < 0.02, theta
    print("ok", theta)

    # Equal block lengths + equal genomes → base-weighted abundance == read fractions.
    read_block = {}
    for r in read_hits:
        read_block[r] = {}
        for s in read_hits[r]:
            read_block[r][s] = 1000
    genome_size = {"A": 1.0, "B": 1.0, "C": 1.0}
    abund = abundance(read_hits, read_block, genome_size, penalty=0.25)
    assert abs(abund["A"] - 0.7) < 0.02, abund
    assert abs(abund["B"] - 0.2) < 0.02, abund
    assert abs(abund["C"] - 0.1) < 0.02, abund

    # Double the aligned blocks of C's reads → C's base share grows.
    for r in read_hits:
        if read_hits[r]["C"] == 0:
            read_block[r]["C"] = 2000
    skewed = abundance(read_hits, read_block, genome_size, penalty=0.25)
    assert skewed["C"] > abund["C"], skewed
    print("ok bases", abund, skewed)

    # FASTA validation (moved here from strain_composition/server.py).
    def _write(d, text):
        p = os.path.join(d, "t.fna")
        with open(p, "w") as fh:
            fh.write(text)
        return {"s1": p}

    def _fails(d, text):
        try:
            parse_refs(_write(d, text))
        except ValueError:
            return True
        return False

    with tempfile.TemporaryDirectory() as d:
        c2s, sizes = parse_refs(_write(d, ">c1 desc\nACGTN\n\nrykm\n"))
        assert c2s == {"s1_c1": "s1"}, c2s
        assert sizes == {"s1": 9}, sizes
        # Same contig name in two refs must not collide.
        assert tag("A", "contig_1") != tag("B", "contig_1")
        assert tag("Strain A", "c1") == "Strain_A_c1"
        assert _fails(d, "ACGT\n"), "missing header not caught"
        assert _fails(d, ">c1\nACGTZ\n"), "bad base not caught"
        assert _fails(d, ">c1\n"), "empty sequence not caught"
        assert _fails(d, ">\nACGT\n"), "empty header not caught"
    print("ok fasta")


if __name__ == "__main__":
    main()
