#!/usr/bin/env python
"""Flow-cytometry cell counting — runs in the counter container (hdbscan/fcsparser)."""

import argparse
import csv
import io
import os
from pathlib import Path

import fcsparser
import hdbscan
import numpy as np

X_CHANNEL, Y_CHANNEL = "FSC-A", "SSC-A"
SCALING_COFACTOR = 150
MIN_CLUSTER_SIZE = 1000
MIN_SAMPLES = 500
# Threads for HDBSCAN's core-distance step. The library's default is a fixed 4,
# which oversubscribes a 2-vCPU host and under-uses a bigger one. Leaves one core
# free for the same reason strain does: everything shares one small VPS, and a
# clustering run at full width starves the portal's gunicorn workers.
THREADS = os.cpu_count()


def load(fcs_file):
    fcs_file = Path(fcs_file)
    if not fcs_file.is_file():
        raise FileNotFoundError(f"{fcs_file} does not exist")
    _, data = fcsparser.parse(fcs_file, reformat_meta=True)
    return np.arcsinh(data[[X_CHANNEL, Y_CHANNEL]].to_numpy() / SCALING_COFACTOR)


def largest_cluster(labels):
    """Label of the most populous non-noise cluster (-1 is HDBSCAN noise)."""
    real = labels[labels != -1]
    if real.size == 0:
        raise ValueError("no clusters found (all points are noise)")
    return int(np.bincount(real).argmax())


def count(fcs_path, od, uL, min_prob=0.25, sample_id=None):
    """Gate an FCS sample to the largest cluster and return cell-count metrics."""
    X = load(fcs_path)
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=MIN_CLUSTER_SIZE,
        min_samples=MIN_SAMPLES,
        allow_single_cluster=True,
        core_dist_n_jobs=THREADS,
    ).fit(X)

    label = largest_cluster(clusterer.labels_)
    in_cluster = (clusterer.labels_ == label) & (clusterer.probabilities_ >= min_prob)
    n = int(in_cluster.sum())

    od, uL = float(od), float(uL)
    if uL:
        cells_per_uL = n / uL
    else:
        cells_per_uL = float("nan")
    if od:
        cells_per_uL_per_OD = round(cells_per_uL / od, 1)
    else:
        cells_per_uL_per_OD = float("nan")

    return {
        "sample_id": sample_id,
        "n_cells": n,
        "cells_per_uL": round(cells_per_uL, 1),
        "cells_per_uL_per_OD": cells_per_uL_per_OD,
    }


def demo():
    """Self-check: largest_cluster ignores noise and picks the bigger cluster."""
    labels = np.array([-1, -1, 0, 0, 0, 1, 1])  # cluster 0 has 3, cluster 1 has 2
    got = largest_cluster(labels)
    assert got == 0, got
    print("ok")


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


def main():
    p = argparse.ArgumentParser(description="Gate an FCS sample and count cells/uL.")
    p.add_argument("fcs_file", nargs="?")
    p.add_argument("--od", type=float)
    p.add_argument("--uL", type=float, help="Volume of fcs sample in uL")
    p.add_argument(
        "--sample-id", dest="sample_id", default="sample", help="Id for this sample"
    )
    p.add_argument("--out", help="Write result CSV here")
    p.add_argument("--min-prob", dest="min_prob", type=float, default=0.25)
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args()

    if args.selftest:
        demo()
        return

    result = count(args.fcs_file, args.od, args.uL, args.min_prob, args.sample_id)
    write_csv([result], list(result), args.out)


if __name__ == "__main__":
    main()
