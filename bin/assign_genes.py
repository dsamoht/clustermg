#!/usr/bin/env python3
"""Hand every gene that matches a final-MAG representative over to it; set the rest aside.

`--hits` is the mmseqs2 search of the genes against the final-MAG representatives, already
filtered on the clustering thresholds. A gene with several hits joins the best one: highest
bitscore, then highest identity, then the smallest representative id, so that ties resolve
the same way on every run.
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="assign genes to their best final-MAG representative.")
    parser.add_argument("--hits", help="mmseqs2 hits: query, target, fident, alnlen, qcov, tcov, evalue, bits.")
    parser.add_argument("faa", nargs="*", help="proteins searched against the final-MAG representatives.")
    return parser.parse_args()


def better(hit, current):
    """hits are (bits, fident, target)"""
    if hit[0] != current[0]:
        return hit[0] > current[0]
    if hit[1] != current[1]:
        return hit[1] > current[1]
    return hit[2] < current[2]


def best_hits(path):
    best = {}
    with open(path) as handle:
        for line in handle:
            query, target, fident, _alnlen, _qcov, _tcov, _evalue, bits = line.rstrip("\n").split("\t")
            hit = (float(bits), float(fident), target)
            if query not in best or better(hit, best[query]):
                best[query] = hit
    return {query: hit[2] for query, hit in best.items()}


def main():
    args = parse_args()

    best = best_hits(args.hits) if args.hits else {}
    seen = set()
    with open("assigned.tsv", "w") as assigned, open("unassigned.faa", "w") as unassigned:
        keep = False
        for path in sorted(args.faa):
            with open(path) as faa:
                for line in faa:
                    if line.startswith(">"):
                        gene = line[1:].split()[0]
                        seen.add(gene)
                        if gene in best:
                            assigned.write(f"{best[gene]}\t{gene}\n")
                            keep = False
                        else:
                            keep = True
                    if keep:
                        unassigned.write(line)

    unknown = set(best) - seen
    if unknown:
        sys.exit(f"error: {len(unknown)} hit(s) for genes absent from the inputs, e.g. {sorted(unknown)[0]}")


if __name__ == "__main__":
    main()
