#!/usr/bin/env python3
"""3-way intersection + binding-site refinement -> predicted ceRNA interactions.

Faithful reimplementation of the legacy ``code/step1..step5`` logic:

1. Take the miRNA-target pairs predicted by all three tools (miRanda,
   TargetScan, RNAhybrid).
2. Keep a miRanda record when its binding site contains a TargetScan site,
   and an RNAhybrid record when its site contains a TargetScan site.
3. For pairs surviving both filters, slide a 6-nt window across the miRanda
   site and emit the interaction when that window falls inside the RNAhybrid
   site.

Outputs (written to --output-dir):
  * three_overlap.txt         lnc  miRNA  score  energy  start  end
  * three_overlap_unique.txt  de-duplicated version
"""

import argparse
import os

from common import interaction_key, load_interactions


def predict(miranda, targetscan, rnahybrid, out_dir):
    m = load_interactions(miranda)
    t = load_interactions(targetscan)
    r = load_interactions(rnahybrid)

    common = set(m) & set(t) & set(r)

    # pairwise containment: TargetScan site within the tool's site
    m_ts = {}   # id -> [(score, energy, start, end)]
    rh_ts = {}  # id -> [(start, end, energy)]
    for id_ in common:
        for (score, energy, s, e) in m[id_]:
            s_i, e_i = int(s), int(e)
            if any(int(ts_s) >= s_i and int(ts_e) <= e_i
                   for (_, _, ts_s, ts_e) in t[id_]):
                m_ts.setdefault(id_, []).append((score, energy, s, e))
        for (_, energy, s, e) in r[id_]:
            s_i, e_i = int(s), int(e)
            if any(int(ts_s) >= s_i and int(ts_e) <= e_i
                   for (_, _, ts_s, ts_e) in t[id_]):
                rh_ts.setdefault(id_, []).append((s, e, energy))

    # sliding 6-nt window of the miRanda site within the RNAhybrid site
    records = []
    for id_ in common:
        if id_ not in m_ts or id_ not in rh_ts:
            continue
        miRNA, lnc = id_.split("!", 1)
        for (score, energy, m_s, m_e) in m_ts[id_]:
            end_final = int(m_e) - 8
            start_tmp = int(m_s) - 1
            for (rh_s, rh_e, _) in rh_ts[id_]:
                rh_s_i, rh_e_i = int(rh_s), int(rh_e)
                for k in range(start_tmp, end_final):
                    new_start = k + 1
                    new_end = new_start + 5
                    if new_start >= rh_s_i and new_end <= rh_e_i:
                        records.append((lnc, miRNA, score, energy, m_s, m_e))
                        break

    seen = set()
    with open(os.path.join(out_dir, "three_overlap.txt"), "w") as fw:
        for rec in records:
            fw.write("\t".join(rec) + "\n")
    with open(os.path.join(out_dir, "three_overlap_unique.txt"), "w") as fw:
        for rec in records:
            if rec not in seen:
                seen.add(rec)
                fw.write("\t".join(rec) + "\n")
    return len(seen)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--miranda", required=True)
    ap.add_argument("--targetscan", required=True)
    ap.add_argument("--rnahybrid", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    n = predict(args.miranda, args.targetscan, args.rnahybrid, args.output_dir)
    print(f"done: {n} unique ceRNA interactions -> {args.output_dir}")


if __name__ == "__main__":
    main()
