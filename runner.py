#!/usr/bin/env python3

import subprocess
import sys

configs = [
    ("NLOBK_MV_smallest", "/Users/hejajama/Downloads/mv_bk.dat", "11.8", "1.04", "23.5", "SMALLEST"),
    ("KCBK_parent", "/Users/hejajama/code/nlodisfit_bayesian/data/pd/bk_map.dat", "663", "1.4", "20.7", "PARENT"),
    ("KCBK_smallest", "/Users/hejajama/code/nlodisfit_bayesian/data/balsd/bk_map.dat", "1.7", "1.25", "8.75", "SMALLEST"),
    ("NLOBK_MVgamma smallest", "/Users/hejajama/Downloads/mvgam_bk.dat", "1314.9257306", "1.2049379", "22.9017918", "SMALLEST"),
    ("NLOBK_MVgamma parent", "/Users/hejajama/Downloads/pd_nlo_bk.dat", str(10**3.88), "1.20", "24.3", "PARENT"),
]

if len(sys.argv) < 2:
    print("Usage: python run_configs.py /path/to/nlodis")
    sys.exit(1)

exe = sys.argv[1]

for idx, cfg in enumerate(configs, start=1):
    args = [exe]
    if idx != 1:
        args.append("--no-header")
    args.extend(cfg)
    proc = subprocess.run(args, text=True, capture_output=True)
    if proc.stderr:
        print(f"Warnings (std::cerr output) from config {idx}:", file=sys.stderr)
        print(proc.stderr, file=sys.stderr)

    if proc.returncode != 0:
        print(proc.stderr, file=sys.stderr)
        sys.exit(proc.returncode)
    with open(f"output_{idx}", "w") as f:
        f.write(proc.stdout)