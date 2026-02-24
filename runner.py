#!/usr/bin/env python3

import subprocess
import sys

configs = [
    {
        "name": "NLOBK_MV_smallest",
        "datafile": "bkdipoles/mv_bk.dat",
        "C2": "11.8",
        "charm_mass": "1.04",
        "proton_area": "23.5",
        "rc_scheme": "SMALLEST",
        "order": "NLO"
    },
    {
        "name": "KCBK_parent",
        "datafile": "bkdipoles/pd_bk.dat",
        "C2": "663",
        "charm_mass": "1.4",
        "proton_area": "20.7",
        "rc_scheme": "PARENT",
        "order": "NLO"
    },
    {
        "name": "KCBK_smallest",
        "datafile": "bkdipoles/bk_map.dat",
        "C2": "1.7",
        "charm_mass": "1.25",
        "proton_area": "8.75",
        "rc_scheme": "SMALLEST",
        "order": "NLO"
    },
    {
        "name": "NLOBK_MVgamma_smallest",
        "datafile": "bkdipoles/mvgam_bk.dat",
        "C2": "1314.9257306",
        "charm_mass": "1.2049379",
        "proton_area": "22.9017918",
        "rc_scheme": "SMALLEST",
        "order": "NLO"
    },
    {
        "name": "NLOBK_MVgamma_parent",
        "datafile": "bkdipoles/pd_nlo_bk.dat",
        "C2": str(10**3.88),
        "charm_mass": "1.20",
        "proton_area": "24.3",
        "rc_scheme": "PARENT",
        "order": "NLO"
    },
    {
        "name": "LO MV1",
        "datafile": "bkdipoles/lo_mv.dat",
        "C2": "14.5",
        "charm_mass": "9999999",
        "proton_area": "18.81",
        "rc_scheme": "PARENT",
        "order": "LO"
    }
]

if len(sys.argv) < 2:
    print("Usage: python runner.py /path/to/nlodis")
    sys.exit(1)

exe = sys.argv[1]

for idx, cfg in enumerate(configs, start=1):
    args = [exe]
    
    if idx != 1:
        args.append("--no-header")
    
    args.extend([
        "--name", cfg["name"],
        "--datafile", cfg["datafile"],
        "--C2", cfg["C2"],
        "--charm_mass", cfg["charm_mass"],
        "--proton_area", cfg["proton_area"],
        "--rc_scheme", cfg["rc_scheme"],
        "--order", cfg["order"],
        "--mcintpoints", str(6e5),
        "--runmode", "HERA_FL"
    ])

    if cfg["order"] == "NLO":
        continue
    
    print(f"Running config {idx}/{len(configs)}: {cfg['name']} with args {args}", file=sys.stderr)
    
    proc = subprocess.run(args, text=True, capture_output=True)
    if proc.returncode != 0:
        print(f"Error running {cfg['name']}:", file=sys.stderr)
        print(proc.stderr, file=sys.stderr)
        sys.exit(proc.returncode)
    
    with open(f"LO_output_{idx}", "w") as f:
        f.write(proc.stdout)
    
    print(f"Output written to output_{idx}", file=sys.stderr)

print("All runs completed successfully!", file=sys.stderr)