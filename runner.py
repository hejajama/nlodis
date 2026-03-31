#!/usr/bin/env python3

import subprocess
import sys



configs = [
    { # 2506.00487  Bal+SD coupling, MAP parameters
         "name": "KCBK_smallest",
         "datafile": "nlobkdatafiles/zenodo.15552940/balsd/bk_map.dat",
         "C2": "1.74",
         "charm_mass": "1.24",
         "proton_area": "9.08", # in mb
         "rc_scheme": "SMALLEST",
         "order": "NLO"
    },
]

if len(sys.argv) < 2:
    print("Usage: python runner.py /path/to/nlodis")
    sys.exit(1)

exe = sys.argv[1]
epsrel = 0.001 # relative accuracy goal for Cuba integration
mcintpoints = "2e6" # default number of Monte Carlo points for Cuba integration
light_quark_mass = 0.005 # in GeV
cubamethod = "vegas"


for idx, cfg in enumerate(configs, start=1):
    args = [exe]

    args.extend([
        "--name", cfg["name"],
        "--datafile", cfg["datafile"],
        "--C2", cfg["C2"],
        "--charm_mass", cfg["charm_mass"],
        "--proton_area", cfg["proton_area"],
        "--rc_scheme", cfg["rc_scheme"],
        "--order", cfg["order"],
        "--mcintpoints", mcintpoints,
        "--runmode", "F2FL_GRID",
        "--cubamethod", cubamethod,
        "--epsrel", str(epsrel),
        "--light_mass", str(light_quark_mass),
        "--nf", "4"
    ])

    print(f"Running config {idx}/{len(configs)}: {cfg['name']} with args {args}", file=sys.stderr)

    fname = (
        f"f2fl_grid_mc_{mcintpoints}_config_{cfg['name'].replace(' ', '_')}"
    )
    print(f"Output will be written to {fname}", file=sys.stderr)
    with open(fname, "w") as f:
        proc = subprocess.Popen(
            args,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=None,
            bufsize=1,
        )

        assert proc.stdout is not None
        for line in proc.stdout:
            f.write(line)
            f.flush()

        proc.wait()

    if proc.returncode != 0:
        print(f"Error running {cfg['name']}. Partial output written to {fname}", file=sys.stderr)
        sys.exit(proc.returncode)

    print(f"Output written to {fname}", file=sys.stderr)

print("All runs completed successfully!", file=sys.stderr)