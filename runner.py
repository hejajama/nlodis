#!/usr/bin/env python3

import subprocess
import sys



configs = [
    # {
    # "name":"NLOBK_MV_smallest_0",
    # "datafile":"bkdata/nlobk_mv_balsd/sample_0",
    # "C2":"11.670",
    # "proton_area":"23.670",
    # "charm_mass":"1.036",
    # "rc_scheme":"SMALLEST",
    # "order":"NLO"
    # },
    {
        "name": "NLOBK_MV_smallest",
        "datafile": "data/balsd/median_bk.dat",
        "C2": "11.834",
        "charm_mass": "1.0399",
        "proton_area": "23.873",
        "rc_scheme": "SMALLEST",
        "order": "NLO"
    },
    # {
    #     "name": "KCBK_parent",
    #     "datafile": "bkdata/pd_bk_map.dat",
    #     "C2": "663",
    #     "charm_mass": "1.4",
    #     "proton_area": "20.7",
    #     "rc_scheme": "PARENT",
    #     "order": "NLO"
    # },
    # {
    #     "name": "KCBK_smallest",
    #     "datafile": "bkdata/balsd_bk_map.dat",
    #     "C2": "1.7",
    #     "charm_mass": "1.25",
    #     "proton_area": "8.75",
    #     "rc_scheme": "SMALLEST",
    #     "order": "NLO"
    # },
    # {
    #     "name": "NLOBK_MVgamma_smallest",
    #     "datafile": "bkdata/mvgam_bk.dat",
    #     "C2": "1314.9257306",
    #     "charm_mass": "1.2049379",
    #     "proton_area": "22.9017918",
    #     "rc_scheme": "SMALLEST",
    #     "order": "NLO"
    # },
    # {
    #     "name": "NLOBK_MVgamma_parent",
    #     "datafile": "bkdata/pd_nlo_bk.dat",
    #     "C2": str(10**3.88),
    #     "charm_mass": "1.20",
    #     "proton_area": "24.3",
    #     "rc_scheme": "PARENT",
    #     "order": "NLO"
    # },
    # {
    #     "name": "LO MV1",
    #     "datafile": "bkdata/lo_mv.dat",
    #     "C2": "14.5",
    #     "charm_mass": "9999999",
    #     "proton_area": "18.81",
    #     "rc_scheme": "PARENT",
    #     "order": "LO"
    # }
]

if len(sys.argv) < 2:
    print("Usage: python runner.py /path/to/nlodis")
    sys.exit(1)

exe = sys.argv[1]

for mcintpoints in ["1e6","2e6"]:
    for cubamethod in ["vegas"]: #"vegas","suave","divonne"]:
        for epsrel in [0.01,0.001,0.0001]:
            print(f"Running with mcintpoints={mcintpoints} and cubamethod={cubamethod} and epsrel={epsrel}", file=sys.stdout)

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
                    "--mcintpoints", mcintpoints,
                    "--runmode", "F2FL_GRID",
                    "--cubamethod", cubamethod,
                    "--epsrel", str(epsrel),
                ])

                #if cfg["order"] == "NLO":
                #    continue
                
                print(f"Running config {idx}/{len(configs)}: {cfg['name']} with args {args}", file=sys.stderr)

                fname = f"convergence/mvmedian/f2fl_mc_{mcintpoints}_{cubamethod}_{cfg['name'].replace(' ', '_')}_epsrel_{str(epsrel)}"
                print(f"Output will be written to {fname}", file=sys.stderr)
                with open(fname, "w") as f:
                    proc = subprocess.Popen(
                        args,
                        text=True,
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