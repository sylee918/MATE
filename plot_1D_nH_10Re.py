#!/usr/bin/env python3
"""
plot_1D_nH_10Re.py
=============================================================================
Purpose: Compare MATE exospheric neutral hydrogen density (nH) simulations
         at 3 Re and 10 Re across DOY 164-174 (2008 storm event).

Compared Runs:
  1. With Charge Exchange (With CX):  MATE_nH_GRCPX1_test_2008*.data
  2. Without Charge Exchange (No CX): MATE_nH_GRCP_RCCX2_2008*.data

Usage:
  python plot_1D_nH_10Re.py
=============================================================================
"""

import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

def main():
    year = 2008
    doys = list(range(164, 175)) # 164 to 174

    runs = [
        {
            "id": "With_CX",
            "name": "With CX (GRCPX1_RCCX2)",
            "prefix": "MATE_nH_GRCPX1_RCCX2_",
            "color": "#0072bd",
            "linestyle": "-",
            "linewidth": 2.0
        },
        {
            "id": "Without_CX",
            "name": "Without CX (GRCP_RCCX2)",
            "prefix": "MATE_nH_GRCP_RCCX2_",
            "color": "#d9531e",
            "linestyle": "--",
            "linewidth": 2.0
        }
    ]

    candidate_dirs = [
        "C:\\Users\\slee122\\OneDrive - NASA\\Desktop\\Work\\git\\MATE",
        "\\\\wsl.localhost\\Ubuntu-22.04\\home\\sylee\\exospherecode\\MATE\\output\\0728",
        "/mnt/c/Users/slee122/OneDrive - NASA/Desktop/Work/git/MATE",
        "/home/sylee/exospherecode/MATE/output/0728",
        "./"
    ]

    sim_data = []

    for r_cfg in runs:
        prefix = r_cfg["prefix"]
        data_dir = None
        for cdir in candidate_dirs:
            if not os.path.isdir(cdir):
                continue
            test_f = os.path.join(cdir, f"{prefix}{year}{doys[0]:03d}.data")
            if os.path.isfile(test_f):
                data_dir = cdir
                break

        if not data_dir:
            print(f"Warning: Could not find files for prefix {prefix}")
            continue

        print(f"\nProcessing {r_cfg['name']} in {data_dir}...")
        all_time, all_3Re, all_10Re = [], [], []

        for doy in doys:
            fpath = os.path.join(data_dir, f"{prefix}{year}{doy:03d}.data")
            if not os.path.isfile(fpath):
                continue
            raw = np.fromfile(fpath, dtype=np.float32)
            n_floats = len(raw)

            if n_floats == 1086912:
                # [17, 72, 37, 24]
                data4D = raw.reshape((24, 37, 72, 17))
                iLat_eq = 18
                nH_3Re = data4D[:, iLat_eq, 0, 2]   # index 2 is 3 Re
                nH_10Re = data4D[:, iLat_eq, 0, 16] # index 16 is 10 Re
                t_vec = doy + np.arange(24) / 24.0
            elif n_floats == (17 * 24):
                data2D = raw.reshape((24, 17))
                nH_3Re = data2D[:, 2]
                nH_10Re = data2D[:, 16]
                t_vec = doy + np.arange(24) / 24.0
            else:
                continue

            all_time.extend(t_vec)
            all_3Re.extend(nH_3Re)
            all_10Re.extend(nH_10Re)

        all_time = np.array(all_time)
        all_3Re = np.array(all_3Re)
        all_10Re = np.array(all_10Re)

        # Exclude 1st time step (initial transient)
        if len(all_time) > 1:
            all_time = all_time[1:]
            all_3Re = all_3Re[1:]
            all_10Re = all_10Re[1:]

        r_cfg["time"] = all_time
        r_cfg["nH_3Re"] = all_3Re
        r_cfg["nH_10Re"] = all_10Re
        sim_data.append(r_cfg)

        print(f"  3 Re : Min={all_3Re.min():.2f}, Max={all_3Re.max():.2f}, Mean={all_3Re.mean():.2f} cm^-3")
        print(f" 10 Re : Min={all_10Re.min():.2f}, Max={all_10Re.max():.2f}, Mean={all_10Re.mean():.2f} cm^-3")

    if not sim_data:
        print("No simulation data loaded!")
        return

    # Load Dst
    has_dst = False
    dst_time, dst_val = None, None
    dst_candidates = [
        "2008164_Dst.txt",
        "C:\\Users\\slee122\\OneDrive - NASA\\Documents\\MATLAB\\2008164_Dst.txt",
        "C:\\Users\\slee122\\OneDrive - NASA\\Desktop\\Work\\git\\MATE\\2008164_Dst.txt",
    ]
    for dpath in dst_candidates:
        if os.path.isfile(dpath):
            try:
                dst_data = np.loadtxt(dpath, skiprows=1)
                dst_time = dst_data[:, 1] + dst_data[:, 2] / 24.0
                dst_val = dst_data[:, 3]
                has_dst = True
                break
            except Exception:
                pass

    # Plot
    fig, axes = plt.subplots(3 if has_dst else 2, 1, figsize=(11, 8.5), sharex=True)
    plt.subplots_adjust(hspace=0.15)

    x_min = np.floor(sim_data[0]["time"].min())
    x_max = sim_data[0]["time"].max()

    ax_idx = 0
    if has_dst:
        axes[ax_idx].plot(dst_time, dst_val, color="#222222", linewidth=2.0)
        axes[ax_idx].axhline(0, color="gray", linestyle="--", alpha=0.5)
        axes[ax_idx].set_ylabel("Dst [nT]", fontsize=12, fontweight="bold")
        axes[ax_idx].set_title(f"Geomagnetic Activity (Dst) & MATE 1D Exospheric H Density: CX Effect ({year})",
                              fontsize=13, fontweight="bold")
        axes[ax_idx].grid(True, linestyle=":", alpha=0.6)
        axes[ax_idx].set_xlim(x_min, x_max)
        ax_idx += 1

    # 3 Re panel
    for s in sim_data:
        axes[ax_idx].plot(s["time"], s["nH_3Re"], label=s["name"],
                          color=s["color"], linestyle=s["linestyle"], linewidth=s["linewidth"])
    axes[ax_idx].set_ylabel(r"$n_H\ \mathrm{[cm^{-3}]}$ at 3 $R_E$", fontsize=12, fontweight="bold")
    axes[ax_idx].legend(loc="upper left", frameon=False, fontsize=10)
    axes[ax_idx].grid(True, linestyle=":", alpha=0.6)
    axes[ax_idx].set_xlim(x_min, x_max)
    ax_idx += 1

    # 10 Re panel
    for s in sim_data:
        axes[ax_idx].plot(s["time"], s["nH_10Re"], label=s["name"],
                          color=s["color"], linestyle=s["linestyle"], linewidth=s["linewidth"])
    axes[ax_idx].set_ylabel(r"$n_H\ \mathrm{[cm^{-3}]}$ at 10 $R_E$", fontsize=12, fontweight="bold")
    axes[ax_idx].set_xlabel(f"Day of Year (DOY in {year})", fontsize=12, fontweight="bold")
    axes[ax_idx].legend(loc="upper left", frameon=False, fontsize=10)
    axes[ax_idx].grid(True, linestyle=":", alpha=0.6)
    axes[ax_idx].set_xlim(x_min, x_max)

    out_fig = "MATE_1D_nH_CX_comparison_3Re_10Re.png"
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    print(f"\nPlot saved to {out_fig}")

if __name__ == "__main__":
    main()
