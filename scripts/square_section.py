"""Render a CVM quaternary square plot from CEWorkbench data.

Reproduces the "square plot" style of Fig. 20, Jindal & Lele (2025),
CALPHAD 89, 102825: a filled contour of G/H/S (or pair SRO) over the unit
(X,Y) square, where (X,Y) parametrizes a 2-D slice of the quaternary
composition simplex (see QuaternarySquareScan.java for the mapping and the
region taxonomy: interior / a true binary solve along each square edge / an
analytic pure-element value at each square corner).

Unlike the ternary isothermal section, this needs no special ternary-axis
library -- the (X,Y) square is already Cartesian, so a plain
matplotlib.pyplot.tricontourf over the scan's actual (possibly irregular)
point set is used directly.

Input mode (fast path only; no slow "compute via API" path -- the square
scan has no natural per-point standalone CLI call the way a single ternary
point does, since a caller must already have decided which two of the four
elements it wants a mole fraction axis over):

  --from-json FILE   Read pre-computed grid data written by the Java
                      `quaternary_square` CLI subcommand. Java computes,
                      this script only renders.

Usage:
    java -cp ... org.ce.ui.cli.Main quaternary_square < request.json > grid.json
    python scripts/square_section.py --from-json grid.json
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_from_json(path):
    """Reads grid data computed by the Java `quaternary_square` CLI
    subcommand. Returns (elements, variant, structure, temperature,
    calculation, x, y, values, corner_labels)."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not data.get("ok"):
        raise RuntimeError(f"quaternary_square response was not ok: {data.get('message')}")

    elements = data["elements"]
    x, y, values = [], [], []
    nan_count = 0
    for p in data["points"]:
        # Jackson serializes Double.NaN as the JSON string "NaN" (raw NaN
        # isn't valid JSON) -- e.g. SRO alpha undefined at a corner. Skip
        # rather than corrupt the array dtype (mirrors isothermal_section.py).
        v = p["value"]
        if isinstance(v, str) or (isinstance(v, float) and v != v):
            nan_count += 1
            continue
        x.append(p["x"])
        y.append(p["y"])
        values.append(v)
    if nan_count:
        print(f"[info] excluded {nan_count} grid points with undefined (NaN) value", file=sys.stderr)
    if data.get("skipped"):
        print(f"[warn] {data['skipped']} grid points skipped upstream (Java-side)", file=sys.stderr)

    return (elements, data["variant"], data["structure"], data["temperature"], data["calculation"],
            np.array(x), np.array(y), np.array(values))


def render_square_plot(elements, variant, structure, temperature, quantity_label,
                        x, y, values, out_path):
    """Renders the filled square contour and saves it to out_path.

    Corner labels follow QuaternarySquareScan's slot convention: (X,Y)=(0,0)
    is pure elements[1] ("Ti"-role slot1), (0,1) is pure elements[0]
    ("Nb"-role slot0), (1,0) is pure elements[2] ("V"-role slot2), (1,1) is
    pure elements[3] ("Zr"-role slot3) -- see QuaternarySquareScan.java's
    toMoleFractions doc for the underlying formula this mirrors.
    """
    if len(x) < 3:
        raise RuntimeError("Not enough grid points to plot.")

    fig, ax = plt.subplots(figsize=(5.5, 4.8))

    vmin, vmax = np.percentile(values, [2, 98])
    # Round contour/colorbar levels to a "nice" step sized to the data's own
    # range -- e.g. multiples of 1000 for a Gibbs-energy plot spanning tens of
    # thousands of J/mol, but multiples of 0.05 for an SRO plot spanning
    # roughly [-1,1]. matplotlib.ticker.MaxNLocator already implements this
    # (steps of 1/2/2.5/5 x10^n, chosen from the data range), so it's used
    # directly rather than hand-rolling a quantity-specific rule.
    levels = plt.MaxNLocator(nbins=10).tick_values(vmin, vmax)
    step = levels[1] - levels[0]

    triangulation = mtri.Triangulation(x, y)
    cs = ax.tricontourf(triangulation, values, levels=levels, cmap="RdYlBu_r", extend="both")
    ax.tricontour(triangulation, values, levels=levels, colors="k", linewidths=0.3, alpha=0.4)

    units = {"GIBBS_ENERGY": "J/mol", "ENTHALPY": "J/mol", "ENTROPY": "J/mol·K"}
    unit = units.get(quantity_label, "")
    unit_suffix = f" ({unit})" if unit else ""
    cbar = fig.colorbar(cs, ax=ax, label=f"{quantity_label}{unit_suffix}", shrink=0.85)
    # Show decimals only if the chosen step actually needs them (e.g. SRO's
    # sub-1 steps), so an energy plot's ticks stay as plain "-2,000" rather
    # than "-2,000.0".
    decimals = max(0, -int(np.floor(np.log10(step))) ) if step < 1 else 0
    cbar.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.{decimals}f}"))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # Corner annotations: which pure element sits at each (X,Y) corner, placed
    # outside the plot area (in points, not data coordinates, so the offset is
    # a fixed screen distance regardless of axis scale) so they never overlap
    # the contour fill or the axis tick labels.
    corner_labels = {
        (0.0, 0.0): elements[1],
        (0.0, 1.0): elements[0],
        (1.0, 0.0): elements[2],
        (1.0, 1.0): elements[3],
    }
    offsets = {
        (0.0, 0.0): (-14, -14, "right", "top"),
        (0.0, 1.0): (-14, 10, "right", "bottom"),
        (1.0, 0.0): (14, -14, "left", "top"),
        (1.0, 1.0): (14, 10, "left", "bottom"),
    }
    for (cx, cy), label in corner_labels.items():
        dx, dy, ha, va = offsets[(cx, cy)]
        ax.annotate(label, (cx, cy), xycoords="data",
                    textcoords="offset points", xytext=(dx, dy),
                    ha=ha, va=va, fontsize=9, fontweight="bold",
                    annotation_clip=False)

    title_variant = "" if variant == "STANDARD" else f" [{variant}]"
    # Extra pad clears the top-corner element labels, which sit above the axes.
    ax.set_title(f"{'-'.join(elements)} {structure} square section at {temperature:g} K"
                 f" — {quantity_label}{title_variant}", pad=26)

    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-json", required=True,
                     help="Path to JSON grid data from the Java `quaternary_square` CLI subcommand.")
    ap.add_argument("--out", default=None,
                     help="output PNG path; defaults to scripts/<elements>_<variant>_<T>K_<quantity>.png")
    args = ap.parse_args()

    elements, variant, structure, temperature, quantity_label, x, y, values = load_from_json(args.from_json)

    out = args.out or os.path.join(
        REPO_ROOT, "scripts", f"{'-'.join(elements)}_{variant}_{temperature:g}K_{quantity_label}.png")
    render_square_plot(elements, variant, structure, temperature, quantity_label, x, y, values, out)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
