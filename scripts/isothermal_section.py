"""Render a CVM ternary isothermal section from CEWorkbench data.

Two input modes:

  --from-json FILE   Read pre-computed grid data written by the Java
                      `ternary_grid` CLI subcommand (fast path — this is what
                      the Swing GUI panel uses; Java computes, this script
                      only renders).

  --elements/--temperature/--n (no --from-json)
                      Compute the grid itself by calling the `api` CLI
                      subcommand once per point (slow path — one JVM start per
                      grid point; useful standalone/for debugging, not what
                      the GUI uses).

Renders a filled contour of the requested quantity on a ternary triangle via
mpltern, in the style of Figs. 15/20 of Jindal & Lele (2025), CALPHAD 89,
102825.

Usage:
    # Fast path (Java already computed the grid):
    java -cp ... org.ce.ui.cli.Main ternary_grid < request.json > grid.json
    python scripts/isothermal_section.py --from-json grid.json --quantity gibbsEnergy

    # Slow path (script computes via repeated `api` calls):
    python scripts/isothermal_section.py --elements Nb-Ti-V --temperature 1273 \
        --quantity gibbsEnergy --n 21 --classpath-jars <jackson jars...>
"""
import argparse
import glob
import json
import os
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt
import mpltern  # noqa: F401  (registers the 'ternary' projection)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN_CLASS = "org.ce.ui.cli.Main"


def find_jackson_jars():
    gradle_cache = os.path.expanduser("~/.gradle/caches/modules-2/files-2.1/com.fasterxml.jackson.core")
    jars = []
    for name in ("jackson-databind", "jackson-core", "jackson-annotations"):
        pattern = os.path.join(gradle_cache, name, "*", "*", f"{name}-*.jar")
        candidates = [p for p in glob.glob(pattern) if "sources" not in p and "javadoc" not in p]
        if not candidates:
            raise FileNotFoundError(f"Could not locate {name} jar under {gradle_cache}")
        jars.append(candidates[0])
    return jars


def build_classpath(extra_jars):
    classes_dir = os.path.join(REPO_ROOT, "build", "classes", "java", "main")
    jars = extra_jars if extra_jars else find_jackson_jars()
    sep = ";" if os.name == "nt" else ":"
    return sep.join([classes_dir] + jars)


def call_api(classpath, request):
    payload = json.dumps(request).encode("utf-8")
    proc = subprocess.run(
        ["java", "-cp", classpath, MAIN_CLASS, "api"],
        input=payload, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        cwd=REPO_ROOT,
    )
    try:
        return json.loads(proc.stdout.decode("utf-8"))
    except json.JSONDecodeError as e:
        raise RuntimeError(
            f"Non-JSON response (exit {proc.returncode}): "
            f"{proc.stdout[:500]!r} / stderr: {proc.stderr[-1000:].decode('utf-8', 'replace')}"
        ) from e


def ternary_grid(n):
    """Barycentric grid over the 2-simplex with n subdivisions per edge."""
    points = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            points.append((i / n, j / n, k / n))
    return points


def nearest_edge_and_interior(fa, fb, fc):
    """For composition (fa, fb, fc), find the nearest edge (component -> 0)
    and a point further along the same ray into the interior.

    Returns (edge_point, interior_point, t) where t in (0, 1) is the
    barycentric fraction of the original point between edge_point (t=0)
    and interior_point (t=1) along the ray from the edge vertex-opposite
    through the original point.
    """
    comps = [fa, fb, fc]
    axis = min(range(3), key=lambda i: comps[i])  # component closest to 0
    t_edge = comps[axis]

    def scale(target_axis_value):
        """Rescale composition so comps[axis] == target_axis_value, keeping
        the ratio of the other two components fixed (moves along the ray
        toward/away from the edge where axis == 0)."""
        others = [i for i in range(3) if i != axis]
        remaining = 1.0 - target_axis_value
        other_sum = comps[others[0]] + comps[others[1]]
        if other_sum <= 0:
            return None
        scaled = list(comps)
        scaled[axis] = target_axis_value
        for i in others:
            scaled[i] = comps[i] / other_sum * remaining
        return tuple(scaled)

    edge_point = scale(0.0)
    # step further into the interior than the failed point itself
    interior_target = min(t_edge + 0.08, 1.0)
    interior_point = scale(interior_target)
    if interior_point is None or interior_target <= t_edge:
        return edge_point, None, None
    t = t_edge / interior_target
    return edge_point, interior_point, t


def compute_via_api(args, elements):
    """Slow path: one `api` subprocess call per grid point."""
    classpath = build_classpath(args.classpath_jars)
    grid = ternary_grid(args.n)

    def request_for(fa, fb, fc):
        composition = {elements[1]: fb, elements[2]: fc}
        return {
            "system": {
                "elements": args.elements,
                "structure": args.structure,
                "model": args.model,
                "engine": args.engine,
            },
            "calculation": args.calculation,
            "conditions": {
                "temperature": args.temperature,
                "composition": composition,
            },
        }

    def query(fa, fb, fc):
        resp = call_api(classpath, request_for(fa, fb, fc))
        if not resp.get("ok"):
            return None
        points = resp.get("points", [])
        if not points or args.quantity not in points[0]:
            return None
        if points[0].get("converged") is False:
            return None
        return points[0][args.quantity]

    ta, la, ra, values = [], [], [], []
    skipped = 0
    interpolated = 0
    for (fa, fb, fc) in grid:
        value = query(fa, fb, fc)
        if value is None:
            edge_pt, interior_pt, t = nearest_edge_and_interior(fa, fb, fc)
            edge_val = query(*edge_pt) if edge_pt else None
            interior_val = query(*interior_pt) if interior_pt else None
            if edge_val is not None and interior_val is not None:
                value = (1 - t) * edge_val + t * interior_val
                interpolated += 1
            else:
                skipped += 1
                continue
        ta.append(fc)
        la.append(fa)
        ra.append(fb)
        values.append(value)

    if interpolated:
        print(f"[info] interpolated {interpolated}/{len(grid)} near-edge grid points", file=sys.stderr)
    if skipped:
        print(f"[warn] skipped {skipped}/{len(grid)} grid points", file=sys.stderr)

    return np.array(ta), np.array(la), np.array(ra), np.array(values)


def load_from_json(path):
    """Fast path: read grid data already computed by the Java `ternary_grid`
    CLI subcommand. Returns (elements, structure, temperature, calculation,
    ta, la, ra, values)."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not data.get("ok"):
        raise RuntimeError(f"ternary_grid response was not ok: {data.get('message')}")

    elements = data["elements"]
    ta, la, ra, values = [], [], [], []
    nan_count = 0
    for p in data["points"]:
        # Jackson serializes Double.NaN as the JSON string "NaN" (raw NaN isn't
        # valid JSON) -- e.g. SRO alpha is undefined where a pair's mole
        # fraction is exactly 0. Skip these rather than corrupt the array dtype.
        v = p["value"]
        if isinstance(v, str) or (isinstance(v, float) and v != v):
            nan_count += 1
            continue
        la.append(p[elements[0]])
        ra.append(p[elements[1]])
        ta.append(p[elements[2]])
        values.append(v)
    if nan_count:
        print(f"[info] excluded {nan_count} grid points with undefined (NaN) value "
              f"(e.g. SRO where a pair's mole fraction is 0)", file=sys.stderr)
    if data.get("skipped"):
        print(f"[warn] {data['skipped']} grid points skipped upstream (Java-side)", file=sys.stderr)

    return (elements, data["structure"], data["temperature"], data["calculation"],
            np.array(ta), np.array(la), np.array(ra), np.array(values))


def render_ternary_plot(elements, structure, temperature, quantity_label,
                         ta, la, ra, values, out_path):
    """Renders the filled ternary contour and saves it to out_path."""
    if len(ta) < 3:
        raise RuntimeError("Not enough grid points to plot.")

    fig = plt.figure(figsize=(5.5, 4.8))
    ax = fig.add_subplot(projection="ternary")

    cs = ax.tricontourf(ta, la, ra, values, levels=20, cmap="RdYlBu_r")
    ax.tricontour(ta, la, ra, values, levels=20, colors="k", linewidths=0.3, alpha=0.4)
    # SRO (Cowley-Warren alpha) is dimensionless; G/H/S are in J/mol.
    unit = "" if quantity_label.startswith("SRO") else " (J/mol)"
    # pad pushes the colorbar clear of the r-axis tick labels, which extend
    # outward past the triangle's right edge.
    fig.colorbar(cs, ax=ax, label=f"{quantity_label}{unit}", shrink=0.8, pad=0.15)

    # mpltern rotates axis-name labels to match each edge's slant by default
    # (label_rotation_mode='axis'); 'horizontal' keeps them plain, upright text.
    ax.taxis.set_label_rotation_mode("horizontal")
    ax.laxis.set_label_rotation_mode("horizontal")
    ax.raxis.set_label_rotation_mode("horizontal")
    ax.set_tlabel(elements[2])
    ax.set_llabel(elements[0])
    ax.set_rlabel(elements[1])
    ax.taxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.laxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.raxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.grid(linewidth=0.4, alpha=0.5)

    ax.set_title(f"{'-'.join(elements)} {structure} isothermal section at {temperature:g} K — {quantity_label}",
                 pad=30)

    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-json", default=None,
                     help="Path to JSON grid data from the Java `ternary_grid` CLI subcommand (fast path).")
    ap.add_argument("--elements", help="e.g. Nb-Ti-V (required if not --from-json)")
    ap.add_argument("--structure", default="BCC_A2")
    ap.add_argument("--model", default="T")
    ap.add_argument("--engine", default="CVM")
    ap.add_argument("--temperature", type=float, help="required if not --from-json")
    ap.add_argument("--calculation", default="GIBBS_ENERGY")
    ap.add_argument("--quantity", default="gibbsEnergy", help="JSON field to plot (slow path only): gibbsEnergy, enthalpy, entropy")
    ap.add_argument("--n", type=int, default=20, help="grid subdivisions per triangle edge (slow path only)")
    ap.add_argument("--classpath-jars", nargs="*", default=None, help="explicit jackson jar paths; auto-detected if omitted")
    ap.add_argument("--out", default=None, help="output PNG path; defaults to scripts/<elements>_<T>K_<quantity>.png")
    args = ap.parse_args()

    if args.from_json:
        elements, structure, temperature, calculation, ta, la, ra, values = load_from_json(args.from_json)
        quantity_label = calculation
    else:
        if not args.elements or args.temperature is None:
            sys.exit("--elements and --temperature are required unless --from-json is given.")
        elements = args.elements.split("-")
        if len(elements) != 3:
            sys.exit(f"--elements must name exactly 3 elements separated by '-', got {args.elements!r}")
        structure = args.structure
        temperature = args.temperature
        quantity_label = args.quantity
        ta, la, ra, values = compute_via_api(args, elements)

    out = args.out or os.path.join(REPO_ROOT, "scripts", f"{'-'.join(elements)}_{temperature:g}K_{quantity_label}.png")
    render_ternary_plot(elements, structure, temperature, quantity_label, ta, la, ra, values, out)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
