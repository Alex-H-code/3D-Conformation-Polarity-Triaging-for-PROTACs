"""
plot_molecular_properties.py
==============================
Renders a multi-panel box-plot + jittered-strip figure (one panel per
property, independent y-axis per panel) as a self-contained SVG file --
the same style used for the compound-shortlist property summary figure.

Usage:
    python plot_molecular_properties.py \\
        --input-csv filtered_compounds_with_adme.csv \\
        --columns degrader_mw degrader_hba degrader_hbd tpsa_2d psa_averaged degrader_rotatable_bonds \\
        --labels MW HBA HBD "TPSA (2D)" "3D-PSA (avg)" NRotB \\
        --output-svg molecular_properties.svg \\
        --title "Molecular Properties"
"""
import argparse
import math
import random


# Fixed-order categorical palette (see dataviz skill palette.md) -- assign in
# this order for up to 8 panels; cycles only if more than 8 columns are given.
PALETTE = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']

INK = '#0b0b0b'
INK_SEC = '#52514e'
MUTED = '#898781'
GRID = '#e1e0d9'
BASELINE = '#c3c2b7'
SURFACE = '#fcfcfb'


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-csv", required=True)
    ap.add_argument("--columns", nargs="+", required=True, help="Column names to plot, one panel each")
    ap.add_argument("--labels", nargs="+", default=None, help="Panel titles (default: column names)")
    ap.add_argument("--output-svg", required=True)
    ap.add_argument("--title", default="Molecular Properties")
    ap.add_argument("--panel-width", type=float, default=175)
    ap.add_argument("--panel-gap", type=float, default=42)
    ap.add_argument("--panel-height", type=float, default=380)
    ap.add_argument("--seed", type=int, default=7, help="Jitter RNG seed, for reproducible figures")
    return ap.parse_args()


def percentile(sorted_vals, p):
    n = len(sorted_vals)
    idx = p * (n - 1)
    lo = int(idx)
    hi = min(lo + 1, n - 1)
    frac = idx - lo
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * frac


def box_stats(vals):
    s = sorted(vals)
    q1 = percentile(s, 0.25)
    med = percentile(s, 0.5)
    q3 = percentile(s, 0.75)
    iqr = q3 - q1
    lo_fence, hi_fence = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    inliers = [v for v in s if lo_fence <= v <= hi_fence]
    lo_whisk = min(inliers) if inliers else s[0]
    hi_whisk = max(inliers) if inliers else s[-1]
    mean = sum(vals) / len(vals)
    return {'q1': q1, 'median': med, 'q3': q3, 'lo_whisk': lo_whisk, 'hi_whisk': hi_whisk,
            'mean': mean, 'min': s[0], 'max': s[-1]}


def nice_ticks(vmin, vmax, count=5):
    """Evenly-spaced whole-number-friendly ticks (1/2/5 x 10^k steps)."""
    span = vmax - vmin
    if span <= 0:
        return [round(vmin)]
    raw_step = span / count
    mag = 10 ** math.floor(math.log10(raw_step))
    err = raw_step / mag
    if err >= 7.5:
        step = 10 * mag
    elif err >= 3.5:
        step = 5 * mag
    elif err >= 1.5:
        step = 2 * mag
    else:
        step = mag
    step = max(step, 1)  # keep whole-number labels
    start = math.ceil(vmin / step) * step
    ticks = []
    v = start
    while v <= vmax + 1e-9:
        ticks.append(round(v, 6))
        v += step
    return ticks


def read_column_values(path, column):
    import csv
    vals = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if column not in (reader.fieldnames or []):
            raise ValueError(f"Column '{column}' not found in {path}. "
                              f"Available columns: {reader.fieldnames}")
        for row in reader:
            if row[column] != "":
                vals.append(float(row[column]))
    if not vals:
        raise ValueError(f"Column '{column}' in {path} had no non-empty values to plot")
    return vals


def main():
    args = parse_args()
    labels = args.labels or args.columns
    if len(labels) != len(args.columns):
        raise ValueError("--labels must match --columns in count if provided")

    panels = []
    for i, (col, label) in enumerate(zip(args.columns, labels)):
        vals = read_column_values(args.input_csv, col)
        panels.append((col, label, PALETTE[i % len(PALETTE)], vals))
        print(f"  {label}: n={len(vals)} min={min(vals):.2f} max={max(vals):.2f}")

    n = len(vals) if panels else 0
    n_panels = len(panels)
    top_pad, bottom_pad, left_margin = 46, 40, 70
    box_half_w, jitter_w = 26, 20
    total_w = left_margin + n_panels * args.panel_width + (n_panels - 1) * args.panel_gap + 30
    total_h = top_pad + args.panel_height + bottom_pad + 20

    svg = [f'<svg viewBox="0 0 {total_w} {total_h}" xmlns="http://www.w3.org/2000/svg" '
           f'font-family="system-ui, -apple-system, Segoe UI, sans-serif">',
           f'<rect x="0" y="0" width="{total_w}" height="{total_h}" fill="{SURFACE}"/>',
           f'<text x="{-(top_pad + args.panel_height/2)}" y="24" text-anchor="middle" transform="rotate(-90)" '
           f'font-size="15" font-weight="600" fill="{INK_SEC}">{args.title} (n={n})</text>']

    random.seed(args.seed)

    for i, (col, label, color, vals) in enumerate(panels):
        st = box_stats(vals)
        px = left_margin + i * (args.panel_width + args.panel_gap)
        cx = px + args.panel_width / 2
        span = st['max'] - st['min']
        pad = span * 0.12 if span > 0 else max(st['max'] * 0.1, 1)
        y_min, y_max = st['min'] - pad, st['max'] + pad

        def yscale(v, y_min=y_min, y_max=y_max):
            return top_pad + args.panel_height * (1 - (v - y_min) / (y_max - y_min))

        svg.append(f'<text x="{cx}" y="{top_pad - 18}" text-anchor="middle" font-size="14" '
                    f'font-weight="600" fill="{INK}">{label}</text>')

        for label_val in nice_ticks(y_min, y_max, count=5):
            ty = yscale(label_val)
            svg.append(f'<line x1="{px}" x2="{px + args.panel_width}" y1="{ty:.1f}" y2="{ty:.1f}" '
                        f'stroke="{GRID}" stroke-width="1"/>')
            svg.append(f'<text x="{px - 8}" y="{ty + 3:.1f}" text-anchor="end" font-size="10" '
                        f'fill="{MUTED}">{round(label_val)}</text>')

        svg.append(f'<line x1="{px}" x2="{px}" y1="{top_pad:.1f}" y2="{top_pad+args.panel_height:.1f}" '
                    f'stroke="{BASELINE}" stroke-width="1"/>')

        y_q1, y_q3, y_med = yscale(st['q1']), yscale(st['q3']), yscale(st['median'])
        y_lo, y_hi, y_mean = yscale(st['lo_whisk']), yscale(st['hi_whisk']), yscale(st['mean'])

        svg.append(f'<line x1="{cx}" x2="{cx}" y1="{y_hi:.1f}" y2="{y_q3:.1f}" stroke="{INK_SEC}" stroke-width="1.5"/>')
        svg.append(f'<line x1="{cx}" x2="{cx}" y1="{y_q1:.1f}" y2="{y_lo:.1f}" stroke="{INK_SEC}" stroke-width="1.5"/>')
        svg.append(f'<line x1="{cx-12}" x2="{cx+12}" y1="{y_hi:.1f}" y2="{y_hi:.1f}" stroke="{INK_SEC}" stroke-width="1.5"/>')
        svg.append(f'<line x1="{cx-12}" x2="{cx+12}" y1="{y_lo:.1f}" y2="{y_lo:.1f}" stroke="{INK_SEC}" stroke-width="1.5"/>')

        box_top, box_h = min(y_q1, y_q3), max(abs(y_q1 - y_q3), 1.5)
        svg.append(f'<rect x="{cx-box_half_w}" y="{box_top:.1f}" width="{box_half_w*2}" height="{box_h:.1f}" '
                    f'fill="{color}" fill-opacity="0.14" stroke="{INK_SEC}" stroke-width="1.5"/>')
        svg.append(f'<line x1="{cx-box_half_w}" x2="{cx+box_half_w}" y1="{y_med:.1f}" y2="{y_med:.1f}" '
                    f'stroke="{INK}" stroke-width="2"/>')

        ms = 5
        svg.append(f'<line x1="{cx-ms}" x2="{cx+ms}" y1="{y_mean-ms:.1f}" y2="{y_mean+ms:.1f}" stroke="{INK}" stroke-width="1.5"/>')
        svg.append(f'<line x1="{cx-ms}" x2="{cx+ms}" y1="{y_mean+ms:.1f}" y2="{y_mean-ms:.1f}" stroke="{INK}" stroke-width="1.5"/>')

        for v in vals:
            jx = cx + (random.random() - 0.5) * jitter_w
            jy = yscale(v)
            svg.append(f'<circle cx="{jx:.1f}" cy="{jy:.1f}" r="3.2" fill="{color}" fill-opacity="0.75" '
                        f'stroke="{SURFACE}" stroke-width="0.6"/>')

    svg.append('</svg>')
    with open(args.output_svg, 'w') as f:
        f.write('\n'.join(svg))
    print(f"Wrote {args.output_svg}")


if __name__ == "__main__":
    main()
