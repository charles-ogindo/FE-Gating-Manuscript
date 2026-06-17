"""Publication figures for a completed MD job.

Renders three matplotlib panels (vector SVG + raster PNG each) from the
existing per-job artifacts under jobs/<md_id>/md/:

  1. RMSD panel  — time on x, RMSD (Å) on y. Three series:
                   pose RMSD (solid, primary — the Q6b receptor-frame
                   displacement that gates the verdict), internal RMSD
                   (dashed, diagnostic — pre-Q6b ligand-on-ligand
                   semantics retained as `rmsd_ligand_internal_*`),
                   backbone Cα (light). Horizontal verdict bands:
                   0–2 Å stable, 2–4 Å drifting, >4 Å unstable. Final
                   pose RMSD + verdict annotated.
  2. H-bond panel — time on x, hbond count on y. Annotates the
                   persistence fraction (hbond_persistence_frac) the
                   classifier reports.
  3. Contacts panel — horizontal bar of top per-residue contact
                   frequency. Labels in AUTHOR numbering when the MD
                   summary records `receptor_renumbered=True` (the
                   Q6b PART 2 relabel against the docking receptor).

Every figure carries the footer caveat
"mechanics fixture; pose RMSD = receptor-frame displacement" so a
figure leaving the bundle retains its provenance context.

The module deliberately uses matplotlib's non-interactive `Agg` backend
so it can render in a worker / FastAPI request thread without a display.
"""

from __future__ import annotations

import csv
import io
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib

# Switch to the headless backend BEFORE importing pyplot. Doing it here
# (module import time) avoids the "backend already set" warning if the
# caller has done `import matplotlib.pyplot as plt` earlier.
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402


CAVEAT_FOOTER = (
    "Mechanics fixture; pose RMSD = receptor-frame displacement. "
    "Internal RMSD is the pre-Q6b ligand-on-ligand value, kept as a "
    "diagnostic. Verdict thresholds are MD-convention defaults, not "
    "literature-calibrated."
)

# Verdict band colors (kept low-saturation so the data lines read first).
_BAND_STABLE   = (0.84, 0.94, 0.84, 0.30)  # green
_BAND_DRIFTING = (0.99, 0.92, 0.74, 0.30)  # amber
_BAND_UNSTABLE = (0.99, 0.83, 0.83, 0.30)  # red


# ---------------------------------------------------------------------
# CSV readers — minimal, no pandas; matches the shape rmsd.csv /
# hbonds.csv / contacts.csv emit (see backend/app/md/analyze.py).
# ---------------------------------------------------------------------

def _read_rmsd_csv(path: Path) -> List[Tuple[float, float, float, float]]:
    """Return [(time_ps, bb_ca, lig_pose, lig_internal), ...]. NaN where
    the CSV cell is empty (pre-Q6b artifacts that never had the column
    OR a frame whose shape diverged from the reference)."""
    out: List[Tuple[float, float, float, float]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                t = float(row["time_ps"])
            except (KeyError, ValueError):
                continue
            bb = _try_float(row.get("rmsd_backbone_ca_A"))
            pose = _try_float(row.get("rmsd_ligand_pose_A"))
            internal = _try_float(row.get("rmsd_ligand_internal_A"))
            # Pre-Q6b column name was rmsd_ligand_heavy_A — accept it as
            # the internal series so the figure still renders on legacy
            # CSVs the operator hasn't re-analyzed.
            if internal is None or math.isnan(internal):
                internal = _try_float(row.get("rmsd_ligand_heavy_A"))
            out.append((
                t,
                bb       if bb is not None       else float("nan"),
                pose     if pose is not None     else float("nan"),
                internal if internal is not None else float("nan"),
            ))
    return out


def _read_hbonds_csv(path: Path) -> List[Tuple[float, int]]:
    """Return [(time_ps, hbond_count), ...]. Stops at the "#summary"
    marker the writer appends — those rows aren't time samples."""
    out: List[Tuple[float, int]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None:
            return out
        for row in reader:
            if not row or row[0].startswith("#"):
                break
            try:
                t = float(row[0])
                c = int(row[1]) if len(row) > 1 and row[1] != "" else 0
            except ValueError:
                continue
            out.append((t, c))
    return out


def _read_contacts_csv(path: Path) -> List[Tuple[str, int, str, float]]:
    """Return [(chain, resseq, resname, frac), ...]. Order preserved from
    the file (analyze.py sorts by frac desc when writing)."""
    out: List[Tuple[str, int, str, float]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                out.append((
                    row.get("chain", "") or "",
                    int(row.get("resseq", 0) or 0),
                    row.get("resname", "") or "",
                    float(row.get("frac_frames_in_contact", 0.0) or 0.0),
                ))
            except (TypeError, ValueError):
                continue
    return out


def _try_float(s: Any) -> Optional[float]:
    if s is None or s == "":
        return None
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------
# Figure renderers
# ---------------------------------------------------------------------

def _figure_to_bytes(fig: Figure) -> Tuple[bytes, bytes]:
    """Render the same figure to SVG + PNG. SVG first (vector — for
    paper figures), PNG second (raster — for slide decks / fallback)."""
    svg_buf = io.BytesIO()
    fig.savefig(svg_buf, format="svg", bbox_inches="tight")
    png_buf = io.BytesIO()
    # 150 dpi is comfortable on retina + print; 220 starts to bloat the
    # bundle without much added clarity at the figure sizes used here.
    fig.savefig(png_buf, format="png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    return svg_buf.getvalue(), png_buf.getvalue()


def _footer(fig: Figure, text: str = CAVEAT_FOOTER) -> None:
    """Stamp a small italic caveat at the bottom of the figure. Y=0 in
    figure coordinates is the very bottom — we pad slightly so it sits
    just below the axes label."""
    fig.text(
        0.5, 0.005, text,
        ha="center", va="bottom",
        fontsize=7, style="italic", color="#444",
        wrap=True,
    )


def render_rmsd_figure(
    rmsd_rows: Sequence[Tuple[float, float, float, float]],
    summary: Dict[str, Any],
) -> Tuple[bytes, bytes]:
    """Build the RMSD panel. Returns (svg_bytes, png_bytes)."""
    times = [r[0] for r in rmsd_rows]
    bb    = [r[1] for r in rmsd_rows]
    pose  = [r[2] for r in rmsd_rows]
    internal = [r[3] for r in rmsd_rows]

    verdict = (summary or {}).get("verdict", "unknown")
    metrics = (summary or {}).get("metrics") or {}
    pose_final = metrics.get("rmsd_ligand_pose_final_a")
    pose_max   = metrics.get("rmsd_ligand_pose_max_a")
    internal_final = metrics.get("rmsd_ligand_internal_final_a")
    md_id = (summary or {}).get("md_job_id", "?")
    ligand = (summary or {}).get("ligand", "?")

    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    # Verdict bands — drawn first so the data lines sit on top.
    y_top = max(
        (v for v in pose + internal + bb if v is not None and not math.isnan(v)),
        default=5.0,
    )
    y_top = max(y_top * 1.15, 5.0)
    ax.axhspan(0.0, 2.0, facecolor=_BAND_STABLE)
    ax.axhspan(2.0, 4.0, facecolor=_BAND_DRIFTING)
    ax.axhspan(4.0, y_top, facecolor=_BAND_UNSTABLE)
    # Right-edge band labels — small, low contrast.
    if times:
        x_max = max(times)
        ax.text(x_max, 1.0, "stable",   ha="right", va="center", fontsize=8, color="#2e6b30", alpha=0.75)
        ax.text(x_max, 3.0, "drifting", ha="right", va="center", fontsize=8, color="#8a5a1a", alpha=0.75)
        ax.text(x_max, max(4.4, y_top * 0.85), "unstable", ha="right", va="center", fontsize=8, color="#9c2a2a", alpha=0.75)

    # Three series. Pose first (primary, solid), internal dashed, backbone light.
    ax.plot(times, pose, color="#1f4ea1", linewidth=2.0, label="Ligand pose RMSD (primary)")
    ax.plot(times, internal, color="#1f4ea1", linewidth=1.2, linestyle="--",
            label="Ligand internal RMSD (diagnostic)")
    ax.plot(times, bb, color="#888", linewidth=1.0, alpha=0.85,
            label="Backbone Cα RMSD")

    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_xlim(left=0, right=max(times) if times else 1.0)
    ax.set_ylim(0.0, y_top)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.4)

    # Title + annotation. Verdict + final pose number first; internal as
    # secondary parenthetical.
    title_bits = [f"MD RMSD — {md_id[:8]} / {ligand}"]
    ax.set_title("  ·  ".join(title_bits), fontsize=11)
    annotation_parts = [f"verdict: {verdict}"]
    if isinstance(pose_final, (int, float)):
        annotation_parts.append(f"pose final: {pose_final:.2f} Å")
    if isinstance(pose_max, (int, float)):
        annotation_parts.append(f"pose max: {pose_max:.2f} Å")
    if isinstance(internal_final, (int, float)):
        annotation_parts.append(f"internal final: {internal_final:.2f} Å")
    ax.text(
        0.01, 0.97, "   ".join(annotation_parts),
        transform=ax.transAxes, fontsize=9, va="top",
        bbox=dict(facecolor="white", edgecolor="#ccc", boxstyle="round,pad=0.3", alpha=0.85),
    )

    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    _footer(fig)
    return _figure_to_bytes(fig)


def render_hbonds_figure(
    hbond_rows: Sequence[Tuple[float, int]],
    summary: Dict[str, Any],
) -> Tuple[bytes, bytes]:
    """Build the H-bond count panel. Returns (svg_bytes, png_bytes)."""
    times = [r[0] for r in hbond_rows]
    counts = [r[1] for r in hbond_rows]
    md_id = (summary or {}).get("md_job_id", "?")
    metrics = (summary or {}).get("metrics") or {}
    persistence = metrics.get("hbond_persistence_frac")

    fig, ax = plt.subplots(figsize=(7.5, 3.8))

    ax.plot(times, counts, color="#2c7a3f", linewidth=1.5,
            marker="o", markersize=2.5, label="Protein–ligand H-bonds")
    if counts:
        mean_count = sum(counts) / len(counts)
        ax.axhline(mean_count, color="#2c7a3f", linewidth=0.8, linestyle="--", alpha=0.5,
                   label=f"mean = {mean_count:.1f}")

    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("H-bond count")
    ax.set_xlim(left=0, right=max(times) if times else 1.0)
    ax.set_ylim(bottom=0)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_title(f"Protein–ligand H-bonds — {md_id[:8]}", fontsize=11)

    if isinstance(persistence, (int, float)):
        ax.text(
            0.01, 0.97,
            f"persistence: {persistence * 100:.0f}% of frames have ≥1 H-bond",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(facecolor="white", edgecolor="#ccc", boxstyle="round,pad=0.3", alpha=0.85),
        )

    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _footer(fig)
    return _figure_to_bytes(fig)


def render_contacts_figure(
    contact_rows: Sequence[Tuple[str, int, str, float]],
    summary: Dict[str, Any],
    top_n: int = 15,
) -> Tuple[bytes, bytes]:
    """Build the per-residue contact-frequency panel. Horizontal bars."""
    rows = sorted(contact_rows, key=lambda r: -r[3])[:top_n]
    # y-axis reads top → bottom (highest frac at the top). matplotlib's
    # barh draws bottom-up; reverse the order so the highest value
    # appears at the top of the chart.
    rows = list(reversed(rows))

    labels = [f"{rn}{rs} {ch}" for (ch, rs, rn, _) in rows]
    values = [frac for (_, _, _, frac) in rows]

    md_id = (summary or {}).get("md_job_id", "?")
    renumbered = bool((summary or {}).get("receptor_renumbered"))
    numbering = "author numbering (Q6b relabel)" if renumbered else "MD numbering"

    fig, ax = plt.subplots(figsize=(7.5, max(4.0, 0.32 * len(rows) + 1.0)))

    bars = ax.barh(range(len(values)), values, color="#cc7a00", height=0.7)
    ax.set_yticks(range(len(values)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Fraction of frames in contact")
    ax.set_xlim(0.0, 1.0)
    ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_title(
        f"Top {len(rows)} per-residue contacts — {md_id[:8]}  ({numbering})",
        fontsize=11,
    )
    # Per-bar fraction tag at the right edge.
    for i, b in enumerate(bars):
        ax.text(
            min(b.get_width() + 0.015, 0.98), b.get_y() + b.get_height() / 2,
            f"{values[i]:.2f}",
            va="center", fontsize=8, color="#444",
        )

    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _footer(fig)
    return _figure_to_bytes(fig)


# ---------------------------------------------------------------------
# Convenience: build all 3 figures off a job dir
# ---------------------------------------------------------------------

def render_all_md_figures(
    md_dir: Path,
    summary: Dict[str, Any],
) -> Dict[str, Tuple[bytes, bytes]]:
    """Render the three panels from on-disk artifacts.

    Returns {"rmsd": (svg, png), "hbonds": (svg, png), "contacts": (svg, png)}.

    Missing CSVs degrade gracefully — that panel's bytes pair will be
    omitted from the dict rather than raising. The report bundler treats
    omitted panels as "this view of the data isn't available for this
    job" rather than failing the whole bundle.
    """
    out: Dict[str, Tuple[bytes, bytes]] = {}
    rmsd_path = md_dir / "rmsd.csv"
    hbonds_path = md_dir / "hbonds.csv"
    contacts_path = md_dir / "contacts.csv"

    if rmsd_path.is_file():
        out["rmsd"] = render_rmsd_figure(_read_rmsd_csv(rmsd_path), summary)
    if hbonds_path.is_file():
        out["hbonds"] = render_hbonds_figure(_read_hbonds_csv(hbonds_path), summary)
    if contacts_path.is_file():
        out["contacts"] = render_contacts_figure(_read_contacts_csv(contacts_path), summary)
    return out
