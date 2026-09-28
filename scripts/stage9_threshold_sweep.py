"""Threshold sweep on the fourteen explicit-solvent benchmark receptors.

Read-only: parses stored frames and rmsd.csv, runs no MD. For each threshold
the CORE is re-derived from the six calibration replicates via the operative
`derive_core_set` / `select_contributors`, then re-applied to every screening
run of that target. Spearman rho is recomputed over all runs and over the
qualifying subset.
"""
from __future__ import annotations
import json, math, statistics, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.app.core.config import JOBS_DIR
from backend.app.free_energy import calibrate as C
from scipy.stats import spearmanr

RUN = Path(__file__).resolve().parent / "run"
TARGETS = ["cmet","tyk2","ptp1b","thrombin","pde2","mcl1","eg5","syk","hif2a","pfkfb3","cdk8","shp2","p38","tnks2"]
THRESHOLDS = [0.60, 0.70, 0.80]


def rho(xs, ys):
    if len(xs) < 3:
        return None, None
    r = spearmanr(xs, ys)
    return (float(r.statistic) if r.statistic == r.statistic else None,
            float(r.pvalue) if r.pvalue == r.pvalue else None)


def main() -> None:
    out = []
    for t in TARGETS:
        cal_meta = json.load(open(RUN / "stage5_calibration" / f"{t}.json"))
        cal_runs = []
        for r in cal_meta["runs"]:
            m = C.measure_run(JOBS_DIR / r["md_job_id"] / "md",
                              label=f"rep{r['replicate']}")
            if m is not None:
                cal_runs.append(m)

        rows = [json.loads(l) for l in
                open(RUN / "stage6_screening" / f"{t}.jsonl")]
        rows = [r for r in rows if r.get("status") == "ok"
                and r.get("dG_pred_kcal") is not None]

        scr = {}
        for r in rows:
            m = C.measure_run(JOBS_DIR / r["md_job_id"] / "md",
                              label=r["ligand_id"])
            if m is not None:
                scr[r["ligand_id"]] = m
        print(f"{t}: {len(cal_runs)} calibration, {len(scr)} screening measured",
              flush=True)

        for th in THRESHOLDS:
            C.CORE_HOLD_FREQUENCY_MIN = th
            good, _ = C.select_contributors(cal_runs)
            core = C.derive_core_set(cal_runs, good)["core"]

            pred_all, exp_all, pred_ok, exp_ok, n_declined = [], [], [], [], 0
            for r in rows:
                m = scr.get(r["ligand_id"])
                if m is None:
                    continue
                pred_all.append(r["dG_pred_kcal"]); exp_all.append(r["dG_exp_kcal"])
                vals = [len(core & f) / len(core)
                        for f in m["frame_sets"][m["win_start"]:]] if core else []
                cr = float(statistics.fmean(vals)) if vals else float("nan")
                if m["landed"] and math.isfinite(cr) and cr >= th:
                    pred_ok.append(r["dG_pred_kcal"]); exp_ok.append(r["dG_exp_kcal"])
                else:
                    n_declined += 1

            ru, pu = rho(pred_all, exp_all)
            rg, pg = rho(pred_ok, exp_ok)
            rec = {"target": t, "threshold": th, "core_size": len(core),
                   "contributors": len(good), "n_total": len(pred_all),
                   "n_qualified": len(pred_ok), "n_declined": n_declined,
                   "rho_ungated": ru, "p_ungated": pu,
                   "rho_gated": rg, "p_gated": pg,
                   "delta_rho": (rg - ru) if (rg is not None and ru is not None) else None}
            out.append(rec)
            print(json.dumps(rec), flush=True)

    C.CORE_HOLD_FREQUENCY_MIN = 0.70
    dest = RUN / "stage9_threshold_sweep.json"
    dest.write_text(json.dumps(out, indent=1))
    print("wrote", dest)


if __name__ == "__main__":
    main()
