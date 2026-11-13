#!/usr/bin/env python3
"""Pre-specified PRJNA310127 exclusion, first completed late on 2026-11-12.

Run: python3 scripts/meta_excl310127.py

Reuse the original pooling statements, including their nll function, without
importing the scripts: their top-level input/output code rewrites old results.
Only the three meta_gradient_excl310127 outputs are written by this script.
"""

import ast
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize_scalar
from scipy.stats import t as tdist


ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-11-12"
TARGET = "China PRJNA310127 16S (SILVA)"
VARIANTS = (
    ("withF2_noSpectrum", "F2 included; Spectrum excluded (headline)", "meta_withF2.py"),
    ("withF2", "F2 included; all datasets", "meta_withF2.py"),
    ("4cohorts_noSpectrum", "Original four-dataset synthesis; Spectrum excluded", "meta4.py"),
    ("4cohorts", "Original four-dataset synthesis; all datasets", "meta4.py"),
)
METRICS = ("k", "total_pairs", "r", "lo", "hi", "P", "I2", "tau2", "Q")
TIMING = (
    "The PRJNA310127 exclusion sensitivity was pre-specified in "
    "plans/analysis_plan_v0.6.md, section 2, row 7 (and manuscript/"
    "protocol_v0.5_addendum.md, section B), but was first computed late on "
    "12 November 2026 (Asia/Shanghai), after the input outcomes were known. "
    "It is a late-completed pre-specified sensitivity, not a new blinded test."
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_pooling_code(path):
    """Compile the original z/v assignment through res, with no AST edits."""
    source = path.read_text()
    body = ast.parse(source, filename=str(path)).body
    start = next(
        i for i, node in enumerate(body)
        if isinstance(node, ast.Assign)
        and ast.unparse(node.targets[0]) == "(z, v)"
    )
    stop = next(
        i for i in range(start, len(body))
        if isinstance(body[i], ast.Assign)
        and ast.unparse(body[i].targets[0]) == "res"
    )
    block = body[start:stop + 1]
    assert any(isinstance(n, ast.FunctionDef) and n.name == "nll" for n in block)
    assert all(isinstance(n, (ast.Assign, ast.FunctionDef)) for n in block)
    code = compile(ast.Module(body=block, type_ignores=[]), str(path), "exec")
    return code, [body[start].lineno, body[stop].end_lineno]


def pool(frame, code):
    namespace = {
        "M": frame.copy(), "np": np, "tdist": tdist,
        "minimize_scalar": minimize_scalar,
    }
    exec(code, namespace)
    # Convert numpy scalars, retaining the original result keys and values.
    return json.loads(json.dumps(namespace["res"], default=float, allow_nan=False))


def validate_inputs(frame):
    assert frame.study.is_unique, "Duplicate dataset labels"
    assert int((frame.study == TARGET).sum()) == 1, "Expected exactly one PRJNA310127 row"
    assert np.isfinite(frame[["pairs", "nonzero", "r", "z", "se"]]).all().all()
    assert (frame.nonzero > 3).all() and (frame.nonzero <= frame.pairs).all()
    assert (frame.pairs % 1 == 0).all() and (frame.nonzero % 1 == 0).all()
    assert (frame.r.abs() < 1).all() and (frame.se > 0).all()
    np.testing.assert_allclose(frame.z, np.arctanh(frame.r), rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(frame.se, 1 / np.sqrt(frame.nonzero - 3), rtol=1e-12)


def result_cell(result):
    return (
        f"{result['r']:.4f} [{result['lo']:.4f}, {result['hi']:.4f}]; "
        f"P={result['P']:.6g}; I²={100 * result['I2']:.2f}%; "
        f"τ²={result['tau2']:.8g}; k={result['k']}; pairs={result['total_pairs']}"
    )


def main():
    source_paths = [
        ROOT / "scripts/meta_withF2.py", ROOT / "scripts/meta4.py",
        ROOT / "plans/analysis_plan_v0.6.md",
        ROOT / "manuscript/protocol_v0.5_addendum.md",
    ]
    for variant, _, _ in VARIANTS:
        source_paths.extend(ROOT / f"results/meta_gradient_{variant}.{ext}" for ext in ("csv", "json"))
    original_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in source_paths}
    analyses = {}
    rows = []
    for variant, label, script in VARIANTS:
        csv_path = ROOT / f"results/meta_gradient_{variant}.csv"
        json_path = csv_path.with_suffix(".json")
        frame = pd.read_csv(csv_path, float_precision="round_trip")
        validate_inputs(frame)
        code, lines = load_pooling_code(ROOT / "scripts" / script)
        with_target = pool(frame, code)
        stored = json.loads(json_path.read_text())
        # Check the unexcluded synthesis against the previously saved result.
        errors = {}
        for key in METRICS:
            np.testing.assert_allclose(with_target[key], stored[key], rtol=1e-9, atol=1e-12,
                                       err_msg=f"Baseline reproduction: {variant}.{key}")
            errors[key] = float(abs(with_target[key] - stored[key]))
        retained = frame.loc[frame.study != TARGET].copy()
        without_target = pool(retained, code)
        removed_pairs = int(frame.loc[frame.study == TARGET, "pairs"].sum())
        assert without_target["k"] == with_target["k"] - 1
        assert without_target["total_pairs"] == with_target["total_pairs"] - removed_pairs
        analyses[variant] = {
            "label": label,
            "input_csv": str(csv_path.relative_to(ROOT)),
            "baseline_json": str(json_path.relative_to(ROOT)),
            "pooling_source": f"scripts/{script}",
            "pooling_source_lines": lines,
            "baseline_check": {"passed": True, "rtol": 1e-9, "atol": 1e-12, "absolute_differences": errors},
            "excluded_study": TARGET,
            "excluded_pairs": removed_pairs,
            "input_studies": frame.to_dict(orient="records"),
            "retained_studies": retained.study.tolist(),
            "with_PRJNA310127": with_target,
            "without_PRJNA310127": without_target,
            "direction_changed": bool((with_target["r"] > 0) != (without_target["r"] > 0)),
            "nominal_P05_status_changed": bool((with_target["P"] < .05) != (without_target["P"] < .05)),
        }
        for included, result in ((True, with_target), (False, without_target)):
            rows.append({
                "analysis": variant, "PRJNA310127_included": included,
                "analysis_date": DATE, "input_csv": str(csv_path.relative_to(ROOT)),
                "pooling_source": f"scripts/{script}",
                **{key: result[key] for key in METRICS},
                "I2_percent": 100 * result["I2"], "method": result["method"],
            })

    moderator_note = (
        "The plan identifies PRJNA310127 as gastric cardia cancer with distant normal "
        "sampling (manuscript/protocol_v0.5_addendum.md:60-63; "
        "plans/analysis_plan_v0.6.md:50,64). A joint cohort-level indicator is codable, "
        "but its positive category contains only PRJNA310127 and is inseparable from "
        "other features of that cohort; it cannot distinguish cardia location from "
        "distant-normal sampling. No moderator model was fitted. The exclusion "
        "implements the pre-specified alternative."
    )
    unchanged = all(
        not a["direction_changed"] and not a["nominal_P05_status_changed"]
        for a in analyses.values()
    )
    conclusion = (
        "All four exclusions retain positive pooled effects, 95% CIs above zero, "
        "and nominal P < 0.05. The normal > tumour gradient conclusion is unchanged. "
        "These sensitivities do not change the frozen confirmatory decision rules."
        if unchanged and all(a["without_PRJNA310127"]["lo"] > 0 for a in analyses.values())
        else "Inspect the paired results: a direction, interval, or nominal P < 0.05 conclusion changed."
    )
    result = {
        "analysis_date": DATE, "timezone": "Asia/Shanghai", "timing_disclosure": TIMING,
        "script": "scripts/meta_excl310127.py", "script_sha256": sha256(Path(__file__)),
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "scipy": scipy.__version__},
        "method": "Unchanged original REML + Hartung-Knapp pooling code, including nll",
        "I2_units": "proportion; I2_percent in CSV and percentages in the summary",
        "tau2_scale": "Fisher z",
        "pair_count_note": "Sum of dataset pair contributions; Spectrum RNA and TCMA WGS overlap and are not independent cohorts.",
        "source_sha256": original_hashes, "analyses": analyses,
        "moderator_feasibility": moderator_note, "conclusion": conclusion,
    }
    assert all(sha256(ROOT / p) == value for p, value in original_hashes.items()), "An original source changed"
    stem = ROOT / "results/meta_gradient_excl310127"
    stem.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    pd.DataFrame(rows).to_csv(stem.with_suffix(".csv"), index=False)
    summary = [
        "# PRJNA310127 exclusion sensitivity", "", TIMING, "",
        "## Method and provenance", "",
        "Each existing per-study CSV was pooled with and without the exact row "
        f"`{TARGET}` (77 pairs). No input effects, standard errors, or thresholds were changed. "
        "The script executes the original pooling statements and `nll` function from "
        "`scripts/meta_withF2.py` or `scripts/meta4.py` without running their input/output code. "
        "All four unexcluded results reproduce the corresponding stored JSON estimates "
        "within relative tolerance 1e-9 and absolute tolerance 1e-12. Original source "
        "SHA-256 hashes, software versions, and per-metric reproduction differences are in the new JSON.", "",
        "REML uses the original bounded (0, 5) optimizer and its existing near-zero τ² rule; "
        "Hartung–Knapp uses the unchanged residual multiplier (floor 1e-12, not 1) and "
        "t distribution with k − 1 degrees of freedom. Pooling uses stored Fisher z "
        "and standard errors, with r-scale 95% CIs. I² is shown as a percentage; "
        "τ² is on the Fisher-z scale. Study-level P values are not pooled.", "",
        "## Results", "",
        "Cells report pooled r [95% CI], two-sided P, I², τ², k and total pairs.", "",
        "| Synthesis | With PRJNA310127 | Without PRJNA310127 |", "|---|---|---|",
    ]
    for a in analyses.values():
        summary.append(f"| {a['label']} | {result_cell(a['with_PRJNA310127'])} | {result_cell(a['without_PRJNA310127'])} |")
    summary.extend([
        "", "k counts contributing datasets. Total pairs are pair contributions, not "
        "necessarily distinct patients across modalities. The all-dataset versions "
        "include overlapping TCGA Spectrum RNA and TCMA WGS; the headline excludes Spectrum. "
        "The `4cohorts` labels identify the original input family, not k after exclusions. "
        "The no-F2, Spectrum-excluded sensitivity contains only two datasets after exclusion.", "",
        "## Interpretation", "", conclusion, "",
        "## Cardia/distant-normal moderator", "", moderator_note, "",
        "## Outputs and reproduction", "",
        "Run `python3 scripts/meta_excl310127.py` from the project directory. The script "
        "writes only `results/meta_gradient_excl310127.json`, "
        "`results/meta_gradient_excl310127.csv`, and this summary. The JSON and CSV "
        "retain unrounded estimates for both inclusion conditions. Existing scripts "
        "and results are read only.", "",
    ])
    Path(str(stem) + "_SUMMARY.md").write_text("\n".join(summary))
    print(pd.DataFrame(rows)[["analysis", "PRJNA310127_included", *METRICS]].to_string(index=False))
    print(conclusion)


if __name__ == "__main__":
    main()
