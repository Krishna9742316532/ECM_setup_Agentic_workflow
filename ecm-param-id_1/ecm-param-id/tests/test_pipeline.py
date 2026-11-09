"""End-to-end: the pipeline runs, picks the right model and is reproducible."""

import json

from ecm.config import Settings
from ecm.pipeline import run


def test_synthetic_pipeline_picks_1rc(synthetic, tmp_path):
    d, truth = synthetic
    s = run(d, tmp_path, Settings(), truth=truth, make_plots=False)
    assert s["model_selection"]["recommended"] == "NLS 1RC"
    assert s["model_selection"]["ocv_model"] == "linear"   # the truth uses the table
    nls1 = next(m for m in s["methods"] if m["method"] == "NLS 1RC")
    assert nls1["level"] == "OK"
    assert nls1["val"]["rmse_mV"] < 11


def test_real_pipeline_is_reproducible(real, tmp_path):
    a = run(real, tmp_path / "a", Settings(), make_plots=False)
    b = run(real, tmp_path / "b", Settings(), make_plots=False)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    assert (tmp_path / "a" / "params.json").read_text() == (tmp_path / "b" / "params.json").read_text()


def test_real_data_does_not_support_rc_branch(real, tmp_path):
    s = run(real, tmp_path, Settings(), make_plots=False)
    assert s["model_selection"]["recommended"] == "NLS R0-only"
    # The polynomial OCV explains the data: Q ~ 60 Ah and no warnings.
    assert s["model_selection"]["ocv_model"] == "poly10"
    assert abs(s["soc_setup"]["Q_Ah"] - 60) < 0.5
    nls0 = next(m for m in s["methods"] if m["method"] == "NLS R0-only")
    assert nls0["level"] == "OK" and nls0["val"]["rmse_mV"] < 11


def test_linear_ocv_only_gives_validation_warning(real, tmp_path):
    s = run(real, tmp_path, Settings(ocv_models=("linear",)), make_plots=False)
    nls0 = next(m for m in s["methods"] if m["method"] == "NLS R0-only")
    assert nls0["level"] == "WARN" and nls0["val"]["rmse_mV"] > 14
