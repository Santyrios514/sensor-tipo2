"""T8 (malla sin factibles) y T9 (reproducibilidad) de los scripts, en una carpeta temporal."""

import importlib.util
import json
import shutil
import subprocess
import sys

import pandas as pd
import pytest
import yaml

from .conftest import ORK, RAIZ, raw_opt


def _script(nombre):
    spec = importlib.util.spec_from_file_location(nombre.replace(".py", ""), RAIZ / "scripts" / nombre)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _config_sin_factibles(tmp_path):
    """Cuerpo cilíndrico con tubo de cola corto y aletas en el calibre: nada llega a SM = 1."""
    raw = raw_opt()
    raw["malla"] = {
        "L_total_mm": [400], "D_mm": [60, 70], "L_n_rel_D": [1.0], "f_cil": [0.75],
        "cola_forma": [{"forma": "conica", "parametro": None}], "k": [0.3, 0.4], "L_tc_mm": [60],
        "mu_cr": [1.0], "gamma_ct": [0.7], "sigma_flecha": [1.0], "r_tip_rel_R": [1.0],
    }
    raw["ejecucion"]["procesos"] = 1
    raw["salida"] = {"dir": str(tmp_path / "data"), "dir_figuras": str(tmp_path / "figs"), "exportar_ork_top": 2}
    raw["openrocket"]["ork_base"] = str(ORK)
    ruta = tmp_path / "config" / "t8.yaml"
    ruta.parent.mkdir(parents=True)
    ruta.write_text(yaml.safe_dump(raw, allow_unicode=True), encoding="utf-8")
    return ruta


def test_T8_malla_sin_factibles(tmp_path):
    ruta = _config_sin_factibles(tmp_path)
    assert _script("01_optimizar_malla.py").main(["--config", str(ruta)]) == 0
    rk = pd.read_csv(tmp_path / "data" / "ranking.csv")
    assert len(rk) > 0 and not rk["factible"].any()
    assert _script("03_frontera_factibilidad.py").main(["--config", str(ruta), "--procesos", "1"]) == 0
    fr = pd.read_csv(tmp_path / "data" / "frontera_factibilidad.csv")
    assert fr.loc[~fr["solo_diagnostico"], "n_factibles"].sum() == 0
    cf = pd.read_csv(tmp_path / "data" / "casi_factibles.csv")
    assert len(cf) > 0 and cf["deficit_SM_cal"].is_monotonic_increasing
    assert {"restriccion_que_falla", "SM_max_alcanzable_cal", "m_en_SM_max_g"} <= set(cf.columns)


@pytest.mark.skipif(shutil.which("java") is None or importlib.util.find_spec("orlab") is None,
                    reason="sin java u orlab")
def test_T8_verificacion_sin_ganador(tmp_path):
    ruta = _config_sin_factibles(tmp_path)
    assert _script("01_optimizar_malla.py").main(["--config", str(ruta), "--sin-figuras"]) == 0
    # en otro proceso: JPype no puede volver a arrancar la JVM dentro de la sesión de pytest
    r = subprocess.run([sys.executable, str(RAIZ / "scripts" / "02_verificar_openrocket.py"), "--config", str(ruta),
                        "--sin-figuras"], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "Sin ganador" in r.stdout
    assert json.loads((tmp_path / "data" / "calibracion.json").read_text())["ganador"] is None


def test_T9_misma_semilla_mismo_ranking(tmp_path):
    from sensor_tipo2.barrido import optimizar
    from sensor_tipo2.config import cargar
    from .conftest import malla_chica
    cfg = cargar(malla_chica(raw_opt()))
    a, b = optimizar(cfg, procesos=1), optimizar(cfg, procesos=1)
    pd.testing.assert_frame_equal(a, b)
