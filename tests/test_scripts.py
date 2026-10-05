"""Scripts de punta a punta en una carpeta temporal: T8 (malla sin factibles), T10 (sin Java),
T11 (planos), T12 (reproducibilidad) y T13 (ningún .ork)."""

import importlib.util
import os
import re
import shutil
import subprocess
import sys

import pandas as pd
import pytest
import yaml

from .conftest import ORK, RAIZ, malla_chica, raw_opt

SCRIPTS_SIN_JVM = ("01_optimizar_malla.py", "03_frontera_factibilidad.py", "04_planos.py")


def _script(nombre):
    spec = importlib.util.spec_from_file_location(nombre.replace(".py", ""), RAIZ / "scripts" / nombre)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _escribir(tmp_path, raw, nombre):
    raw["ejecucion"]["procesos"] = 1
    raw["salida"] = {"dir": str(tmp_path / "data"), "dir_figuras": str(tmp_path / "figs"),
                     "dir_planos": str(tmp_path / "planos"), "N_planos": 2, "dpi_planos": 60}
    raw["openrocket"]["ork_base"] = str(ORK)
    ruta = tmp_path / "config" / nombre
    ruta.parent.mkdir(parents=True, exist_ok=True)
    ruta.write_text(yaml.safe_dump(raw, allow_unicode=True), encoding="utf-8")
    return ruta


def _config_sin_factibles(tmp_path):
    """Cuerpo cilíndrico con tubo de cola corto y aletas en el calibre: nada llega a SM = 1."""
    raw = raw_opt()
    raw["malla"] = {
        "L_total_mm": [400], "D_mm": [60, 70], "L_n_rel_D": [1.0], "f_cil": [0.75],
        "cola_forma": [{"forma": "conica", "parametro": None}], "k": [0.3, 0.4], "L_tc_mm": [60],
        "mu_cr": [1.0], "gamma_ct": [0.7], "sigma_flecha": [1.0], "r_tip_rel_R": [1.0],
    }
    return _escribir(tmp_path, raw, "t8.yaml")


def _config_chica(tmp_path):
    raw = malla_chica(raw_opt())
    raw["ejecucion"]["refinamiento"] = {"activar": False}
    return _escribir(tmp_path, raw, "chica.yaml")


def _sin_ork(tmp_path):
    """T13: ningún .ork en ninguna carpeta de salida."""
    assert not list(tmp_path.rglob("*.ork"))


def test_T8_malla_sin_factibles(tmp_path):
    ruta = _config_sin_factibles(tmp_path)
    assert _script("01_optimizar_malla.py").main(["--config", str(ruta), "--sin-figuras"]) == 0
    rk = pd.read_csv(tmp_path / "data" / "ranking.csv")
    assert len(rk) > 0 and not rk["factible"].any()
    assert _script("03_frontera_factibilidad.py").main(["--config", str(ruta), "--procesos", "1"]) == 0
    fr = pd.read_csv(tmp_path / "data" / "frontera_factibilidad.csv")
    assert fr.loc[~fr["solo_diagnostico"], "n_factibles"].sum() == 0
    cf = pd.read_csv(tmp_path / "data" / "casi_factibles.csv")
    assert len(cf) > 0 and cf["deficit_SM_cal"].is_monotonic_increasing
    assert {"restriccion_que_falla", "SM_max_alcanzable_cal", "m_en_SM_max_g"} <= set(cf.columns)
    assert _script("04_planos.py").main(["--config", str(ruta)]) == 0
    png = sorted((tmp_path / "planos").glob("casi*.png"))
    assert 1 <= len(png) <= 3
    assert all(p.with_name(p.name[:-4] + ".pdf").exists() for p in png)
    assert not list((tmp_path / "planos").glob("rank*"))
    # el plano está rotulado INFACTIBLE
    from sensor_tipo2.barrido import evaluar_cuerpo, specs_de_fila
    from sensor_tipo2.config import cargar
    from sensor_tipo2.planos import plano
    cfg = cargar(ruta)
    c, a = specs_de_fila(cf.iloc[0])
    fila = pd.Series(evaluar_cuerpo(cfg, c, [a])[0])
    assert fila["cand_id"] == cf["cand_id"].iloc[0] and not fila["factible"]
    assert plano(cfg, fila, tmp_path / "x" / "casi", falla="SM", dpi=40).infactible
    _sin_ork(tmp_path)


def test_T10_sin_java_ni_orlab(tmp_path):
    """T10: los scripts 01, 03 y 04 corren completos sin Java ni orlab, y ningún módulo salvo
    verificacion.py importa JPype u orlab."""
    patron = re.compile(r"^\s*(import|from)\s+(jpype|orlab)\b", re.M)
    culpables = [p.relative_to(RAIZ).as_posix() for p in list((RAIZ / "src").rglob("*.py")) + list((RAIZ / "scripts").glob("*.py"))
                 if patron.search(p.read_text(encoding="utf-8")) and p.name != "verificacion.py"]
    assert culpables == []
    falsos = tmp_path / "sin_jvm"
    for mod in ("jpype", "orlab"):
        (falsos / mod).mkdir(parents=True)
        (falsos / mod / "__init__.py").write_text(f"raise ImportError('{mod} bloqueado por el test T10')\n")
    env = {**os.environ, "PYTHONPATH": os.pathsep.join([str(falsos), str(RAIZ / "src")]), "PATH": "/usr/bin:/bin",
           "JAVA_HOME": ""}
    ruta = _config_chica(tmp_path)
    for s in SCRIPTS_SIN_JVM:
        args = ["--procesos", "1"] if s.startswith("03") else []
        r = subprocess.run([sys.executable, str(RAIZ / "scripts" / s), "--config", str(ruta), *args],
                           capture_output=True, text=True, timeout=900, env=env)
        assert r.returncode == 0, f"{s}: {r.stderr[-2000:]}"
    assert list((tmp_path / "planos").glob("rank01_*.pdf"))
    assert (tmp_path / "planos" / "comparativo_top.png").exists()
    _sin_ork(tmp_path)


def test_T11_planos_y_cotas(tmp_path):
    """T11: cada plano existe en PNG y PDF y las cotas dibujadas coinciden con el ranking (±0.05 mm)."""
    from sensor_tipo2.barrido import optimizar
    from sensor_tipo2.config import cargar
    from sensor_tipo2.planos import COL_RANKING, plano
    cfg = cargar(malla_chica(raw_opt()))
    rk = optimizar(cfg, refinar=False, procesos=1)
    fac = rk[rk["factible"]]
    # el mejor (cuerpo abombado) y uno con cuerpo cilíndrico y tapón trasero
    filas = [fac.iloc[0]] + [f for _, f in fac[fac["f_cil"] > 0].head(1).iterrows()]
    assert len(filas) == 2
    for f in filas:
        pl = plano(cfg, f, tmp_path / f"p_{f['puesto']:.0f}", puesto=f["puesto"], dpi=40)
        assert all(r.exists() and r.stat().st_size > 0 for r in pl.rutas)
        assert {r.suffix for r in pl.rutas} == {".png", ".pdf"}
        assert set(COL_RANKING) <= set(pl.cotas)
        for k, col in COL_RANKING.items():
            v, texto = pl.cotas[k]
            esperado = 0.0 if pd.isna(f[col]) else float(f[col])
            assert abs(float(texto) - esperado) <= 0.05 + 1e-9, (k, texto, esperado)
    _sin_ork(tmp_path)


def test_T12_misma_semilla_mismo_ranking():
    from sensor_tipo2.barrido import optimizar
    from sensor_tipo2.config import cargar
    cfg = cargar(malla_chica(raw_opt()))
    a, b = optimizar(cfg, procesos=1), optimizar(cfg, procesos=2)
    pd.testing.assert_frame_equal(a, b)


@pytest.mark.skipif(shutil.which("java") is None or importlib.util.find_spec("orlab") is None,
                    reason="sin java u orlab")
def test_script_02_valida_sin_ork(tmp_path):
    """01 → 02 → 04: validacion_or.csv, ganador validado en el ranking, cajetín con la validación
    y ningún .ork (T13). Con la malla sin factibles, 02 termina sin ganador."""
    ruta = _config_chica(tmp_path)
    assert _script("01_optimizar_malla.py").main(["--config", str(ruta), "--sin-figuras"]) == 0
    # en otro proceso: JPype no puede volver a arrancar la JVM dentro de la sesión de pytest
    r = subprocess.run([sys.executable, str(RAIZ / "scripts" / "02_validar_ganadores.py"), "--config", str(ruta)],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stderr[-2000:]
    v = pd.read_csv(tmp_path / "data" / "validacion_or.csv")
    assert len(v) == 5 and v["validado"].any() and (v["dx_cp_mm"].abs() <= 2.0).all()
    rk = pd.read_csv(tmp_path / "data" / "ranking.csv")
    assert rk["ganador"].sum() == 1 and bool(rk.loc[rk["ganador"], "validado_or"].iloc[0])
    assert _script("04_planos.py").main(["--config", str(ruta)]) == 0
    _sin_ork(tmp_path)
    t8 = tmp_path / "t8"
    ruta8 = _config_sin_factibles(t8)
    assert _script("01_optimizar_malla.py").main(["--config", str(ruta8), "--sin-figuras"]) == 0
    r = subprocess.run([sys.executable, str(RAIZ / "scripts" / "02_validar_ganadores.py"), "--config", str(ruta8)],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0 and "Sin ganador" in r.stdout, r.stderr[-2000:]
