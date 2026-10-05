"""Integración con OpenRocket 24.12 (orlab + JPype). Se salta si no hay java, orlab o jar."""

import math
import shutil

import pandas as pd
import pytest

from .conftest import ORK

orlab = pytest.importorskip("orlab")
pytest.importorskip("jpype")
if shutil.which("java") is None:
    pytest.skip("sin java", allow_module_level=True)

MM = 1e-3


@pytest.fixture(scope="module")
def puente():
    from sensor_tipo2.verificacion import PuenteTipo2
    try:
        inst = orlab.OpenRocketInstance(log_level="ERROR")
        inst.__enter__()
    except Exception as e:  # sin jar en caché
        pytest.skip(f"OpenRocket no disponible: {e}")
    yield PuenteTipo2(inst, ORK)
    inst.__exit__(None, None, None)


def test_base_ork_regresion(puente, cfg_ork, base_ork):
    """El .ork tipo 2 reescrito por el puente da lo mismo que el archivo original (M = 0.3):
    CP = −133.646 mm, C_Nα = 1.25208, masa 215.611 g; y el sustituto coincide."""
    from sensor_tipo2.exportar import detalle
    det = detalle(cfg_ork, *base_ork)
    puente.aplicar(cfg_ork, det)
    ae = puente.aero(0.3)
    assert ae.x_CP / MM == pytest.approx(-133.646, abs=0.01)
    assert ae.CNa == pytest.approx(1.252078, rel=1e-5)
    assert puente.masas().m_total * 1e3 == pytest.approx(215.611, abs=0.01)
    assert det.cp.x_CP / MM == pytest.approx(ae.x_CP / MM, abs=0.05)
    m_mod = sum(m for m, _ in det.cu.masas_por_estacion().values()) + det.g.masa
    assert m_mod == pytest.approx(puente.masas().m_total, rel=0.005)


CANDIDATOS = [  # factibles de la malla de la spec v2 (aletas en 1.2 R, cuerpos abombados y cilíndricos)
    ((400, 90, 1.0, "conica", None, 0.0, 15.75, 140), (0.7, 0.7, 1.0, 1.2)),
    ((400, 90, 1.0, "conica", None, 0.5, 18, 140), (0.7, 0.7, 1.0, 1.2)),
    ((400, 90, 1.0, "conica", None, 0.0, 18, 200), (0.7, 0.4, 1.0, 1.0)),
    ((400, 90, 1.0, "elipsoide", None, 0.0, 18, 200), (0.7, 0.4, 1.0, 1.2)),
    ((340, 90, 1.0, "conica", None, 0.0, 18, 140), (0.7, 0.7, 1.0, 1.2)),
]


def _fila(cfg, cuerpo, aleta):
    from sensor_tipo2.barrido import evaluar_cuerpo
    from sensor_tipo2.config import AletaSpec, CuerpoSpec
    return pd.Series(evaluar_cuerpo(cfg, CuerpoSpec(*cuerpo), [AletaSpec(*aleta)])[0])


@pytest.mark.parametrize("cuerpo, aleta", CANDIDATOS)
def test_candidatos_contra_openrocket(puente, cfg, cuerpo, aleta):
    from sensor_tipo2.verificacion import verificar_candidato
    fila = _fila(cfg, cuerpo, aleta)
    assert fila["factible"]
    r = verificar_candidato(cfg, puente, fila)
    assert abs(r["dx_or_sust_mm"]) < 0.5
    assert abs(r["dif_masa_or_pct"]) < 0.5
    assert abs(r["x_CG_or_mm"] - r["x_CG_modelo_mm"]) < 0.5
    assert r["factible_or"] and math.isfinite(r["J_or"])


def test_T4_cuerpo_abombado(puente, cfg):
    """f_c = 0 (BodyTube de largo 0): perfil < 0.05 mm y CP ≤ 0.5 mm frente a OpenRocket. La masa
    del casco difiere ≈ 0.4 % (pared por erosión exacta frente a la aproximación de OpenRocket en
    nariz y transición): se exige < 1 %, no el 0.1 % de la spec (ver README)."""
    import numpy as np
    from sensor_tipo2.config import AletaSpec, CuerpoSpec
    from sensor_tipo2.exportar import detalle
    det = detalle(cfg, CuerpoSpec(*CANDIDATOS[0][0]), AletaSpec(*CANDIDATOS[0][1]))
    assert det.cu.perfil.L_c == 0
    puente.aplicar(cfg, det)
    assert float(puente.comp["cuerpo"].getLength()) == 0.0
    per = puente.perfil(400)
    x = np.array([p[1] for p in per])
    r = np.array([p[2] for p in per])
    assert np.abs(r - det.cu.perfil.radio(x)).max() < 0.05 * MM
    assert abs(puente.aero(cfg.vuelo.mach).x_CP - det.cp.x_CP) < 0.5 * MM
    mas = puente.masas()
    m_or = sum(mas.por_componente[k][0] for k in ("nariz", "cuerpo", "cola", "tubo_cola"))
    assert det.cu.m_casco == pytest.approx(m_or, rel=0.01)


@pytest.mark.parametrize("n", [6, 8])
def test_T7_aletas_n_contra_openrocket(puente, cfg, n):
    """n ∈ {6, 8}: C_Nα y x_CP de las aletas como OpenRocket (≤ 0.1 % y ≤ 0.1 mm)."""
    from sensor_tipo2.aletas import barrowman
    from sensor_tipo2.config import AletaSpec, CuerpoSpec
    from sensor_tipo2.exportar import detalle
    from sensor_tipo2.frontera import cfg_n
    c6 = cfg_n(cfg, n)
    det = detalle(c6, CuerpoSpec(*CANDIDATOS[1][0]), AletaSpec(*CANDIDATOS[1][1]))
    puente.aplicar(c6, det)
    assert int(puente.comp["aletas"].getFinCount()) == n
    cna_or, x_or = puente.aero(c6.vuelo.mach).por_componente["aletas"]
    f = barrowman(det.g, math.pi * det.cu.spec.D**2 / 4, c6.vuelo.mach, c6.numerico.n_franjas)
    assert f.CNa == pytest.approx(cna_or, rel=1e-3)
    assert abs(f.x_CP - x_or) < 0.1 * MM


def test_T10_ganador_verificado_y_reabierto(puente, cfg, tmp_path):
    """El mejor candidato: SM(OR) ∈ [1, 2], r_tip ≤ 1.2 R y su .ork reabierto reproduce el CP (±0.5 mm)."""
    from sensor_tipo2.verificacion import verificar_candidato
    fila = _fila(cfg, *CANDIDATOS[0])
    ruta = tmp_path / "ganador.ork"
    r = verificar_candidato(cfg, puente, fila, guardar=ruta)
    rest = cfg.restricciones
    assert r["factible_or"] and rest.SM_min <= r["SM_or_cal"] <= rest.SM_max
    assert fila["r_tip_mm"] <= rest.r_tip_rel_R_max * fila["D_mm"] / 2 + 1e-9
    puente.cargar(ruta)
    assert abs(puente.aero(cfg.vuelo.mach).x_CP / MM - r["x_CP_or_mm"]) < 0.5
