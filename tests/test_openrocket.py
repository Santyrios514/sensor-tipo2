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
    """Criterio de validación del §7: |Δx_CP| ≤ tol_cp_mm, |Δm| ≤ 0.5 % y SM_OR ∈ [1, 2]."""
    from sensor_tipo2.verificacion import validar_candidato
    fila = _fila(cfg, cuerpo, aleta)
    assert fila["factible"]
    r = validar_candidato(cfg, puente, fila, tol_cp=2.0 * MM, tol_masa_pct=0.5)
    assert r["validado"], r["motivo"]
    assert abs(r["dx_cp_mm"]) < 0.5
    assert abs(r["dif_masa_or_pct"]) < 0.5
    assert abs(r["dx_cg_mm"]) < 0.5


def _masa_esquema_or(forma, R, L, param, t, rho, N=128):
    """Réplica de SymmetricComponent.calculateProperties (OpenRocket 24.12): N troncos de cono con la
    pared medida en vertical sobre la secante, r_i = r − t·hyp/l (≥ 0)."""
    import numpy as np
    from sensor_tipo2.perfiles import forma_or
    x = np.linspace(0.0, L, N + 1)
    r = forma_or(x, forma, R, L, param)
    l = L / N
    h = t * np.hypot(np.diff(r), l) / l
    r1, r2 = r[:-1], r[1:]
    i1, i2 = np.maximum(r1 - h, 0), np.maximum(r2 - h, 0)
    return rho * np.pi / 3 * l * ((r1 * r1 + r1 * r2 + r2 * r2) - (i1 * i1 + i1 * i2 + i2 * i2)).sum()


def test_T4_cuerpo_abombado(puente, cfg):
    """f_c = 0 (BodyTube de largo 0): perfil < 0.05 mm, CP ≤ 0.5 mm y masas del casco < 0.1 % en
    transición y tubo de cola. La nariz difiere ≈ 0.5 %: OpenRocket mide la pared en vertical sobre
    la secante (t/cos θ), que en una nariz roma y curva da menos pared que el espesor normal uniforme
    de la erosión; se comprueba que el valor de OpenRocket es exactamente su esquema replicado sobre
    el mismo perfil (≤ 0.01 %), es decir, que la geometría coincide y la diferencia es el esquema."""
    import numpy as np
    from sensor_tipo2.config import AletaSpec, CuerpoSpec
    from sensor_tipo2.exportar import detalle
    from sensor_tipo2.verificacion import densidad_pared
    from sensor_tipo2.config import espesor_total
    det = detalle(cfg, CuerpoSpec(*CANDIDATOS[0][0]), AletaSpec(*CANDIDATOS[0][1]))
    p = det.cu.perfil
    assert p.L_c == 0
    puente.aplicar(cfg, det)
    assert float(puente.comp["cuerpo"].getLength()) == 0.0
    per = puente.perfil(400)
    x = np.array([q[1] for q in per])
    r = np.array([q[2] for q in per])
    assert np.abs(r - p.radio(x)).max() < 0.05 * MM
    assert abs(puente.aero(cfg.vuelo.mach).x_CP - det.cp.x_CP) < 0.5 * MM
    m_or = {k: v[0] for k, v in puente.masas().por_componente.items()}
    m_mod = {k: v[0] for k, v in det.cu.masas_por_estacion().items()}
    for k in ("cola", "tubo_cola"):
        assert m_mod[k] == pytest.approx(m_or[k], rel=1e-3), k
    capas = cfg.pared["nariz"]
    m_rep = _masa_esquema_or(p.forma_n, p.R, p.Ln, p.param_n, espesor_total(capas),
                             densidad_pared(det.cu.cav, capas, "nariz"))
    assert m_or["nariz"] == pytest.approx(m_rep, rel=1e-4)
    assert m_mod["nariz"] == pytest.approx(m_or["nariz"], rel=0.01)


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


def test_validar_ganadores_sin_ork(puente, cfg, tmp_path):
    """§7: el mejor validado es el ganador; un tol_cp_mm imposible detiene la validación sin ganador
    (no se calibra). No se escribe ningún .ork."""
    from dataclasses import replace
    from sensor_tipo2.barrido import completar
    from sensor_tipo2.verificacion import validar_ganadores
    rk = completar(cfg, pd.DataFrame([_fila(cfg, *c) for c in CANDIDATOS[:3]]))
    res = validar_ganadores(cfg, puente, rk, log=lambda *_: None)
    assert res.ganador == rk.iloc[0]["cand_id"] and not res.cp_fuera_de_tolerancia
    rest = cfg.restricciones
    v = res.tabla.iloc[0]
    assert rest.SM_min <= v["SM_or_cal"] <= rest.SM_max
    estricta = replace(cfg, ejecucion={**cfg.ejecucion, "validacion_or": {"N_ganadores": 2, "tol_cp_mm": 1e-6}})
    res2 = validar_ganadores(estricta, puente, rk, log=lambda *_: None)
    assert res2.ganador is None and len(res2.cp_fuera_de_tolerancia) == 2
    assert not list(tmp_path.rglob("*.ork"))
