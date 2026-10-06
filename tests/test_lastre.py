"""Llenado de lastre: restricciones del criterio de dbf-sensor sobre candidatos reales."""

import copy
import math

import pytest

from sensor_tipo2.aletas import construir
from sensor_tipo2.config import AletaSpec, ConfigError, CuerpoSpec, cargar
from sensor_tipo2.geometria import construir_cuerpo
from sensor_tipo2.lastre import llenar
from sensor_tipo2.sustituto import cp_sustituto

from .conftest import FACTIBLE

CUERPO, ALETA = CuerpoSpec(*FACTIBLE[0]), AletaSpec(*FACTIBLE[1])


def _llenar(cfg, a=ALETA, c=CUERPO, exigir=True):
    cu = construir_cuerpo(cfg, c)
    g, mot = construir(cfg, cu.perfil, a)
    assert not (exigir and mot), mot
    cp, m = cp_sustituto(cu, g, cfg.vuelo.mach, 48, cfg.restricciones.eps_CN)
    return cu, cp, llenar(cfg, cu, g.masa, g.x_cg, cp.x_CP, cp.CNa)


def test_cumple_restricciones(cfg):
    cu, cp, Ll = _llenar(cfg)
    assert Ll.ok, Ll.banderas
    f = Ll.fila
    r = cfg.restricciones
    assert r.SM_min - 1e-6 <= f["SM_cal"] <= r.SM_max
    assert f["m_total_g"] <= (cfg.m_max + cfg.numerico.tol_masa_max) * 1e3
    assert f["tol_amarre_mm"] >= cfg.remolque.tol_min * 1e3 - 1e-3
    if Ll.ell2 > 0:  # el tapón trasero solo aparece con el delantero lleno
        assert Ll.ell == pytest.approx(cu.lim.ell_geo, abs=1e-4)


@pytest.mark.parametrize("m_max", [2000, 3000])
def test_tope_de_masa(raw, m_max):
    raw["masa"]["m_max_g"] = m_max
    raw["remolque"]["tol_amarre_min_mm"] = None
    raw["restricciones"]["SM_max_cal"] = 50  # con poco plomo estas aletas dan SM > 2
    _, _, Ll = _llenar(cargar(raw))
    assert Ll.ok and Ll.restriccion_activa == "masa"
    assert Ll.fila["m_total_g"] == pytest.approx(m_max, abs=0.1)


def test_mas_masa_permitida_no_reduce_la_masa(raw):
    masas = []
    for m_max in (2000, 4000, 8000, 20000):
        raw["masa"]["m_max_g"] = m_max
        Ll = _llenar(cargar(raw))[2]
        masas.append(Ll.fila["m_total_g"] if Ll.ok else 0.0)
    assert masas == sorted(masas)


def test_sin_tolerancia_de_amarre_limita_SM_o_volumen(raw):
    raw["masa"]["m_max_g"] = 1e6
    raw["remolque"]["tol_amarre_min_mm"] = None
    raw["restricciones"]["SM_max_cal"] = 50
    _, _, Ll = _llenar(cargar(raw))
    assert Ll.ok and Ll.restriccion_activa in ("SM_min", "volumen")


def test_cuerpo_del_ork_con_aletas_en_el_calibre_no_estabiliza(cfg_libre, base_ork):
    """Con el cuerpo del .ork (tubo de 50 mm) y r_tip = R el CP queda delante de la nariz."""
    c, a = base_ork
    cu, cp, Ll = _llenar(cfg_libre, a, c)
    assert cp.x_CP < 0
    assert not Ll.ok and "SM_inalcanzable_por_geometria" in Ll.banderas
    assert Ll.fila["SM_inf_cal"] < cfg_libre.restricciones.SM_min


def test_tubo_largo_estabiliza_con_aletas_en_el_calibre(cfg):
    """Spec v2 §2: la transición (−) y las aletas (+) separadas por un tubo largo forman un par que
    mueve el CP atrás aunque su fuerza neta siga siendo negativa. Con r_tip = R hay llenado factible."""
    cu, cp, Ll = _llenar(cfg, AletaSpec(0.7, 0.7, 1.0, 1.0), CUERPO)
    neta = sum(p.CNa for p in cp.partes if p.nombre in ("cola", "aletas"))
    assert neta < 0 and cp.x_CP > 0.4 * cu.perfil.L
    assert Ll.ok


# --------------------------------------------------------------------------- T2: llenado idéntico a v1


def test_T2_llenado_identico_a_la_version_anterior():
    """230 candidatos con L_c > 0 calculados con el código anterior a la spec v2 (4affe21):
    con electronica.r_min_mm = null el llenado da lo mismo."""
    import pandas as pd
    from pathlib import Path
    from .conftest import raw_libre
    ref = pd.read_csv(Path(__file__).parent / "datos" / "llenado_v1.csv", keep_default_na=True)
    raw = raw_libre()
    raw["masa"]["m_max_g"] = 10150  # la masa de la configuración v1 con la que se generaron
    raw["masas_puntuales"] = [{"nombre": "herraje_remolque", "masa_g": 15, "x_mm": 40}]  # como en v1
    raw["restricciones"]["cola_base_roma"]["infactible"] = False
    cfg = cargar(raw)
    assert cfg.electronica.r_min is None
    for _, r in ref.iterrows():
        Ln = r.L_n_rel_D * r.D_mm
        L_disp = r.L_mm - Ln - r.L_tc_mm
        par = None if pd.isna(r.parametro) else float(r.parametro)
        c = CuerpoSpec(r.L_mm, r.D_mm, r.L_n_rel_D, r.forma, par, (L_disp - r.L_t_mm) / L_disp, r.d_tc_mm, r.L_tc_mm)
        cu, cp, Ll = _llenar(cfg, AletaSpec(r.mu_cr, r.gamma_ct, r.sigma_flecha, r.r_tip_rel_R), c, exigir=False)
        assert cp.x_CP * 1e3 == pytest.approx(r.x_CP_mm, abs=1e-6)
        assert Ll.ok == r.ok and Ll.restriccion_activa == r.restriccion
        assert ";".join(Ll.banderas) == ("" if pd.isna(r.banderas) else r.banderas)
        for k in ("m_total_g", "SM_cal", "x_CG_mm", "tol_amarre_mm", "ell_mm", "ell_trasero_mm", "SM_inf_cal"):
            v = Ll.fila.get(k, math.nan)
            # SM_inf (informativo) sale de minimize_scalar con xatol = 0.01 mm: un redondeo de 1e-16 en
            # L_t (ahora desde f_c) cambia su ruta; el llenado en sí coincide a 1e-7
            tol = {"rel": 1e-5, "abs": 1e-4} if k == "SM_inf_cal" else {"rel": 1e-7, "abs": 1e-6}
            assert (math.isnan(v) and math.isnan(r[k])) or v == pytest.approx(r[k], **tol), k


# --------------------------------------------------------------------------- T5: electrónica con r_min


@pytest.mark.parametrize("r_min_mm", [20.0, 30.0, 40.0])
def test_T5_electronica_en_tramo_con_radio_suficiente(raw, r_min_mm):
    raw["electronica"]["r_min_mm"] = r_min_mm
    cfg = cargar(raw)
    for fc in (0.0, 0.5):
        c = CuerpoSpec(400, 90, 1.0, "conica", None, fc, 15.75, 140)
        cu = construir_cuerpo(cfg, c)
        assert cu.ok, cu.motivos
        g, _ = construir(cfg, cu.perfil, ALETA)
        cp, _ = cp_sustituto(cu, g, cfg.vuelo.mach, 48, cfg.restricciones.eps_CN)
        Ll = llenar(cfg, cu, g.masa, g.x_cg, cp.x_CP, cp.CNa)
        if not math.isfinite(Ll.ell):
            continue
        x_e = float(Ll.modelo.x_e(Ll.ell))
        tramo = (cu.cav.x >= x_e) & (cu.cav.x <= x_e + cfg.electronica.Le)
        assert cu.cav.r_i[tramo].min() >= r_min_mm * 1e-3 - 1e-9


def test_T5_electronica_que_no_cabe(raw):
    raw["electronica"]["r_min_mm"] = 60.0  # más que el radio interior de cualquier cuerpo de D = 90
    cfg = cargar(raw)
    cu = construir_cuerpo(cfg, CuerpoSpec(400, 90, 1.0, "conica", None, 0.0, 15.75, 140))
    assert cu.motivos == ["electronica_no_cabe"]
    from sensor_tipo2.barrido import evaluar_cuerpo
    f = evaluar_cuerpo(cfg, cu.spec, [ALETA])[0]
    assert not f["factible"] and f["motivos"] == "electronica_no_cabe"


# --------------------------------------------------------------------------- herraje de remolque en el CG


def test_herraje_en_el_CG(raw):
    """El herraje (x_mm: cg) va en el CG del sensor: suma masa pero no mueve el CG, y en las masas que
    se cargan en OpenRocket queda exactamente en el x_CG final (el amarre)."""
    import numpy as np
    from sensor_tipo2.exportar import detalle, masas_del_llenado
    from .conftest import FACTIBLE
    assert raw["masas_puntuales"][0]["x_mm"] == "cg"
    cfg = cargar(raw)
    sin = copy.deepcopy(raw)
    sin["masas_puntuales"] = []
    cfg_sin = cargar(sin)
    c, a = CuerpoSpec(*FACTIBLE[0]), AletaSpec(*FACTIBLE[1])
    d, d_sin = detalle(cfg, c, a), detalle(cfg_sin, c, a)
    m, m_sin = d.llenado.modelo, d_sin.llenado.modelo
    ells = np.linspace(0.0, d.cu.lim.ell_geo, 7)
    assert np.allclose(m.x_CG(ells), m_sin.x_CG(ells), rtol=0, atol=1e-12)
    assert np.allclose(m.m(ells) - m_sin.m(ells), 0.015)
    assert np.allclose(m.x_CG_con_trasero(ells[-1], 0.02), m_sin.x_CG_con_trasero(ells[-1], 0.02), rtol=0, atol=1e-12)
    # en la lista de masas para OpenRocket el herraje está en el CG final, y el CG del conjunto cierra
    f = d.llenado.fila
    masas = {n: (mm, x) for n, mm, x in masas_del_llenado(cfg, d)}
    assert masas["herraje_remolque"][1] == pytest.approx(f["x_CG_mm"] * 1e-3, abs=1e-12)
    M = d.cu.M_casco + d.g.masa * d.g.x_cg + sum(mm * x for mm, x in masas.values())
    mt = d.cu.m_casco + d.g.masa + sum(mm for mm, _ in masas.values())
    assert mt * 1e3 == pytest.approx(f["m_total_g"], rel=1e-9)
    assert M / mt * 1e3 == pytest.approx(f["x_CG_mm"], abs=1e-6)


def test_masa_puntual_x_invalida(raw):
    raw["masas_puntuales"] = [{"nombre": "x", "masa_g": 10, "x_mm": "centro"}]
    with pytest.raises(ConfigError, match="masas_puntuales"):
        cargar(raw)
