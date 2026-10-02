"""Llenado de lastre: restricciones del criterio de dbf-sensor sobre candidatos reales."""


import pytest

from sensor_tipo2.aletas import construir
from sensor_tipo2.config import AletaSpec, CuerpoSpec, cargar
from sensor_tipo2.geometria import construir_cuerpo
from sensor_tipo2.lastre import llenar
from sensor_tipo2.sustituto import cp_sustituto

from .conftest import FACTIBLE

CUERPO, ALETA = CuerpoSpec(*FACTIBLE[0]), AletaSpec(*FACTIBLE[1])


def _llenar(cfg, a=ALETA, c=CUERPO):
    cu = construir_cuerpo(cfg, c)
    g, mot = construir(cfg, cu.perfil, a)
    assert not mot
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


def test_aletas_dentro_del_diametro_no_estabilizan(cfg, base_ork):
    """Con r_tip = R el CP queda delante de la nariz: ni un lastre infinitamente denso da SM ≥ 1."""
    c, a = base_ork
    cu, cp, Ll = _llenar(cfg, a, c)
    assert cp.x_CP < 0
    assert not Ll.ok and "SM_inalcanzable_por_geometria" in Ll.banderas
    assert Ll.fila["SM_inf_cal"] < cfg.restricciones.SM_min
