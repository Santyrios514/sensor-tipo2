"""Perfiles, cavidad, masas de pared y Barrowman del cuerpo contra soluciones analíticas."""

import math

import numpy as np
import pytest

from sensor_tipo2.config import CuerpoSpec
from sensor_tipo2.geometria import (construir_cuerpo, cuerpo_revolucion, erosionar, fraccion_base_roma,
                                    k_efectivo, malla)
from sensor_tipo2.perfiles import FORMAS, Perfil, forma_or, longitud_recorte

MM = 1e-3


def _perfil(forma="elipsoide", recortada=True, param=None, k=30 / 70):
    return Perfil(L=0.4, D=0.07, Ln=0.065, forma_n="elipsoide", param_n=None, Lt=0.065, forma_t=forma,
                  param_t=param, recortada=recortada, d_tc=k * 0.07, L_tc=0.05)


@pytest.mark.parametrize("forma", FORMAS)
@pytest.mark.parametrize("recortada", [True, False])
def test_transicion_continua_y_monotona(forma, recortada):
    p = _perfil(forma, recortada, 0.5 if forma == "potencia" else None)
    x = np.linspace(p.x_t0, p.x_tc0, 2001)
    r = p.radio(x)
    assert p.radio(np.array([p.x_t0]))[0] == pytest.approx(p.R)
    assert p.radio_transicion(np.array([p.Lt]))[0] == pytest.approx(p.r_tc, abs=1e-9)
    assert np.all(np.diff(r[1:-1]) <= 1e-12)
    assert p.radio(np.array([0.3999]))[0] == pytest.approx(p.r_tc)


def test_nariz_elipsoide_y_recorte():
    p = _perfil()
    x = np.array([0.0, p.Ln / 2, p.Ln])
    assert p.radio(x) == pytest.approx([0.0, p.R * math.sqrt(0.75), p.R])
    c = longitud_recorte("elipsoide", p.R, p.r_tc, p.Lt, None)
    assert float(forma_or(np.array([c]), "elipsoide", p.R, c + p.Lt, None)[0]) == pytest.approx(p.r_tc)
    assert list(p.estacion(np.array([0.0, p.Ln, p.x_t0, p.x_tc0]))) == [0, 1, 2, 3]


def test_erosion_compone():
    x = malla(0.2, 1e-4)
    r = np.where(x < 0.05, 0.03 * np.sqrt(np.clip(2 * x / 0.05 - (x / 0.05) ** 2, 0, None)), 0.03)
    a = erosionar(erosionar(r, 0.3e-3, 1e-4), 1.5e-3, 1e-4)
    b = erosionar(r, 1.8e-3, 1e-4)
    lejos = x > 0.005  # en la punta (pendiente infinita) la erosión discreta no compone exacto
    assert np.max(np.abs(a - b)[lejos]) < 2e-5


def test_cavidad_y_masa_del_cilindro(cfg):
    spec = CuerpoSpec(400, 80, 1.0, "conica", None, 0.75, 32, 50)
    cu = construir_cuerpo(cfg, spec)
    p, cav = cu.perfil, cu.cav
    T = sum(c.t for c in cfg.pared["cuerpo"])
    a, b = p.Ln + 0.01, p.x_t0 - 0.01
    assert float(cav.vol(a, b)) == pytest.approx(math.pi * (p.R - T) ** 2 * (b - a), rel=1e-6)
    assert float(cav.mom(a, b)) == pytest.approx(math.pi * (p.R - T) ** 2 * (b**2 - a**2) / 2, rel=1e-6)
    # masa de la pared del cuerpo (las ventanas de las uniones cambian < 1 %)
    m_an = sum(c.rho * c.phi for c in cfg.pared["cuerpo"][:1]) * math.pi * (p.R**2 - (p.R - 0.3e-3) ** 2) * p.L_c
    m_cap = next(c.m for c in cu.capas if c.estacion == "cuerpo" and c.capa_idx == 1)
    assert m_cap == pytest.approx(m_an, rel=0.01)


def test_barrowman_nariz_elipsoide_y_transicion_conica(cfg):
    spec = CuerpoSpec(400, 80, 1.0, "conica", None, 0.75, 32, 50)
    cu = construir_cuerpo(cfg, spec)
    p = cu.perfil
    nariz, cola = cu.partes
    assert nariz.CNa == pytest.approx(2.0) and nariz.x == pytest.approx(p.Ln / 3, rel=1e-4)
    k = p.k
    assert cola.CNa == pytest.approx(2 * (k**2 - 1))
    assert cola.x == pytest.approx(p.x_t0 + p.Lt / 3 * (1 + 2 * k) / (1 + k), rel=1e-5)


def test_cuerpo_revolucion_sin_cambio_de_area():
    x = np.linspace(0, 1, 101)
    assert cuerpo_revolucion(x, np.full_like(x, 0.1), 0.2, 0.8, 1.0, "cilindro") is None


def test_base_roma_y_k_efectivo():
    assert fraccion_base_roma(3.5) == 0 and fraccion_base_roma(0.5) == 1
    assert fraccion_base_roma(2.0) == pytest.approx(0.5)
    assert k_efectivo(0.4, 0.0) == pytest.approx(0.4) and k_efectivo(0.4, 1.0) == pytest.approx(1.0)


def test_convencion_del_ork(base_ork):
    """Cuerpo central ≥ transición + tubo de cola (v1, apagada por defecto en v2): el .ork la
    cumple (220 ≥ 115); una cola larga no."""
    from sensor_tipo2.config import cargar
    from .conftest import raw_libre
    raw = raw_libre()
    raw["restricciones"]["L_c_min_rel_cola"] = 1.0
    cfg = cargar(raw)
    assert construir_cuerpo(cfg, base_ork[0]).ok
    larga = CuerpoSpec(400, 90, 0.3, "conica", None, 0.3, 20, 180)
    assert "cuerpo_central_corto" in construir_cuerpo(cfg, larga).motivos


def test_motivos_de_cuerpo(cfg):
    def motivos(*a):
        return construir_cuerpo(cfg, CuerpoSpec(*a)).motivos
    assert "d_tubo_cola_bajo_minimo" in motivos(400, 60, 1, "conica", None, 0.5, 12, 80)
    assert "cola_base_roma" in motivos(400, 80, 1, "conica", None, 0.9, 32, 50)
    assert "L_disp_no_positivo" in motivos(400, 90, 1.5, "conica", None, 0.0, 36, 300)
    assert "nariz_bajo_minimo" in motivos(400, 80, 0.75, "conica", None, 0.5, 24, 100)
    assert "k_bajo_minimo" in motivos(400, 90, 1.0, "conica", None, 0.5, 13, 100)


# --------------------------------------------------------------------------- T3 y cuerpo abombado (spec v2)


def test_T3_largos_de_la_malla(cfg):
    """L_n + L_c + L_t + L_tc = L y L_c ≥ 0 en toda la malla; L_c = f_c L_disp."""
    for c in cfg.cuerpos():
        assert c.Ln + c.L_c + c.Lt + c.L_tc == pytest.approx(c.L, abs=1e-9)
        assert c.L_c >= 0
        if c.L_disp > 0:
            assert c.L_c == pytest.approx(c.f_cil * c.L_disp)


def test_cuerpo_abombado(cfg):
    """f_c = 0: la nariz y la transición se unen en el diámetro máximo; el cuerpo cilíndrico no existe."""
    c = CuerpoSpec(400, 90, 1.0, "conica", None, 0.0, 15.75, 140)
    cu = construir_cuerpo(cfg, c)
    assert cu.ok and cu.perfil.L_c == 0 and cu.perfil.x_t0 == pytest.approx(cu.perfil.Ln)
    assert cu.masas_por_estacion()["cuerpo"][0] == 0
    x = cu.cav.x
    assert cu.cav.r_e[np.argmin(np.abs(x - cu.perfil.Ln))] == pytest.approx(cu.perfil.R, abs=1e-6)
    assert cu.cav.r_e.max() <= cu.perfil.R + 1e-12
