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
    spec = CuerpoSpec(400, 80, 1.0, "conica", None, 0.8, 0.4, 50)
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
    spec = CuerpoSpec(400, 80, 1.0, "conica", None, 0.8, 0.4, 50)
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


def test_motivos_de_cuerpo(cfg):
    assert "d_tubo_cola_bajo_minimo" in construir_cuerpo(cfg, CuerpoSpec(400, 60, 1, "conica", None, 1, 0.3, 50)).motivos
    assert "cola_base_roma" in construir_cuerpo(cfg, CuerpoSpec(400, 80, 1, "conica", None, 0.4, 0.4, 50)).motivos
    assert "L_c_no_positivo" in construir_cuerpo(cfg, CuerpoSpec(400, 90, 1.5, "conica", None, 1.25, 0.4, 160)).motivos
