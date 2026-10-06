"""Aleta sobre el tubo de cola: geometría, masa, Barrowman (contra OpenRocket 24.12) y flutter."""

import math

import pytest

from sensor_tipo2.aletas import area_centroide, barrowman, construir, envolvente, velocidad_flutter
from sensor_tipo2.config import AletaSpec
from sensor_tipo2.geometria import construir_cuerpo

MM = 1e-3


def test_trapecio_area_centroide_masa(cfg_libre, base_ork):
    cfg = cfg_libre
    c, a = base_ork
    cu = construir_cuerpo(cfg, c)
    g, mot = construir(cfg, cu.perfil, a)
    assert not mot
    c_r, c_t, x_s, h = 0.05, 0.035, 0.015, 0.02
    assert (g.c_r, g.c_t, g.x_s, g.h) == pytest.approx((c_r, c_t, x_s, h))
    assert g.x_LE == pytest.approx(0.35) and g.x_TE == pytest.approx(0.4) and g.r_tip == pytest.approx(0.035)
    assert g.area == pytest.approx((c_r + c_t) / 2 * h)
    x_bar = (c_r**2 + c_r * c_t + c_t**2 + x_s * (c_r + 2 * c_t)) / (3 * (c_r + c_t))
    assert g.x_cg == pytest.approx(g.x_LE + x_bar)
    p = cfg.aleta
    assert g.masa == pytest.approx(p.n * p.rho * p.phi * p.t * g.area)
    assert area_centroide(g.poligono)[0] == pytest.approx(g.area)


def test_barrowman_igual_a_openrocket(cfg_ork, base_ork):
    """OpenRocket 24.12, .ork tipo 2 a M = 0.3: aletas C_Nα = 0.884731 en x = 367.780 mm."""
    c, a = base_ork
    cu = construir_cuerpo(cfg_ork, c)
    g, _ = construir(cfg_ork, cu.perfil, a)
    f = barrowman(g, math.pi * 0.07**2 / 4, 0.3, 48)
    assert f.CNa == pytest.approx(0.884731, rel=1e-4)
    assert f.x_CP / MM == pytest.approx(367.780, abs=0.05)


def test_validaciones(cfg_libre, base_ork):
    cfg = cfg_libre
    c, _ = base_ork
    p = construir_cuerpo(cfg, c).perfil
    _, mot = construir(cfg, p, AletaSpec(0.1, 0.4, 1.0, 1.5))
    assert "c_t_menor_minimo" in mot
    _, mot = construir(cfg, p, AletaSpec(1.0, 0.4, 1.0, 0.5))
    assert "h_menor_minimo" in mot


def test_tope_del_D_acostado(base_ork):
    """4 aletas guardadas a 45°: D_acostado = max(D, √2 r_tip); r_tip/R = √2 es el límite."""
    from sensor_tipo2.config import cargar
    from .conftest import raw_libre
    raw = raw_libre()
    raw["restricciones"]["D_acostado_max_rel_D"] = 1.0
    cfg = cargar(raw)
    c, _ = base_ork
    p = construir_cuerpo(cfg, c).perfil
    g, mot = construir(cfg, p, AletaSpec(1.0, 0.7, 1.0, 1.41421356))
    assert not mot and max(envolvente(g.r_tip, p.R, 4, math.radians(45))) == pytest.approx(p.D)
    _, mot = construir(cfg, p, AletaSpec(1.0, 0.7, 1.0, 1.45))
    assert mot == ["D_acostado_mayor_maximo"]


def test_tope_de_r_tip_y_cuerda_minima(cfg):
    """Spec v2: r_tip ≤ 1.2 R (r_tip_sobre_tope) y c_r ≥ 50 mm (c_r_menor_minimo); AR = h²/A_f."""
    from sensor_tipo2.config import CuerpoSpec
    p = construir_cuerpo(cfg, CuerpoSpec(400, 90, 1.0, "conica", None, 0.0, 15.75, 140)).perfil
    g, mot = construir(cfg, p, AletaSpec(0.7, 0.7, 1.0, 1.2))
    assert not mot and g.AR == pytest.approx(g.h**2 / g.area)
    assert construir(cfg, p, AletaSpec(0.7, 0.7, 1.0, 1.25))[1] == ["r_tip_sobre_tope"]
    assert "c_r_menor_minimo" in construir(cfg, p, AletaSpec(0.3, 0.7, 1.0, 1.2))[1]  # 42 mm < 50 mm


def test_T6_malla_sobre_el_tope_falla_al_cargar(raw):
    from sensor_tipo2.config import ConfigError, cargar
    raw["malla"]["r_tip_rel_R"] = [1.0, 1.2, 1.3]
    with pytest.raises(ConfigError, match="r_tip_rel_R = 1.3 supera el tope"):
        cargar(raw)


def test_envolvente_y_flutter(cfg_libre, base_ork):
    cfg = cfg_libre
    assert envolvente(0.05, 0.035, 4, math.radians(45)) == pytest.approx((2 * 0.05 * math.cos(math.pi / 4),) * 2)
    assert envolvente(0.03, 0.035, 4, 0.0) == pytest.approx((0.07, 0.07))
    c, a = base_ork
    p = construir_cuerpo(cfg, c).perfil
    g1, _ = construir(cfg, p, AletaSpec(1.0, 0.7, 1.0, 2.0))
    V = velocidad_flutter(g1, cfg.vuelo.altitud)
    assert V > 3 * cfg.vuelo.V
