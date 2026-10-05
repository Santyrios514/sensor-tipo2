"""Objetivo lexicográfico ponderado, tolerancias y Pareto."""

import numpy as np
import pytest

from sensor_tipo2.objetivo import J, f_valores, pareto, tolerancias
from sensor_tipo2.sustituto import MOTIVO_CN, ajustar, combinar
from sensor_tipo2.geometria import Contribucion


def test_masa_manda_sobre_diametro(cfg):
    """Por encima de la tolerancia implícita (~101 g) la masa decide aunque el otro sea mucho más chico."""
    F = f_valores(cfg, [cfg.m_max, cfg.m_max - 0.2], [0.20, 0.07], [0.9, 0.3], [1.5, 1.5])
    j = J(cfg, F)
    assert j[0] < j[1]


def test_cotas_del_diametro(cfg, raw):
    """Spec v2: criterio 2 = D aparente; con r_tip ≤ 1.2 R el rango es [min D, 1.2 max D]."""
    from sensor_tipo2.config import cargar
    assert cfg.objetivo.col_diametro == "D_ap_mm"
    assert cfg.cotas["D_hi"] == pytest.approx(1.2 * max(cfg.malla["D_mm"]) * 1e-3)
    assert cfg.cotas["k_lo"] == pytest.approx(min(cfg.malla["k"]))
    raw["objetivo"]["diametro"] = "acostado"
    c2 = cargar(raw)
    assert c2.objetivo.col_diametro == "D_acostado_mm"
    assert c2.cotas["D_hi"] == pytest.approx(max(raw["malla"]["D_mm"]) * 1e-3)  # √2 · 1.2 R < D


def test_f_en_rango(cfg):
    F = f_valores(cfg, [0, 2 * cfg.m_max], [0.0, 1.0], [0.0, 1.0], [-5, 5])
    assert np.all((F >= 0) & (F <= 1))


def test_tolerancias(cfg):
    t = tolerancias(cfg)
    assert t["eps"]["f1_masa"] == pytest.approx((1e4 + 1e2 + 1) / 1e6)
    assert t["fisicas"]["masa_g"] == pytest.approx(t["eps"]["f1_masa"] * cfg.m_max * 1e3)


def test_pareto():
    F = np.array([[0, 1], [1, 0], [1, 1], [0.5, 0.5], [0, 1]])
    assert list(pareto(F)) == [True, True, False, True, True]


def test_guarda_del_CP():
    res, m = combinar((Contribucion("nariz", 2.0, 0.02), Contribucion("cola", -2.0, 0.3),
                       Contribucion("aletas", 0.1, 0.38)), 0.5)
    assert res is None and m == MOTIVO_CN


def test_calibracion_lineal():
    x = np.linspace(0.1, 0.3, 10)
    cal = ajustar(x, 0.002 + 1.01 * x)
    assert cal.alpha == pytest.approx(0.002) and cal.beta == pytest.approx(1.01) and cal.R2 == pytest.approx(1.0)
