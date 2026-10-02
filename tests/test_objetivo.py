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
