"""Objetivo lexicográfico ponderado (T9), tolerancias y Pareto."""

import numpy as np
import pytest

from sensor_tipo2.objetivo import J, columnas, f_valores, pareto, tolerancias
from sensor_tipo2.sustituto import MOTIVO_CN, combinar
from sensor_tipo2.geometria import Contribucion


def _J(cfg, m, D, k, SM):
    return J(cfg, f_valores(cfg, m, D, k, SM))


def test_orden_por_defecto(cfg):
    assert cfg.objetivo.orden == ("masa", "diametro", "SM", "taper")
    assert columnas(cfg) == ["f1_masa", "f2_diametro", "f3_SM", "f4_taper"]
    assert cfg.objetivo.f_SM_modo == "max"


def test_T9_lexicografico(cfg):
    """T9: si A mejora a B en un criterio por más de su tolerancia implícita, J_A < J_B sean cuales
    sean los criterios inferiores (en el peor caso: B es perfecto en todos ellos)."""
    t = tolerancias(cfg)["fisicas"]
    c, r = cfg.cotas, cfg.restricciones
    mejor_D, peor_D = c["D_lo"], c["D_hi"]
    mejor_k, peor_k = c["k_lo"], 1.0
    # masa
    dm = (t["masa"] + 1) * 1e-3
    j = _J(cfg, [cfg.m_max, cfg.m_max - dm], [peor_D, mejor_D], [peor_k, mejor_k], [r.SM_min, r.SM_max])
    assert j[0] < j[1]
    # diámetro, con la misma masa
    dD = (t["diametro"] + 0.01) * 1e-3
    j = _J(cfg, [cfg.m_max] * 2, [mejor_D, mejor_D + dD], [peor_k, mejor_k], [r.SM_min, r.SM_max])
    assert j[0] < j[1]
    # SM (modo max: se prefiere el mayor), con la misma masa y diámetro
    dS = t["SM"] + 1e-3
    j = _J(cfg, [cfg.m_max] * 2, [mejor_D] * 2, [peor_k, mejor_k], [r.SM_max, r.SM_max - dS])
    assert j[0] < j[1]
    # dentro de la tolerancia, el criterio inferior sí puede decidir
    j = _J(cfg, [cfg.m_max] * 2, [mejor_D] * 2, [peor_k, mejor_k], [r.SM_max, r.SM_max - 0.5 * t["SM"]])
    assert j[0] > j[1]


def test_orden_dbf_sensor(raw):
    """objetivo.orden permite volver al orden de dbf-sensor (taper antes que SM)."""
    from sensor_tipo2.config import ConfigError, cargar
    raw["objetivo"]["orden"] = ["masa", "diametro", "taper", "SM"]
    c = cargar(raw)
    assert columnas(c) == ["f1_masa", "f2_diametro", "f3_taper", "f4_SM"]
    j = _J(c, [c.m_max] * 2, [c.cotas["D_lo"]] * 2, [0.2, 0.3], [1.0, 2.0])
    assert j[0] < j[1]  # ahora el taper manda sobre el SM
    raw["objetivo"]["orden"] = ["masa", "diametro", "SM"]
    with pytest.raises(ConfigError, match="permutación"):
        cargar(raw)


def test_masa_manda_sobre_diametro(cfg):
    """Por encima de la tolerancia implícita (~116 g) la masa decide aunque el otro sea mucho más chico."""
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
    assert t["fisicas"]["masa"] == pytest.approx(t["eps"]["f1_masa"] * cfg.m_max * 1e3)
    assert t["fisicas"]["SM"] == pytest.approx(t["eps"]["f3_SM"] * (cfg.restricciones.SM_max - cfg.restricciones.SM_min))


def test_pareto():
    F = np.array([[0, 1], [1, 0], [1, 1], [0.5, 0.5], [0, 1]])
    assert list(pareto(F)) == [True, True, False, True, True]


def test_guarda_del_CP():
    res, m = combinar((Contribucion("nariz", 2.0, 0.02), Contribucion("cola", -2.0, 0.3),
                       Contribucion("aletas", 0.1, 0.38)), 0.5)
    assert res is None and m == MOTIVO_CN
