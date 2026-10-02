"""Barrido: reproducibilidad, orden del ranking, refinamiento y configuración."""

import pandas as pd
import pytest

from sensor_tipo2.barrido import grupos_malla, grupos_refinamiento, optimizar, specs_de_fila
from sensor_tipo2.config import ConfigError, cand_id, cargar

from .conftest import malla_chica


@pytest.fixture(scope="module")
def chica():
    from .conftest import raw_opt
    return cargar(malla_chica(raw_opt()))


@pytest.fixture(scope="module")
def ranking(chica):
    return optimizar(chica, refinar=True, procesos=1)


def test_reproducible_en_paralelo(chica, ranking):
    df2 = optimizar(chica, refinar=True, procesos=2)
    pd.testing.assert_series_equal(ranking["cand_id"], df2["cand_id"])
    pd.testing.assert_series_equal(ranking["J"], df2["J"])


def test_orden_y_factibles(chica, ranking):
    fac = ranking[ranking["factible"]]
    assert len(fac) > 0
    assert list(fac.index) == list(range(len(fac)))
    assert fac["J"].is_monotonic_increasing
    r = chica.restricciones
    assert fac["SM_cal"].between(r.SM_min - 1e-6, r.SM_max).all()
    assert (fac["m_total_g"] <= (chica.m_max + chica.numerico.tol_masa_max) * 1e3).all()
    assert not ranking.loc[ranking["r_tip_rel_R"] == 1.0, "factible"].any()
    assert not ranking.loc[ranking["r_tip_rel_R"] == 1.6, "factible"].any()  # se sale del D acostado
    assert (fac["D_acostado_mm"] <= fac["D_mm"] + 1e-6).all()
    assert ranking["cand_id"].is_unique


def test_refinamiento_no_repite(chica, ranking):
    base = ranking[~ranking["refinamiento"]]
    gr = grupos_refinamiento(chica, base, 2)
    ids = {cand_id(c, a) for c, al in gr.items() for a in al}
    assert ids and not ids & set(base["cand_id"])


def test_specs_de_fila_ida_y_vuelta(chica, ranking):
    for _, f in ranking.head(5).iterrows():
        c, a = specs_de_fila(f)
        assert cand_id(c, a) == f["cand_id"]


def test_config_invalida(raw):
    raw["malla"]["L_total_mm"] = [420]
    with pytest.raises(ConfigError, match="L_max"):
        cargar(raw)
    raw2 = malla_chica(raw, d_tc_mm=[-5])
    raw2["malla"]["L_total_mm"] = [400]
    with pytest.raises(ConfigError, match="malla.d_tc_mm"):
        cargar(raw2)


def test_D_acostado_invalido(raw):
    """El D acostado nunca es menor que D: un tope < 1 es un error de configuración."""
    raw["restricciones"]["D_acostado_max_rel_D"] = 0.9
    with pytest.raises(ConfigError, match="D_acostado_max_rel_D"):
        cargar(raw)


def test_tamano_de_la_malla(cfg):
    gr = grupos_malla(cfg)
    assert len(gr) * len(next(iter(gr.values()))) == cfg.n_evaluaciones
