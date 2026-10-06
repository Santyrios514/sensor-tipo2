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
    assert ranking["cand_id"].is_unique


def test_T3_factibles_cumplen_la_spec(chica, ranking):
    """En todo candidato factible: L_n + L_c + L_t + L_tc = L, L_c ≥ 0, L_n/D ≥ mínimo,
    r_tip ≤ 1.2 R y c_r ≥ c_r,min."""
    fac = ranking[ranking["factible"]]
    r = chica.restricciones
    suma = fac["L_n_mm"] + fac["L_c_mm"] + fac["L_t_mm"] + fac["L_tc_mm"]
    assert (abs(suma - fac["L_mm"]) < 1e-6).all() and (fac["L_c_mm"] >= 0).all()
    assert (fac["L_n_rel_D"] >= r.L_n_rel_D_min - 1e-12).all()
    assert (fac["r_tip_mm"] <= r.r_tip_rel_R_max * fac["D_mm"] / 2 + 1e-6).all()
    assert (fac["c_r_mm"] >= r.c_r_min * 1e3 - 1e-6).all()
    assert (fac["f_cil"] == 0).any()  # hay cuerpos abombados factibles
    # banderas informativas de la spec v2
    for col in ("frac_h_fuera_sombra", "aletas_en_estela", "tubo_esbelto", "esbeltez_tubo", "AR_aleta"):
        assert fac[col].notna().all(), col
    assert (fac["frac_h_fuera_sombra"].between(0, 1)).all()


def test_prefiltro(chica):
    """Los cuerpos y aletas imposibles se descartan antes de evaluar (spec v2 §5, paso 1)."""
    from sensor_tipo2.barrido import grupos_malla, n_candidatos
    total = len(chica.cuerpos()) * len(chica.aletas())
    gr = grupos_malla(chica)
    assert 0 < n_candidatos(gr) < total
    assert all(a.mu_cr * c.L_tc >= chica.restricciones.c_r_min - 1e-12 for c, al in gr.items() for a in al)
    assert all(c.d_tc >= chica.restricciones.d_tc_min - 1e-12 for c in gr)


def test_refinamiento_no_repite(chica, ranking):
    base = ranking[~ranking["refinamiento"]]
    gr = grupos_refinamiento(chica, base, 2)
    ids = {cand_id(c, a) for c, al in gr.items() for a in al}
    assert ids and not ids & set(base["cand_id"])


def test_refinamiento_solo_grupos_A_B_C(chica, ranking):
    """Spec v2 rev. 2 §5: el refinamiento no toca μ, γ, σ, la forma de la transición ni L."""
    base = ranking[~ranking["refinamiento"]]
    semillas = base[base["factible"]].head(2)
    fijos = {(r.mu_cr, r.gamma_ct, r.sigma_flecha, r.cola_forma, r.L_mm) for r in semillas.itertuples()}
    gr = grupos_refinamiento(chica, base, 2)
    for c, al in gr.items():
        for a in al:
            assert (a.mu_cr, a.gamma_ct, a.sigma_flecha, c.forma, c.L_mm) in fijos
            assert a.r_tip_rel_R <= chica.restricciones.r_tip_rel_R_max + 1e-12


def test_specs_de_fila_ida_y_vuelta(chica, ranking):
    for _, f in ranking.head(5).iterrows():
        c, a = specs_de_fila(f)
        assert cand_id(c, a) == f["cand_id"]


def test_config_invalida(raw):
    raw["malla"]["L_total_mm"] = [420]
    with pytest.raises(ConfigError, match="L_max"):
        cargar(raw)
    raw2 = malla_chica(raw, k=[1.2])
    raw2["malla"]["L_total_mm"] = [400]
    with pytest.raises(ConfigError, match="malla.k"):
        cargar(raw2)


def test_D_acostado_invalido(raw):
    """El D acostado nunca es menor que D: un tope < 1 es un error de configuración."""
    raw["restricciones"]["D_acostado_max_rel_D"] = 0.9
    with pytest.raises(ConfigError, match="D_acostado_max_rel_D"):
        cargar(raw)


def test_tamano_de_la_malla(cfg):
    from sensor_tipo2.barrido import n_candidatos
    assert n_candidatos(grupos_malla(cfg, prefiltrar=False)) == cfg.n_evaluaciones
    assert n_candidatos(grupos_malla(cfg)) <= int(cfg.ejecucion["max_evaluaciones"])
