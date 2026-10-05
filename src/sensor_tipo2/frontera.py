"""Estudio de factibilidad (spec v2 §6): ¿hay diseño estable con las aletas dentro del calibre?

Para cada combinación de tope de r_tip/R y número de aletas n se cuenta cuántos candidatos son
factibles y cuál admite más masa. Cada cuerpo se construye una vez y se evalúa con todas las n
(la geometría del cuerpo no depende de las aletas). Las filas con tope > `r_tip_rel_R_max` o
n ≠ `geometria_fija.aletas.n` son **solo diagnóstico**: no entran al ranking ni pueden ganar.

Si no hay factibles con las restricciones reales, `casi_factibles` devuelve los candidatos más
cercanos, ordenados por su déficit de estabilidad.
"""

from __future__ import annotations

import multiprocessing as mp
from dataclasses import replace

import numpy as np
import pandas as pd

from .barrido import aleta_imposible, evaluar_cuerpo
from .config import ConfigOpt
from .geometria import construir_cuerpo, motivos_previos
from .objetivo import J, f_valores

MM, G = 1e-3, 1e-3
TOPES = (1.0, 1.1, 1.2, 1.3, 1.4, 1.6)
NS = (4, 6, 8)
COLUMNAS = ["cand_id", "n_aletas", "factible", "motivos", "m_total_g", "D_ap_mm", "D_acostado_mm", "k_efectivo",
            "SM_cal", "tol_amarre_mm", "r_tip_rel_R", "D_mm", "SM_max_alcanzable_cal", "m_en_SM_max_g",
            "tol_en_SM_max_mm", "restriccion_activa"]


def cfg_n(cfg: ConfigOpt, n: int) -> ConfigOpt:
    """La configuración con n aletas y sin tope de r_tip (para el diagnóstico)."""
    return replace(cfg, aleta=replace(cfg.aleta, n=int(n)),
                   restricciones=replace(cfg.restricciones, r_tip_rel_R_max=None))


def malla_frontera(cfg: ConfigOpt, topes=TOPES, gruesa: bool = False) -> dict:
    """La malla de la configuración con r_tip/R extendido a los topes. Con `gruesa`, un valor de
    cada dos en D, el tubo y L_tc (siempre con el mayor), como el paso (2) de la spec v2 §5."""
    m = dict(cfg.malla)
    m["r_tip_rel_R"] = sorted(set(float(v) for v in m["r_tip_rel_R"]) | set(float(t) for t in topes))
    if gruesa:
        for k in ("D_mm", cfg.var_tubo, "L_tc_mm"):
            v = sorted(float(x) for x in m[k])
            m[k] = sorted(set(v[::2]) | {v[-1]})
    return m


_CFG: ConfigOpt | None = None
_NS: tuple[int, ...] = NS
_TOPES: tuple[float, ...] = TOPES
N_CASI = 50


def _init(cfg: ConfigOpt, ns: tuple[int, ...], topes: tuple[float, ...]):
    global _CFG, _NS, _TOPES
    _CFG, _NS, _TOPES = cfg, ns, topes


def _con_J(cfg: ConfigOpt, d: pd.DataFrame) -> pd.DataFrame:
    F = f_valores(cfg, d["m_total_g"].astype(float) * G, d[cfg.objetivo.col_diametro].astype(float) * MM,
                  d["k_efectivo"].astype(float), d["SM_cal"].astype(float))
    return d.assign(J=np.where(np.isfinite(d["m_total_g"].astype(float)), J(cfg, np.nan_to_num(F)), np.nan))


def _resumen(cfg: ConfigOpt, d: pd.DataFrame, ns, topes) -> tuple[dict, pd.DataFrame]:
    """{(n, tope): (n_factibles, mejor fila)} y los casi factibles (restricciones reales) de d."""
    d = _con_J(cfg, d)
    d["factible"] = d["factible"].fillna(False).astype(bool)
    stats = {}
    for n in ns:
        dn = d[(d["n_aletas"] == n) & d["factible"]]
        for t in topes:
            dt = dn[dn["r_tip_rel_R"] <= t + 1e-9]
            mejor = dt.sort_values(["m_total_g", "J"], ascending=[False, True]).head(1) if len(dt) else None
            stats[(n, t)] = (len(dt), mejor)
    return stats, casi_factibles(cfg, d, N_CASI)


def _tarea(args):
    c, aletas = args
    cu = construir_cuerpo(_CFG, c)
    filas = []
    for n in _NS:
        for f in evaluar_cuerpo(cfg_n(_CFG, n), c, aletas, cu=cu):
            f["n_aletas"] = n
            filas.append({k: f.get(k) for k in COLUMNAS})
    return _resumen(_CFG, pd.DataFrame(filas, columns=COLUMNAS), _NS, _TOPES)


def evaluar(cfg: ConfigOpt, topes=TOPES, ns=NS, gruesa: bool = False, procesos: int | None = None,
            progreso=None) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """(tabla por (n, tope), casi factibles, número de evaluaciones). Cada proceso resume sus
    cuerpos (conteos, mejor candidato y los N_CASI más cercanos) para no guardar millones de filas."""
    m = malla_frontera(cfg, topes, gruesa)
    base = cfg_n(cfg, cfg.aleta.n)
    aletas = cfg.aletas(m)
    tareas = []
    for c in cfg.cuerpos(m):
        if motivos_previos(base, c):
            continue
        al = [a for a in aletas if not aleta_imposible(base, a, c.L_tc)]
        if al:
            tareas.append((c, al))
    n_eval = sum(len(al) for _, al in tareas) * len(ns)
    procesos = procesos or cfg.procesos
    cuentas = {(n, t): 0 for n in ns for t in topes}
    mejores: dict = {}
    casi = []

    def acumular(res):
        stats, cf = res
        for key, (cnt, fila) in stats.items():
            cuentas[key] += cnt
            if fila is not None:
                mejores.setdefault(key, []).append(fila)
        if len(cf):
            casi.append(cf)

    if procesos <= 1 or len(tareas) <= 1:
        _init(cfg, tuple(ns), tuple(topes))
        for i, t in enumerate(tareas):
            acumular(_tarea(t))
            if progreso:
                progreso(i + 1, len(tareas))
    else:
        with mp.get_context("fork").Pool(procesos, initializer=_init, initargs=(cfg, tuple(ns), tuple(topes))) as pool:
            for i, r in enumerate(pool.imap(_tarea, tareas, chunksize=4)):
                acumular(r)
                if progreso:
                    progreso(i + 1, len(tareas))
    tope_real = cfg.restricciones.r_tip_rel_R_max
    filas = []
    for n in ns:
        for t in topes:
            fila = {"n_aletas": n, "tope_r_tip_rel_R": t, "n_factibles": cuentas[(n, t)],
                    "solo_diagnostico": bool(n != cfg.aleta.n or (tope_real is not None and t > tope_real + 1e-9))}
            if mejores.get((n, t)):
                g = pd.concat(mejores[(n, t)]).sort_values(["m_total_g", "J"], ascending=[False, True]).iloc[0]
                fila.update({"m_total_max_g": g["m_total_g"], "D_ap_mm": g["D_ap_mm"], "D_mm": g["D_mm"],
                             "r_tip_rel_R": g["r_tip_rel_R"], "SM_cal": g["SM_cal"],
                             "tol_amarre_mm": g["tol_amarre_mm"], "restriccion_activa": g["restriccion_activa"],
                             "cand_id": g["cand_id"]})
            filas.append(fila)
    cf = pd.concat(casi, ignore_index=True) if casi else casi_factibles(cfg, pd.DataFrame(columns=COLUMNAS))
    if len(cf):
        cf = cf.sort_values(["deficit_SM_cal", "deficit_tol_amarre_mm", "m_en_SM_max_g"],
                            ascending=[True, True, False]).head(N_CASI).reset_index(drop=True)
    return pd.DataFrame(filas), cf, n_eval


def _falla(motivos: str) -> str:
    m = set(filter(None, (motivos or "").split(";")))
    if m & {"SM_inalcanzable", "SM_inalcanzable_por_geometria", "CN_total_pequeno"}:
        return "SM"
    if "tol_amarre_inalcanzable" in m:
        return "tolerancia_amarre"
    if "electronica_no_cabe" in m:
        return "electronica"
    if m & {"inviable_geo", "sin_region_lastre", "SM_inalcanzable_con_masa_max"}:
        return "volumen"
    if "SM_sobre_max" in m:
        return "SM_max"
    return ";".join(sorted(m)) or "?"


def casi_factibles(cfg: ConfigOpt, df: pd.DataFrame, n_max: int = 50) -> pd.DataFrame:
    """Los n_max candidatos más cercanos a ser factibles con las restricciones reales (n de la
    config y r_tip ≤ tope), ordenados por déficit de SM y luego de tolerancia de amarre."""
    rest = cfg.restricciones
    d = df[(df["n_aletas"] == cfg.aleta.n) & ~df["factible"].fillna(False).astype(bool)]
    if rest.r_tip_rel_R_max is not None:
        d = d[d["r_tip_rel_R"] <= rest.r_tip_rel_R_max + 1e-9]
    d = d[np.isfinite(pd.to_numeric(d["SM_max_alcanzable_cal"], errors="coerce").astype(float))].copy()
    if d.empty:
        return pd.DataFrame(columns=["cand_id", "restriccion_que_falla", "deficit_SM_cal", "deficit_tol_amarre_mm"])
    d["restriccion_que_falla"] = d["motivos"].map(_falla)
    d["deficit_SM_cal"] = np.maximum(rest.SM_min - d["SM_max_alcanzable_cal"].astype(float), 0.0)
    tol_min = (cfg.remolque.tol_min or 0.0) / MM
    d["deficit_tol_amarre_mm"] = np.maximum(tol_min - d["tol_en_SM_max_mm"].astype(float).fillna(0.0), 0.0)
    d = d.sort_values(["deficit_SM_cal", "deficit_tol_amarre_mm", "m_en_SM_max_g"], ascending=[True, True, False])
    cols = ["cand_id", "restriccion_que_falla", "deficit_SM_cal", "deficit_tol_amarre_mm", "SM_max_alcanzable_cal",
            "m_en_SM_max_g", "tol_en_SM_max_mm", "D_mm", "r_tip_rel_R", "D_ap_mm", "motivos"]
    return d[cols].head(n_max).reset_index(drop=True)


def fig_masa_vs_tope(tab: pd.DataFrame, cfg: ConfigOpt, ruta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .exportar import AZUL, NARANJA, AQUA, TINTA2, REJILLA
    colores = {4: AZUL, 6: NARANJA, 8: AQUA}
    fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
    for n, d in tab.groupby("n_aletas"):
        ax.plot(d["tope_r_tip_rel_R"], d["m_total_max_g"].astype(float) / 1000, marker="o", lw=1.8,
                color=colores.get(n, TINTA2), label=f"n = {n}" + ("" if n == cfg.aleta.n else " (diagnóstico)"))
    ax.axhline(cfg.m_max / G / 1000, color=TINTA2, ls=":", lw=1.0)
    ax.annotate(f"m_max = {cfg.m_max / G / 1000:g} kg", (tab["tope_r_tip_rel_R"].max(), cfg.m_max / G / 1000),
                textcoords="offset points", xytext=(-4, -12), ha="right", fontsize=8, color=TINTA2)
    t = cfg.restricciones.r_tip_rel_R_max
    if t is not None:
        ax.axvline(t, color=TINTA2, ls="--", lw=1.0)
        ax.annotate(f"tope real r_tip/R = {t:g}", (t, ax.get_ylim()[0]), textcoords="offset points",
                    xytext=(4, 6), fontsize=8, color=TINTA2)
    ax.set_xlabel("tope de r_tip / R")
    ax.set_ylabel("masa total máxima factible [kg]")
    ax.set_title("Masa factible frente al tope de las aletas", loc="left", fontsize=10)
    ax.grid(color=REJILLA, lw=0.6)
    ax.legend(frameon=False, loc="lower right")
    ruta.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta

