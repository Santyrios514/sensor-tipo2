"""Salidas: CSV, JSON, reconstrucción de un candidato con todo el detalle y figuras."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from .aletas import GeomAleta, construir  # noqa: E402
from .barrido import RANKING, completar, evaluar_cuerpo, specs_de_fila  # noqa: E402
from .config import AletaSpec, ConfigOpt, CuerpoSpec  # noqa: E402
from .geometria import Cuerpo, construir_cuerpo  # noqa: E402
from .lastre import Llenado, llenar  # noqa: E402
from .objetivo import tolerancias  # noqa: E402
from .sustituto import ResultadoCP, cp_sustituto  # noqa: E402

MM = 1e-3
AZUL, NARANJA, AQUA, AMARILLO = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
TINTA, TINTA2, REJILLA, SUPERFICIE = "#1f1f1e", "#5f5e58", "#e4e3df", "#fcfcfb"
CAPAS = ["#eda100", "#f5b89c", "#a8dfc7", "#e87ba4"]


def escribir_csv(df: pd.DataFrame, ruta: Path, columnas: list[str] | None = None):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    if columnas is not None:
        df = df.reindex(columns=list(dict.fromkeys(columnas + [c for c in df.columns if c not in columnas])))
    df.to_csv(ruta, index=False, float_format="%.6g")


def escribir_json(obj: dict, ruta: Path):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    ruta.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=float), encoding="utf-8")


def leer_ranking(ruta: Path) -> pd.DataFrame:
    df = pd.read_csv(ruta, keep_default_na=True, low_memory=False)
    for c in ("factible", "en_pareto", "validado_or", "flutter_margen_bajo", "refinamiento", "ganador",
              "tubo_esbelto", "aletas_en_estela"):
        if c in df:
            df[c] = df[c].astype(str).str.lower().isin(["true", "1"])
    df["motivos"] = df["motivos"].fillna("")
    return df


def exportar_ranking(cfg: ConfigOpt, df: pd.DataFrame, dir_datos: Path):
    escribir_csv(df, dir_datos / "ranking.csv", RANKING)
    escribir_csv(df[df["en_pareto"].astype(bool)], dir_datos / "pareto.csv", RANKING)
    escribir_csv(frontera(df), dir_datos / "frontera_masa_D.csv")
    escribir_json(tolerancias(cfg), dir_datos / "tolerancias_implicitas.json")


# --------------------------------------------------------------------------- detalle de un candidato


@dataclass
class Detalle:
    cu: Cuerpo
    g: GeomAleta | None
    cp: ResultadoCP | None
    llenado: Llenado | None
    motivos: list[str]


def detalle(cfg: ConfigOpt, c: CuerpoSpec, a: AletaSpec, x_CP: float | None = None,
            CNa: float | None = None) -> Detalle:
    """Reconstruye el candidato con todo el detalle. x_CP y CNa fuerzan el CP en lugar del del
    modelo propio (diagnóstico)."""
    cu = construir_cuerpo(cfg, c)
    if not cu.ok:
        return Detalle(cu, None, None, None, cu.motivos)
    g, mot = construir(cfg, cu.perfil, a)
    if mot:
        return Detalle(cu, g, None, None, mot)
    cp, m = cp_sustituto(cu, g, cfg.vuelo.mach, cfg.numerico.n_franjas, cfg.restricciones.eps_CN)
    if cp is None:
        return Detalle(cu, g, None, None, [m])
    x = cp.x_CP if x_CP is None else x_CP
    Ll = llenar(cfg, cu, g.masa, g.x_cg, x, cp.CNa if CNa is None else CNa)
    return Detalle(cu, g, cp.con_x(x), Ll, Ll.motivos)


def masas_del_llenado(cfg: ConfigOpt, det: Detalle) -> list[tuple[str, float, float]]:
    """(nombre, m, x_CG) del lastre delantero y trasero, la electrónica y las masas puntuales."""
    Ll = det.llenado
    mod = Ll.modelo
    out = []
    if Ll.ell > 0:
        out.append(("lastre_delantero", float(mod.m_b(Ll.ell)), float(mod.xbar_b(Ll.ell))))
    if Ll.ell2 > 0:
        a = float(mod.x_r0(Ll.ell))
        V2 = float(mod.fll * mod.cav.vol(a, a + Ll.ell2))
        M2 = float(mod.fll * mod.cav.mom(a, a + Ll.ell2))
        out.append(("lastre_trasero", mod.rho_b * V2, M2 / V2))
    out.append(("electronica", cfg.electronica.me, float(cfg.electronica.x_cg(mod.x_e(Ll.ell)))))
    x_cg = float(mod.x_CG_con_trasero(Ll.ell, Ll.ell2))  # el herraje (x_mm: cg) va en el CG final
    out += [(p.nombre, p.m, x_cg if p.en_cg else p.x) for p in cfg.puntuales]
    return out


# --------------------------------------------------------------------------- figuras


def _estilo(ax):
    ax.grid(color=REJILLA, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def perfil_candidato(cfg: ConfigOpt, fila: pd.Series) -> pd.DataFrame:
    """Perfil de un candidato para CAD/CFD: x, r_e y r_i del casco y el polígono de una aleta (mm)."""
    c, a = specs_de_fila(fila)
    cu = construir_cuerpo(cfg, c)
    g, _ = construir(cfg, cu.perfil, a)
    cav = cu.cav
    paso = max(1, int(round(0.5e-3 / cav.dx)))  # cada 0.5 mm
    casco = pd.DataFrame({"elemento": "casco", "x_mm": cav.x[::paso] / MM, "r_e_mm": cav.r_e[::paso] / MM,
                          "r_i_mm": cav.r_i[::paso] / MM})
    P = g.poligono / MM
    aleta = pd.DataFrame({"elemento": "aleta", "x_mm": P[:, 0], "r_e_mm": P[:, 1]})
    return pd.concat([casco, aleta], ignore_index=True)


def fig_pareto(df: pd.DataFrame, ruta: Path, col: str = "D_acostado_mm"):
    """Frente de Pareto en (diámetro del criterio 2, k efectivo, SM)."""
    p = df[df["en_pareto"].astype(bool)]
    if p.empty:
        return None
    etq = "D acostado [mm]" if col == "D_acostado_mm" else "D aparente [mm]"
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    sc = axs[0].scatter(p[col], p["k_efectivo"], c=p["SM_cal"], cmap="Blues", s=36, edgecolor=TINTA2, linewidth=0.5)
    axs[0].set_xlabel(etq)
    axs[0].set_ylabel("k efectivo")
    axs[0].set_title("Frente de Pareto (color = SM)", loc="left", fontsize=10)
    fig.colorbar(sc, ax=axs[0], label="SM [cal]")
    axs[1].scatter(p[col], p["SM_cal"], s=30, color=AZUL, edgecolor="white", linewidth=1)
    axs[1].set_xlabel(etq)
    axs[1].set_ylabel("SM [cal]")
    axs[2].scatter(p["k_efectivo"], p["SM_cal"], s=30, color=NARANJA, edgecolor="white", linewidth=1)
    axs[2].set_xlabel("k efectivo")
    axs[2].set_ylabel("SM [cal]")
    for ax in axs:
        _estilo(ax)
    g = df[df["factible"].astype(bool)].head(1)
    if not g.empty:
        axs[0].plot(g[col], g["k_efectivo"], marker="*", ms=16, color=AMARILLO, mec=TINTA, label="ganador")
        axs[0].legend(frameon=False, loc="upper right")
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


def fig_factibilidad(df: pd.DataFrame, ruta: Path):
    m = df.loc[~df["factible"].astype(bool), "motivos"].fillna("").str.split(";").explode()
    m = m[m != ""].value_counts().sort_values()
    if m.empty:
        return None
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(m) + 1.4), constrained_layout=True)
    ax.barh(m.index, m.values, color=AZUL, height=0.6)
    for y, v in enumerate(m.values):
        ax.text(v, y, f" {v:,}", va="center", fontsize=8, color=TINTA)
    ax.set_xlabel("candidatos descartados (un candidato puede tener varios motivos)")
    ax.set_title(f"Motivos de descarte · {int((~df['factible'].astype(bool)).sum()):,} de {len(df):,} candidatos",
                 loc="left", fontsize=10)
    _estilo(ax)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


_CUERPO = ("D_mm", "L_n_rel_D", "f_cil", "k", "d_tc_mm", "L_tc_mm")


def variables(cfg: ConfigOpt) -> list[tuple[str, str]]:
    tubo = ("k", "k = d_tc / D") if cfg.var_tubo == "k" else ("d_tc_mm", "d_tc [mm]")
    return [("D_mm", "D [mm]"), ("L_n_rel_D", "L_n / D"), ("f_cil", "f_c"), tubo, ("L_tc_mm", "L_tc [mm]"),
            ("mu_cr", "μ"), ("gamma_ct", "γ"), ("sigma_flecha", "σ"), ("r_tip_rel_R", "r_tip / R")]


def _variar(cfg: ConfigOpt, c0: CuerpoSpec, var: str, v: float) -> CuerpoSpec:
    """El cuerpo con una variable cambiada; con la malla en k, el tubo conserva su k al cambiar D."""
    from dataclasses import replace
    if var == "k":
        return replace(c0, d_tc_mm=round(v * c0.D_mm, 9))
    if var == "D_mm" and cfg.var_tubo == "k":
        return replace(c0, D_mm=v, d_tc_mm=round(c0.k * v, 9))
    return replace(c0, **{var: v})


def sensibilidad(cfg: ConfigOpt, fila: pd.Series) -> pd.DataFrame:
    """J y SM variando una variable a la vez (valores de la malla) alrededor de `fila`."""
    from dataclasses import replace
    c0, a0 = specs_de_fila(fila)
    filas = []
    for var, _ in variables(cfg):
        for v in sorted(set(float(x) for x in cfg.malla[var]) | {float(fila[var])}):
            c = _variar(cfg, c0, var, v) if var in _CUERPO else c0
            a = replace(a0, **{var: v}) if var not in _CUERPO else a0
            for f in evaluar_cuerpo(cfg, c, [a]):
                f["variable"], f["valor"] = var, v
                filas.append(f)
    return completar(cfg, pd.DataFrame(filas))


def fig_sensibilidad(cfg: ConfigOpt, fila: pd.Series, ruta: Path):
    s = sensibilidad(cfg, fila)
    VARIABLES = variables(cfg)
    fig, axs = plt.subplots(2, len(VARIABLES), figsize=(2.0 * len(VARIABLES), 5.2), sharey="row",
                            constrained_layout=True)
    J0 = float(fila["J"])
    for i, (var, etq) in enumerate(VARIABLES):
        d = s[s["variable"] == var].sort_values("valor")
        ok = d["factible"].astype(bool)
        for ax, col, color in ((axs[0, i], "J", AZUL), (axs[1, i], "SM_cal", NARANJA)):
            y = d[col] - J0 if col == "J" else d[col]
            ax.plot(d["valor"], y, color=color, lw=1.5, zorder=1)
            ax.scatter(d.loc[ok, "valor"], y[ok], s=22, color=color, zorder=2)
            ax.scatter(d.loc[~ok, "valor"], y[~ok], s=22, facecolor="white", edgecolor=color, zorder=2)
            ax.axvline(float(fila[var]), color=TINTA2, lw=0.8, ls=":")
            _estilo(ax)
        axs[1, i].set_xlabel(etq)
    axs[0, 0].set_ylabel("J − J_ganador")
    axs[0, 0].set_yscale("symlog", linthresh=1.0)
    axs[1, 0].set_ylabel("SM [cal]")
    for ax in axs[1]:
        ax.axhspan(cfg.restricciones.SM_min, cfg.restricciones.SM_max, color=REJILLA, alpha=0.5, lw=0, zorder=0)
    fig.suptitle("Sensibilidad alrededor del ganador (punto hueco = infactible)", x=0.01, ha="left", fontsize=10)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


def frontera(df: pd.DataFrame) -> pd.DataFrame:
    """Por cada D, el mejor candidato factible (menor J, que es la mayor masa salvo empates): la masa
    máxima que admite cada diámetro. Sirve para elegir qué tan esbelto hacer el sensor."""
    fac = df[df["factible"].astype(bool)]
    if fac.empty:
        return pd.DataFrame()
    cols = ["D_mm", "D_acostado_mm", "D_ap_mm", "m_total_g", "m_lastre_g", "SM_cal", "tol_amarre_mm",
            "restriccion_activa", "r_tip_rel_R", "L_n_mm", "L_t_mm", "d_tc_mm", "L_tc_mm", "cola_forma",
            "theta_eq_deg", "cand_id"]
    idx = fac.sort_values(["D_mm", "m_total_g", "J"], ascending=[True, False, True]).groupby("D_mm").head(1).index
    return fac.loc[idx, cols].reset_index(drop=True)


def fig_frontera(fr: pd.DataFrame, m_max_g: float, ruta: Path):
    if fr.empty:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 4), constrained_layout=True)
    ax.plot(fr["D_mm"], fr["m_total_g"] / 1000, color=AZUL, marker="o", lw=1.8)
    for _, r in fr.iterrows():
        ax.annotate(r["restriccion_activa"], (r["D_mm"], r["m_total_g"] / 1000), textcoords="offset points",
                    xytext=(0, 7), ha="center", fontsize=7, color=TINTA2)
    ax.axhline(m_max_g / 1000, color=NARANJA, lw=1.0, ls="--", label=f"m_max = {m_max_g / 1000:g} kg")
    ax.set_xlabel("D del cuerpo [mm]")
    ax.set_ylabel("masa total máxima factible [kg]")
    ax.set_title("Frontera masa–diámetro (etiqueta = restricción activa)", loc="left", fontsize=10)
    ax.legend(frameon=False, loc="upper left")
    _estilo(ax)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


def figuras_barrido(cfg: ConfigOpt, df: pd.DataFrame, dir_fig: Path) -> list[Path]:
    """Figuras de análisis (pareto, factibilidad, frontera masa–D, sensibilidad). Los dibujos de los
    candidatos son los planos (planos.py, script 04)."""
    dir_fig.mkdir(parents=True, exist_ok=True)
    rutas = [fig_pareto(df, dir_fig / "pareto.png", cfg.objetivo.col_diametro),
             fig_factibilidad(df, dir_fig / "factibilidad.png"),
             fig_frontera(frontera(df), cfg.m_max * 1e3, dir_fig / "frontera_masa_D.png")]
    fac = df[df["factible"].astype(bool)]
    if not fac.empty:
        rutas.append(fig_sensibilidad(cfg, fac.iloc[0], dir_fig / "sensibilidad.png"))
    return [r for r in rutas if r is not None]
