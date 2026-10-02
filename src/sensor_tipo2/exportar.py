"""Salidas: CSV, JSON, reconstrucción de un candidato con todo el detalle y figuras."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .aletas import GeomAleta, construir, envolvente  # noqa: E402
from .barrido import RANKING, completar, evaluar_cuerpo, specs_de_fila  # noqa: E402
from .config import AletaSpec, ConfigOpt, CuerpoSpec  # noqa: E402
from .geometria import Cuerpo, construir_cuerpo  # noqa: E402
from .lastre import Llenado, llenar  # noqa: E402
from .objetivo import tolerancias  # noqa: E402
from .perfiles import ESTACIONES  # noqa: E402
from .sustituto import IDENTIDAD, Calibracion, ResultadoCP, cp_sustituto  # noqa: E402

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
    for c in ("factible", "en_pareto", "verificado_or", "flutter_margen_bajo", "refinamiento", "ganador"):
        if c in df:
            df[c] = df[c].astype(str).str.lower().isin(["true", "1"])
    df["motivos"] = df["motivos"].fillna("")
    return df


def exportar_ranking(cfg: ConfigOpt, df: pd.DataFrame, dir_datos: Path):
    escribir_csv(df, dir_datos / "ranking.csv", RANKING)
    escribir_csv(df[df["en_pareto"].astype(bool)], dir_datos / "pareto.csv", RANKING)
    escribir_json(tolerancias(cfg), dir_datos / "tolerancias_implicitas.json")


# --------------------------------------------------------------------------- detalle de un candidato


@dataclass
class Detalle:
    cu: Cuerpo
    g: GeomAleta | None
    cp: ResultadoCP | None
    llenado: Llenado | None
    motivos: list[str]


def detalle(cfg: ConfigOpt, c: CuerpoSpec, a: AletaSpec, cal: Calibracion = IDENTIDAD,
            x_CP: float | None = None, CNa: float | None = None) -> Detalle:
    """Reconstruye el candidato con todo el detalle. x_CP y CNa fuerzan el CP (p. ej. el de
    OpenRocket) en lugar del sustituto o calibrado."""
    cu = construir_cuerpo(cfg, c)
    if not cu.ok:
        return Detalle(cu, None, None, None, cu.motivos)
    g, mot = construir(cfg, cu.perfil, a)
    if mot:
        return Detalle(cu, g, None, None, mot)
    cp, m = cp_sustituto(cu, g, cfg.vuelo.mach, cfg.numerico.n_franjas, cfg.restricciones.eps_CN)
    if cp is None:
        return Detalle(cu, g, None, None, [m])
    x = float(cal.aplicar(cp.x_CP)) if x_CP is None else x_CP
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
    out += [(p.nombre, p.m, p.x) for p in cfg.puntuales]
    return out


# --------------------------------------------------------------------------- figuras


def _estilo(ax):
    ax.grid(color=REJILLA, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_dibujo(cfg: ConfigOpt, det: Detalle, ruta: Path, titulo: str = "") -> Path | None:
    """Vista lateral (perfil, capas de pared, aletas, lastre, electrónica, CG y CP) y vista
    frontal del sensor guardado acostado con las aletas en diagonal."""
    if det.llenado is None or not math.isfinite(det.llenado.ell):
        return None
    cu, g, Ll = det.cu, det.g, det.llenado
    cav, p = cu.cav, cu.perfil
    x = cav.x / MM
    fig, (al, af) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [3.3, 1]},
                                 constrained_layout=True)
    # aletas proyectadas según su rotación de vuelo (detrás del cuerpo)
    P = g.poligono / MM
    for phi in g.params.rotacion + 2 * np.pi * np.arange(g.params.n) / g.params.n:
        c = np.cos(phi)
        if abs(c) < 1e-6:
            continue
        al.fill(P[:, 0], P[:, 1] * c, color=AMARILLO, alpha=0.45, lw=0, zorder=0.5)
        al.plot(np.r_[P[:, 0], P[0, 0]], np.r_[P[:, 1], P[0, 1]] * c, color=TINTA2, lw=0.8, zorder=0.5)
    # pared por capas
    etiquetas = set()
    for i, s in enumerate(ESTACIONES):
        m = cav.estacion == i
        fr = cav.fronteras[s]
        for k, capa in enumerate(cfg.pared[s], start=1):
            lab = capa.material if capa.material not in etiquetas else None
            etiquetas.add(capa.material)
            for sg in (1, -1):
                al.fill_between(x[m], sg * fr[k][m] / MM, sg * fr[k - 1][m] / MM, color=CAPAS[(k - 1) % len(CAPAS)],
                                lw=0, label=lab if sg == 1 else None)
    al.plot(x, cav.r_e / MM, color=TINTA, lw=1.1)
    al.plot(x, -cav.r_e / MM, color=TINTA, lw=1.1)
    mod = Ll.modelo
    tramos = [(cu.x_b0, cu.x_b0 + Ll.ell)]
    if Ll.ell2 > 0:
        a2 = float(mod.x_r0(Ll.ell))
        tramos.append((a2, a2 + Ll.ell2))
    for k, (t0, t1) in enumerate(tramos):
        m = (cav.x >= t0) & (cav.x <= t1)
        al.fill_between(x[m], -cav.r_i[m] / MM, cav.r_i[m] / MM, color=AZUL, alpha=0.7, lw=0,
                        label="lastre de plomo" if k == 0 else None)
    x_e = float(mod.x_e(Ll.ell))
    m = (cav.x >= x_e) & (cav.x <= x_e + cfg.electronica.Le)
    al.fill_between(x[m], -cav.r_i[m] / MM, cav.r_i[m] / MM, color=AQUA, alpha=0.55, lw=0, label="electrónica")
    f = Ll.fila
    xcg, xcp = f["x_CG_mm"], det.cp.x_CP / MM
    al.plot(xcg, 0, "o", ms=9, mfc=SUPERFICIE, mec=TINTA, mew=1.5, zorder=5)
    al.annotate(f"CG {xcg:.0f}", (xcg, 0), textcoords="offset points", xytext=(0, -16), ha="center", fontsize=8)
    al.plot(xcp, 0, "D", ms=7, mfc=NARANJA, mec=TINTA, mew=0.8, zorder=5)
    al.annotate(f"CP {xcp:.0f}", (xcp, 0), textcoords="offset points", xytext=(0, 9), ha="center", fontsize=8)
    ymax = max(p.R, g.r_tip) / MM * 1.25
    for xu in p.uniones:
        al.axvline(xu / MM, color=REJILLA, lw=0.8, zorder=0)
    al.annotate("", (0, -ymax * 0.95), (p.L / MM, -ymax * 0.95),
                arrowprops=dict(arrowstyle="<->", color=TINTA, lw=0.9, shrinkA=0, shrinkB=0))
    al.text(p.L / MM / 2, -ymax * 0.95, f"L = {p.L / MM:.0f} mm · D = {p.D / MM:.0f} · d_tc = {p.d_tc / MM:.1f}",
            ha="center", va="bottom", fontsize=8, bbox=dict(fc=SUPERFICIE, ec="none", pad=0.5))
    al.set_xlim(-8, p.L / MM + 8)
    al.set_ylim(-ymax * 1.08, ymax)
    al.set_aspect("equal")
    al.grid(False)
    al.set_xlabel("x [mm] (desde la punta)")
    al.set_title(f"Vista lateral · aletas a {math.degrees(g.params.rotacion):.0f}°", loc="left", fontsize=9)
    al.legend(frameon=False, loc="upper left", fontsize=7, ncol=3)

    # vista frontal: guardado acostado
    R, r_tip, t = p.R / MM, g.r_tip / MM, g.params.t / MM
    h_g, w_g = (v / MM for v in envolvente(g.r_tip, p.R, g.params.n, cfg.rot_guardado))
    for phi in cfg.rot_guardado + 2 * np.pi * np.arange(g.params.n) / g.params.n:
        u = np.array([np.sin(phi), np.cos(phi)])
        nrm = np.array([u[1], -u[0]])
        r0 = p.r_tc / MM
        cc = np.array([r0 * u + t / 2 * nrm, r_tip * u + t / 2 * nrm, r_tip * u - t / 2 * nrm, r0 * u - t / 2 * nrm])
        af.fill(cc[:, 0], cc[:, 1], color=AMARILLO, alpha=0.8, ec=TINTA2, lw=0.8)
    af.add_patch(plt.Circle((0, 0), R, fc=REJILLA, ec=TINTA, lw=1.1))
    af.add_patch(plt.Circle((0, 0), p.r_tc / MM, fc="none", ec=TINTA2, lw=0.6, ls=":"))
    D_ap = 2 * max(R, r_tip)
    af.add_patch(plt.Circle((0, 0), D_ap / 2, fc="none", ec=NARANJA, lw=1.0, ls="--"))
    af.add_patch(plt.Rectangle((-w_g / 2, -h_g / 2), w_g, h_g, fc="none", ec=AZUL, lw=1.2, ls="--"))
    lim = max(D_ap, h_g, w_g) / 2 * 1.3
    af.text(0, lim * 0.95, f"D aparente = {D_ap:.1f} mm", ha="center", va="top", fontsize=8, color=NARANJA)
    af.text(0, -lim * 0.95, f"D acostado = {max(h_g, w_g):.1f} mm", ha="center", va="bottom", fontsize=8, color=AZUL)
    af.set_xlim(-lim, lim)
    af.set_ylim(-lim, lim)
    af.set_aspect("equal")
    af.axis("off")
    af.set_title(f"Guardado acostado · aletas a {math.degrees(cfg.rot_guardado):.0f}°", loc="left", fontsize=9)
    fig.suptitle(titulo or (f"{cu.spec.id} · masa {f['m_total_g'] / 1000:.2f} kg (plomo {f['m_lastre_g'] / 1000:.2f} kg)"
                            f" · SM {f['SM_cal']:.2f} · tol. amarre {f['tol_amarre_mm']:.1f} mm"),
                 x=0.01, ha="left", fontsize=10)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    return ruta


def fig_ganador(cfg: ConfigOpt, fila: pd.Series, dir_fig: Path, cal: Calibracion = IDENTIDAD,
                x_CP: float | None = None, CNa: float | None = None) -> list[Path]:
    c, a = specs_de_fila(fila)
    det = detalle(cfg, c, a, cal, x_CP, CNa)
    r = fig_dibujo(cfg, det, dir_fig / "ganador_dibujo.png",
                   titulo=f"{fila['cand_id']} · CP de {'OpenRocket' if x_CP is not None else 'Barrowman interno'}")
    return [r] if r else []


def fig_pareto(df: pd.DataFrame, ruta: Path):
    p = df[df["en_pareto"].astype(bool)]
    if p.empty:
        return None
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    sc = axs[0].scatter(p["D_ap_mm"], p["k_efectivo"], c=p["SM_cal"], cmap="Blues", s=36, edgecolor=TINTA2, linewidth=0.5)
    axs[0].set_xlabel("D aparente [mm]")
    axs[0].set_ylabel("k efectivo")
    axs[0].set_title("Frente de Pareto (color = SM)", loc="left", fontsize=10)
    fig.colorbar(sc, ax=axs[0], label="SM [cal]")
    axs[1].scatter(p["D_ap_mm"], p["SM_cal"], s=30, color=AZUL, edgecolor="white", linewidth=1)
    axs[1].set_xlabel("D aparente [mm]")
    axs[1].set_ylabel("SM [cal]")
    axs[2].scatter(p["k_efectivo"], p["SM_cal"], s=30, color=NARANJA, edgecolor="white", linewidth=1)
    axs[2].set_xlabel("k efectivo")
    axs[2].set_ylabel("SM [cal]")
    for ax in axs:
        _estilo(ax)
    g = df[df["factible"].astype(bool)].head(1)
    if not g.empty:
        axs[0].plot(g["D_ap_mm"], g["k_efectivo"], marker="*", ms=16, color=AMARILLO, mec=TINTA, label="ganador")
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


VARIABLES = [("D_mm", "D [mm]"), ("L_n_rel_D", "L_n / D"), ("L_t_rel_D", "L_t / D"), ("k", "k = d_tc / D"),
             ("L_tc_mm", "L_tc [mm]"), ("mu_cr", "μ"), ("gamma_ct", "γ"), ("sigma_flecha", "σ"),
             ("r_tip_rel_R", "r_tip / R")]
_CUERPO = ("D_mm", "L_n_rel_D", "L_t_rel_D", "k", "L_tc_mm")


def sensibilidad(cfg: ConfigOpt, fila: pd.Series, cal: Calibracion = IDENTIDAD) -> pd.DataFrame:
    """J y SM variando una variable a la vez (valores de la malla) alrededor de `fila`."""
    from dataclasses import replace
    c0, a0 = specs_de_fila(fila)
    filas = []
    for var, _ in VARIABLES:
        for v in sorted(set(float(x) for x in cfg.malla[var]) | {float(fila[var])}):
            c = replace(c0, **{var: v}) if var in _CUERPO else c0
            a = replace(a0, **{var: v}) if var not in _CUERPO else a0
            for f in evaluar_cuerpo(cfg, c, [a], cal):
                f["variable"], f["valor"] = var, v
                filas.append(f)
    return completar(cfg, pd.DataFrame(filas))


def fig_sensibilidad(cfg: ConfigOpt, fila: pd.Series, ruta: Path, cal: Calibracion = IDENTIDAD):
    s = sensibilidad(cfg, fila, cal)
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


def fig_calibracion(ver: pd.DataFrame, cal: Calibracion, ruta: Path):
    ver = ver[np.isfinite(ver.get("x_CP_or_mm", pd.Series(dtype=float)).astype(float))]
    if ver.empty:
        return None
    fig, ax = plt.subplots(figsize=(5.6, 5), constrained_layout=True)
    for rol, color, etq in (("verificacion", AZUL, "mejores"), ("calibracion", NARANJA, "muestra estratificada")):
        d = ver[ver["rol"].str.contains(rol)]
        ax.scatter(d["x_CP_sust_mm"], d["x_CP_or_mm"], s=30, color=color, edgecolor="white", linewidth=1, label=etq)
    xs = np.linspace(ver["x_CP_sust_mm"].min(), ver["x_CP_sust_mm"].max(), 50)
    ax.plot(xs, cal.alpha / MM + cal.beta * xs, color=TINTA, lw=1.2,
            label=f"x_OR = {cal.alpha / MM:.2f} + {cal.beta:.4f} x_sust  (R² = {cal.R2:.4f})")
    ax.plot(xs, xs, color=TINTA2, lw=0.8, ls="--", label="identidad")
    ax.set_xlabel("x_CP sustituto [mm]")
    ax.set_ylabel("x_CP OpenRocket [mm]")
    ax.set_title(f"Calibración del CP · residuo máx. {cal.res_max / MM:.2f} mm", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    _estilo(ax)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)
    return ruta


def figuras_barrido(cfg: ConfigOpt, df: pd.DataFrame, dir_fig: Path, cal: Calibracion = IDENTIDAD) -> list[Path]:
    dir_fig.mkdir(parents=True, exist_ok=True)
    rutas = [fig_pareto(df, dir_fig / "pareto.png"), fig_factibilidad(df, dir_fig / "factibilidad.png")]
    fac = df[df["factible"].astype(bool)]
    if not fac.empty:
        g = fac.iloc[0]
        rutas += fig_ganador(cfg, g, dir_fig, cal)
        rutas.append(fig_sensibilidad(cfg, g, dir_fig / "sensibilidad.png", cal))
    return [r for r in rutas if r is not None]
