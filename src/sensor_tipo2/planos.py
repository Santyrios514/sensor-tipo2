"""Planos de los candidatos (spec v2 rev. 2 §8), sin JVM. Reemplazan a los `.ork`.

Cada plano es una hoja A3 apaisada a escala normalizada (1:x, la mayor que cabe), con:

- vista lateral en corte: perfil, capas de pared, tapones de plomo, electrónica, herraje, tubo de
  cola y aletas en verdadera magnitud; CG y CP con el brazo SM·D acotado;
- vista posterior: cuerpo, tubo de cola, aletas con su rotación y el círculo del D aparente;
- cotas en mm (L, L_n, L_c, L_t, L_tc, D, d_tc, D_ap, c_r, c_t, x_s, h, espesor de aleta, tapones de
  plomo, posición y largo de la electrónica), barra de escala y cajetín con masas, estabilidad,
  banderas y la línea de validación con OpenRocket (o "sin validar").

`plano` devuelve las cotas tal como se dibujaron (valor y texto) para poder contrastarlas con el
ranking (T11). Un candidato infactible se dibuja en su estado de SM máximo con el rótulo
INFACTIBLE y la restricción que falla.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Circle, Polygon, Rectangle  # noqa: E402

from .barrido import specs_de_fila  # noqa: E402
from .config import ConfigOpt  # noqa: E402
from .exportar import (AMARILLO, AQUA, AZUL, CAPAS, NARANJA, REJILLA, SUPERFICIE, TINTA, TINTA2,  # noqa: E402
                       Detalle, detalle)
from .perfiles import ESTACIONES  # noqa: E402

MM, G = 1e-3, 1e-3
HOJA_MM = (420.0, 297.0)  # A3 apaisado
MARCO_MM = 10.0
ESCALAS = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0)  # 1:x normalizadas
ROJO = "#c8322b"
FUENTE = 6.5  # pt de las cotas

# cota del plano → columna del ranking con el mismo valor (T11)
COL_RANKING = {"L": "L_mm", "L_n": "L_n_mm", "L_c": "L_c_mm", "L_t": "L_t_mm", "L_tc": "L_tc_mm", "D": "D_mm",
               "d_tc": "d_tc_mm", "D_ap": "D_ap_mm", "c_r": "c_r_mm", "c_t": "c_t_mm", "x_s": "x_s_mm", "h": "h_mm",
               "plomo_delantero": "ell_mm", "plomo_trasero": "ell_trasero_mm", "x_electronica": "x_electronica_mm"}


@dataclass
class Plano:
    rutas: list[Path]
    cotas: dict[str, tuple[float, str]] = field(default_factory=dict)  # nombre → (valor mm, texto dibujado)
    infactible: bool = False


# --------------------------------------------------------------------------- cotas


def _txt(v: float) -> str:
    """Cota en mm con hasta dos decimales (15.75, 46.13, 90)."""
    t = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if t in ("-0", "") else t


class _Cotas:
    """Dibuja cotas y guarda (valor, texto) de cada una."""

    def __init__(self, ax):
        self.ax = ax
        self.d: dict[str, tuple[float, str]] = {}

    def _reg(self, nombre, v, etiqueta):
        t = f"{etiqueta} {_txt(v)}" if etiqueta else _txt(v)
        self.d[nombre] = (float(v), _txt(v))
        return t

    def h(self, nombre, x0, x1, y, etiqueta="", y_ext=(None, None), arriba=True):
        """Cota horizontal entre x0 y x1 a la altura y; y_ext: desde dónde salen las líneas de referencia."""
        ax = self.ax
        t = self._reg(nombre, x1 - x0, etiqueta)
        for x, ye in zip((x0, x1), y_ext):
            if ye is not None:
                ax.plot([x, x], [ye, y + (1.5 if y > ye else -1.5)], color=TINTA2, lw=0.35)
        ax.annotate("", (x0, y), (x1, y), arrowprops=dict(arrowstyle="<|-|>", color=TINTA, lw=0.5,
                                                          mutation_scale=5, shrinkA=0, shrinkB=0))
        ax.text((x0 + x1) / 2, y + (0.8 if arriba else -0.8), t, ha="center", va="bottom" if arriba else "top",
                fontsize=FUENTE, color=TINTA, bbox=dict(fc="white", ec="none", pad=0.3))

    def v(self, nombre, x, y0, y1, etiqueta="", x_ext=(None, None), izquierda=True):
        ax = self.ax
        t = self._reg(nombre, y1 - y0, etiqueta)
        for y, xe in zip((y0, y1), x_ext):
            if xe is not None:
                ax.plot([xe, x + (-1.5 if x < xe else 1.5)], [y, y], color=TINTA2, lw=0.35)
        ax.annotate("", (x, y0), (x, y1), arrowprops=dict(arrowstyle="<|-|>", color=TINTA, lw=0.5,
                                                          mutation_scale=5, shrinkA=0, shrinkB=0))
        ax.text(x + (-0.8 if izquierda else 0.8), (y0 + y1) / 2, t, ha="right" if izquierda else "left",
                va="center", rotation=90, fontsize=FUENTE, color=TINTA, bbox=dict(fc="white", ec="none", pad=0.3))

    def valor(self, nombre, v):
        """Cota anotada como texto (no como línea): D_ap, espesor, posición de la electrónica."""
        self.d[nombre] = (float(v), _txt(v))
        return _txt(v)


# --------------------------------------------------------------------------- estado del llenado


@dataclass
class Estado:
    """Lo que se dibuja del lastre: tapones, electrónica, CG y CP (factible o en su SM máximo)."""

    ell: float
    ell2: float
    x_r0: float
    x_e: float
    x_CG: float
    x_CP: float
    m_total: float
    factible: bool
    fila: dict


def estado(cfg: ConfigOpt, det: Detalle) -> Estado | None:
    Ll = det.llenado
    if Ll is None or det.cp is None:
        return None
    mod = Ll.modelo
    if math.isfinite(Ll.ell):
        ell, ell2 = Ll.ell, Ll.ell2
        xcg, m = float(mod.x_CG_con_trasero(ell, ell2)), float(mod.m_con_trasero(ell, ell2))
        return Estado(ell, ell2, float(mod.x_r0(ell)), float(mod.x_e(ell)), xcg, det.cp.x_CP, m, Ll.ok, Ll.fila)
    # infactible: estado de SM máximo (solo tapón delantero)
    ell_geo = det.cu.lim.ell_geo
    if not ell_geo > 0:
        return None
    ells = np.linspace(0.0, ell_geo, cfg.numerico.n_ell)
    ell = float(ells[int(np.argmin(mod.x_CG(ells)))])
    return Estado(ell, 0.0, float(mod.x_r0(ell)), float(mod.x_e(ell)), float(mod.x_CG(ell)), det.cp.x_CP,
                  float(mod.m(ell)), False, Ll.fila)


# --------------------------------------------------------------------------- vistas


def _vista_lateral(ax, cfg: ConfigOpt, det: Detalle, est: Estado | None, cot: _Cotas, Rmax: float):
    cu, g = det.cu, det.g
    cav, p = cu.cav, cu.perfil
    x = cav.x / MM
    R, r_tc, L = p.R / MM, p.r_tc / MM, p.L / MM
    # aletas en verdadera magnitud (arriba y abajo)
    if g is not None:
        P = g.poligono / MM
        for s in (1, -1):
            ax.add_patch(Polygon(np.c_[P[:, 0], s * P[:, 1]], closed=True, fc=AMARILLO, alpha=0.55, ec=TINTA,
                                 lw=0.6, zorder=1))
    # pared por capas
    for i, s in enumerate(ESTACIONES):
        m = cav.estacion == i
        fr = cav.fronteras[s]
        for k, _ in enumerate(cfg.pared[s], start=1):
            for sg in (1, -1):
                ax.fill_between(x[m], sg * fr[k][m] / MM, sg * fr[k - 1][m] / MM, color=CAPAS[(k - 1) % len(CAPAS)],
                                lw=0, zorder=2)
    ax.plot(x, cav.r_e / MM, color=TINTA, lw=0.8, zorder=3)
    ax.plot(x, -cav.r_e / MM, color=TINTA, lw=0.8, zorder=3)
    ax.plot(x, cav.r_i / MM, color=TINTA2, lw=0.3, zorder=3)
    ax.plot(x, -cav.r_i / MM, color=TINTA2, lw=0.3, zorder=3)
    for xu in p.uniones:
        r = float(np.interp(xu, cav.x, cav.r_e)) / MM
        ax.plot([xu / MM] * 2, [-r, r], color=TINTA2, lw=0.4, ls=(0, (4, 2)), zorder=3)
    ax.plot([-6, L + 6], [0, 0], color=TINTA2, lw=0.4, ls=(0, (10, 2, 2, 2)), zorder=3)  # eje
    # lastre y electrónica
    Le = cfg.electronica.Le / MM
    if est is not None:
        tramos = [(cu.x_b0, cu.x_b0 + est.ell)] + ([(est.x_r0, est.x_r0 + est.ell2)] if est.ell2 > 0 else [])
        for t0, t1 in tramos:
            m = (cav.x >= t0) & (cav.x <= t1)
            ax.fill_between(x[m], -cav.r_i[m] / MM, cav.r_i[m] / MM, color=AZUL, alpha=0.55, lw=0, hatch="////",
                            ec="white", zorder=2.5)
        m = (cav.x >= est.x_e) & (cav.x <= est.x_e + cfg.electronica.Le)
        ax.fill_between(x[m], -cav.r_i[m] / MM, cav.r_i[m] / MM, color=AQUA, alpha=0.6, lw=0, zorder=2.5)
        ax.text((est.x_e / MM + Le / 2), 0, "electrónica", rotation=90, ha="center", va="center", fontsize=5.5,
                color=TINTA, zorder=6)
    for pm in cfg.puntuales:
        ax.add_patch(Rectangle((pm.x / MM - 2, -2), 4, 4, fc=TINTA, ec="none", zorder=6))
        ax.annotate(pm.nombre.replace("_", " "), (pm.x / MM, -2), xytext=(pm.x / MM, -R * 0.55),
                    fontsize=5.5, ha="center", color=TINTA, arrowprops=dict(arrowstyle="-", lw=0.4, color=TINTA2),
                    zorder=6, bbox=dict(fc="white", ec="none", pad=0.2))
    # CG, CP y brazo SM·D
    if est is not None:
        xcg, xcp = est.x_CG / MM, est.x_CP / MM
        ax.plot(xcg, 0, "o", ms=6, mfc="white", mec=TINTA, mew=1.0, zorder=7)
        ax.plot(xcg, 0, marker=(2, 0, 45), ms=6, color=TINTA, mew=0.8, zorder=7)
        ax.text(xcg, -3, f"CG {xcg:.1f}", ha="center", va="top", fontsize=FUENTE, zorder=7,
                bbox=dict(fc="white", ec="none", pad=0.3))
        ax.plot(xcp, 0, "D", ms=5, mfc=NARANJA, mec=TINTA, mew=0.6, zorder=7)
        ax.text(xcp, -3, f"CP {xcp:.1f}", ha="center", va="top", fontsize=FUENTE, zorder=7,
                bbox=dict(fc="white", ec="none", pad=0.3))
        SM = (est.x_CP - est.x_CG) / p.D
        cot.h("SM_D", xcg, xcp, 6.0, "SM·D =")
        ax.texts[-1].set_text(f"SM·D = {_txt(xcp - xcg)} ({SM:.2f} cal)")
    # cotas del cuerpo (abajo)
    re = lambda xx: float(np.interp(xx, cav.x / MM, cav.r_e / MM))  # noqa: E731
    y1, y2 = -(Rmax + 10), -(Rmax + 20)
    Ln, xt0, xtc0 = p.Ln / MM, p.x_t0 / MM, p.x_tc0 / MM
    cot.h("L_n", 0, Ln, y1, "L_n", (0, -re(Ln)))
    if p.L_c > 0:
        cot.h("L_c", Ln, xt0, y1, "L_c", (None, -re(xt0)))
    else:
        cot.valor("L_c", 0.0)
        ax.text(Ln, y1 + 0.8, "L_c = 0 (abombado)", ha="center", va="bottom", fontsize=FUENTE, color=TINTA2,
                bbox=dict(fc="white", ec="none", pad=0.3))
    cot.h("L_t", xt0, xtc0, y1, "L_t", (None, -re(xtc0)))
    cot.h("L_tc", xtc0, L, y1, "L_tc", (None, -re(L)))
    cot.h("L", 0, L, y2, "L", (y1, y1))
    # D y d_tc (verticales)
    cot.v("D", -14, -R, R, "D", (Ln, Ln))
    cot.v("d_tc", L + 24, -r_tc, r_tc, "d_tc", (L, L), izquierda=False)
    # aleta (arriba)
    if g is not None:
        xLE, cr, ct, xs = g.x_LE / MM, g.c_r / MM, g.c_t / MM, g.x_s / MM
        r_tip = g.r_tip / MM
        ya, yb = r_tip + 7, r_tip + 15
        if xs > 0:
            cot.h("x_s", xLE, xLE + xs, ya, "x_s", (r_tc, r_tip))
        else:
            cot.valor("x_s", 0.0)
        cot.h("c_t", xLE + xs, xLE + xs + ct, ya, "c_t", (None if xs > 0 else r_tip, r_tip))
        cot.h("c_r", xLE, xLE + cr, yb, "c_r", (r_tc, r_tc))
        cot.v("h", L + 10, r_tc, r_tip, "h", (L, xLE + xs + ct), izquierda=False)
    # lastre y electrónica (arriba del cuerpo)
    if est is not None:
        yp = R + 7
        b0 = cu.x_b0 / MM
        cot.h("plomo_delantero", b0, b0 + est.ell / MM, yp, "plomo", (re(b0), re(b0 + est.ell / MM)))
        xe = est.x_e / MM
        cot.h("electronica", xe, xe + Le, R + 15, "", (re(xe), re(xe + Le)))
        ax.texts[-1].set_text(f"electrónica {_txt(Le)} @ x = {cot.valor('x_electronica', xe)}")
        if est.ell2 > 0:
            a = est.x_r0 / MM
            cot.h("plomo_trasero", a, a + est.ell2 / MM, yp, "plomo", (re(a), re(a + est.ell2 / MM)))
        else:
            cot.valor("plomo_trasero", 0.0)
    # barra de escala
    yb_ = -(Rmax + 31)
    for i in range(5):
        ax.add_patch(Rectangle((i * 10, yb_), 10, 2.2, fc=TINTA if i % 2 == 0 else "white", ec=TINTA, lw=0.4))
    for v in (0, 50):
        ax.text(v, yb_ - 1, f"{v}", ha="center", va="top", fontsize=FUENTE)
    ax.text(52, yb_ + 1.1, "mm", ha="left", va="center", fontsize=FUENTE)


def _vista_posterior(ax, cfg: ConfigOpt, det: Detalle, cot: _Cotas, lim: float):
    p, g = det.cu.perfil, det.g
    R, r_tc = p.R / MM, p.r_tc / MM
    ax.add_patch(Circle((0, 0), R, fc=REJILLA, ec=TINTA, lw=0.8))
    ax.add_patch(Circle((0, 0), r_tc, fc=SUPERFICIE, ec=TINTA, lw=0.6))
    ax.plot([-lim * 0.9, lim * 0.9], [0, 0], color=TINTA2, lw=0.3, ls=(0, (10, 2, 2, 2)))
    ax.plot([0, 0], [-lim * 0.9, lim * 0.9], color=TINTA2, lw=0.3, ls=(0, (10, 2, 2, 2)))
    if g is None:
        return
    t, r_tip = g.params.t / MM, g.r_tip / MM
    for phi in g.params.rotacion + 2 * np.pi * np.arange(g.params.n) / g.params.n:
        u = np.array([np.sin(phi), np.cos(phi)])
        nrm = np.array([u[1], -u[0]])
        cc = np.array([r_tc * u + t / 2 * nrm, r_tip * u + t / 2 * nrm, r_tip * u - t / 2 * nrm, r_tc * u - t / 2 * nrm])
        ax.add_patch(Polygon(cc, closed=True, fc=AMARILLO, ec=TINTA, lw=0.5, zorder=3))
    D_ap = 2 * max(R, r_tip)
    ax.add_patch(Circle((0, 0), D_ap / 2, fc="none", ec=NARANJA, lw=0.8, ls=(0, (5, 3))))
    ax.text(0, D_ap / 2 + 2, f"D_ap = {cot.valor('D_ap', D_ap)}", ha="center", va="bottom", fontsize=FUENTE,
            color=NARANJA)
    ax.text(0, -lim + 2, f"{g.params.n} aletas a {math.degrees(g.params.rotacion):g}° · t = "
            f"{cot.valor('t_aleta', t)} · r_tip = {r_tip:.1f}", ha="center", va="bottom", fontsize=FUENTE)


# --------------------------------------------------------------------------- cajetín


def _filas_cajetin(cfg: ConfigOpt, fila: pd.Series, det: Detalle, est: Estado | None, puesto, val: dict | None,
                   escala: float, falla: str) -> list[tuple[str, str]]:
    g, cu = det.g, det.cu
    f = est.fila if est is not None else {}
    me = cfg.electronica.me / G
    m_punt = sum(pm.m for pm in cfg.puntuales) / G
    rho = cfg.lastre.rho_b
    mod = det.llenado.modelo if det.llenado is not None else None
    m_del = rho * float(mod.V_b(est.ell)) / G if est is not None else math.nan
    m_tras = rho * float(mod.V_tras(est.ell, est.ell2)) / G if est is not None and est.ell2 > 0 else 0.0
    m_tot = est.m_total / G if est is not None else math.nan
    SM = (est.x_CP - est.x_CG) / cu.spec.D if est is not None else math.nan
    banderas = [b for b in ("aletas_en_estela", "tubo_esbelto", "flutter_margen_bajo") if bool(fila.get(b, False))]
    if val is None:
        linea_val = "sin validar"
    else:
        linea_val = (f"OpenRocket: x_CP {val['x_CP_or_mm']:.2f} mm (Δ {val['dx_cp_mm']:+.2f}), SM {val['SM_or_cal']:.3f} "
                     f"(Δ {val['SM_or_cal'] - fila['SM_cal']:+.3f}), masa Δ {val['dif_masa_or_pct']:+.2f} % · "
                     + ("VALIDADO" if val.get("validado") else "NO VALIDADO"))
    estado_txt = "factible" if est is not None and est.factible else f"INFACTIBLE ({falla})"
    sufijo = "" if est is None or est.factible else " (en su SM máximo)"
    return [
        ("Candidato", str(fila["cand_id"])),
        ("Puesto en el ranking", "—" if puesto is None or pd.isna(puesto) else f"{int(puesto)} · {estado_txt}"),
        ("Masa total" + sufijo, f"{m_tot:.1f} g"),
        ("  casco / aletas", f"{cu.m_casco / G:.1f} g / {(g.masa / G if g is not None else math.nan):.1f} g"),
        ("  plomo delantero / trasero", f"{m_del:.1f} g / {m_tras:.1f} g"),
        ("  electrónica / herraje", f"{me:.1f} g / {m_punt:.1f} g"),
        ("x_CG / x_CP (desde la punta)", f"{(est.x_CG / MM if est else math.nan):.1f} / "
                                         f"{(est.x_CP / MM if est else math.nan):.1f} mm"),
        ("SM · C_Nα", f"{SM:.3f} cal · {fila.get('CN_alpha_total', math.nan):.3f}"),
        ("Tolerancia de amarre", f"{f.get('tol_amarre_mm', math.nan):.2f} mm"),
        ("θ_eq transición · AR aleta", f"{math.degrees(cu.theta_eq):.1f}° · {(g.AR if g is not None else math.nan):.2f}"),
        ("V_flutter", f"{fila.get('V_flutter_m_s', math.nan):.0f} m/s"),
        ("Restricción activa", str(fila.get("restriccion_activa", "")) if est is not None and est.factible else falla),
        ("Banderas", ", ".join(banderas) or "ninguna"),
        ("Validación", linea_val),
        ("Escala · unidades", f"1:{escala:g} en A3 · cotas en mm"),
    ]


def _cajetin(ax, filas: list[tuple[str, str]], titulo: str):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 1, 1, fc="white", ec=TINTA, lw=0.8))
    n = len(filas) + 1
    alto = 1 / n
    ax.add_patch(Rectangle((0, 1 - alto), 1, alto, fc=REJILLA, ec=TINTA, lw=0.6))
    ax.text(0.015, 1 - alto / 2, titulo, va="center", fontsize=7.5, weight="bold")
    for i, (k, v) in enumerate(filas, start=1):
        y = 1 - (i + 0.5) * alto
        ax.plot([0, 1], [1 - i * alto] * 2, color=TINTA2, lw=0.3)
        ax.text(0.015, y, k, va="center", fontsize=6.2, color=TINTA2)
        ax.text(0.33, y, v, va="center", fontsize=6.2, color=TINTA)
    ax.plot([0.32, 0.32], [0, 1 - alto], color=TINTA2, lw=0.3)


# --------------------------------------------------------------------------- plano


def _escala(ancho_mm: float, alto_mm: float, ancho_disp: float, alto_disp: float) -> float:
    for e in ESCALAS:
        if ancho_mm / e <= ancho_disp and alto_mm / e <= alto_disp:
            return e
    return ESCALAS[-1]


def plano(cfg: ConfigOpt, fila: pd.Series, ruta_base: Path, puesto=None, validacion: dict | None = None,
          falla: str = "", dpi: int = 300) -> Plano:
    """Dibuja el plano de `fila` (una fila del ranking) en ruta_base.png y ruta_base.pdf."""
    c, a = specs_de_fila(fila)
    det = detalle(cfg, c, a)
    est = estado(cfg, det)
    cu, g = det.cu, det.g
    if not cu.ok:
        raise ValueError(f"{fila['cand_id']}: cuerpo inválido ({';'.join(cu.motivos)})")
    p = cu.perfil
    L, R = p.L / MM, p.R / MM
    Rmax = max(R, g.r_tip / MM if g is not None else R)
    lat = (-34.0, L + 40.0, -(Rmax + 37.0), Rmax + 24.0)
    lim = Rmax + 14.0
    W, H = HOJA_MM
    ancho_disp = W - 2 * MARCO_MM - 26.0
    alto_disp = H - 2 * MARCO_MM - 120.0
    esc = _escala((lat[1] - lat[0]) + 2 * lim, max(lat[3] - lat[2], 2 * lim), ancho_disp, alto_disp)
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    fig.patch.set_facecolor("white")

    def ejes(x0, y0, w, h):
        return fig.add_axes([x0 / W, y0 / H, w / W, h / H])

    marco = ejes(0, 0, W, H)
    marco.set_xlim(0, W)
    marco.set_ylim(0, H)
    marco.axis("off")
    marco.add_patch(Rectangle((MARCO_MM, MARCO_MM), W - 2 * MARCO_MM, H - 2 * MARCO_MM, fc="none", ec=TINTA, lw=1.0))
    w_lat, h_lat = (lat[1] - lat[0]) / esc, (lat[3] - lat[2]) / esc
    y_vistas = H - MARCO_MM - 14.0 - max(h_lat, 2 * lim / esc)
    al = ejes(MARCO_MM + 8.0, y_vistas, w_lat, h_lat)
    al.set_xlim(lat[0], lat[1])
    al.set_ylim(lat[2], lat[3])
    al.set_aspect("equal")
    al.axis("off")
    cot = _Cotas(al)
    cot.valor("L_electronica", cfg.electronica.Le / MM)
    _vista_lateral(al, cfg, det, est, cot, Rmax)
    x_post = MARCO_MM + 8.0 + w_lat + 10.0
    ap = ejes(x_post, y_vistas + (h_lat - 2 * lim / esc) / 2, 2 * lim / esc, 2 * lim / esc)
    ap.set_xlim(-lim, lim)
    ap.set_ylim(-lim, lim)
    ap.set_aspect("equal")
    ap.axis("off")
    cot.ax = ap
    _vista_posterior(ap, cfg, det, cot, lim)
    y_tit = y_vistas + max(h_lat, 2 * lim / esc) + 3.0
    marco.text(MARCO_MM + 8.0, y_tit, "VISTA LATERAL EN CORTE (aletas en verdadera magnitud)", fontsize=8, weight="bold")
    marco.text(x_post, y_tit, "VISTA POSTERIOR (desde popa)", fontsize=8, weight="bold")
    infactible = est is None or not est.factible
    if infactible:
        marco.text(MARCO_MM + 8.0 + w_lat / 2, y_vistas + h_lat / 2, "INFACTIBLE", fontsize=46, color=ROJO,
                   alpha=0.35, ha="center", va="center", rotation=12, weight="bold", zorder=20)
        marco.text(MARCO_MM + 8.0, y_vistas - 4.0, f"INFACTIBLE · falla: {falla or 'restricción'}", fontsize=9,
                   color=ROJO, weight="bold", va="top")
    # tabla de cotas (las mismas que en las vistas)
    nombres = {"L": "L total", "L_n": "L_n nariz", "L_c": "L_c cuerpo cilíndrico", "L_t": "L_t transición",
               "L_tc": "L_tc tubo de cola", "D": "D cuerpo", "d_tc": "d_tc tubo de cola", "D_ap": "D_ap aparente",
               "c_r": "c_r cuerda de raíz", "c_t": "c_t cuerda de punta", "x_s": "x_s flecha", "h": "h envergadura",
               "t_aleta": "t espesor de aleta", "plomo_delantero": "tapón de plomo delantero",
               "plomo_trasero": "tapón de plomo trasero", "x_electronica": "x inicio de la electrónica",
               "L_electronica": "largo de la electrónica", "SM_D": "brazo SM·D (CG → CP)"}
    items = [(nombres[k], cot.d[k][1]) for k in nombres if k in cot.d]
    y_tab = y_vistas - 12.0
    marco.text(MARCO_MM + 6.0, y_tab, "COTAS [mm]", fontsize=7.5, weight="bold", va="top")
    por_col = math.ceil(len(items) / 3)
    for i, (k, v) in enumerate(items):
        xc = MARCO_MM + 6.0 + (i // por_col) * 62.0
        yc = y_tab - 6.0 - (i % por_col) * 4.6
        marco.text(xc, yc, k, fontsize=6.4, color=TINTA2, va="top")
        marco.text(xc + 52.0, yc, v, fontsize=6.4, color=TINTA, va="top", ha="right")
    # cajetín y notas
    filas = _filas_cajetin(cfg, fila, det, est, puesto, validacion, esc, falla)
    ancho_caj, alto_caj = 205.0, 92.0
    caj = ejes(W - MARCO_MM - ancho_caj, MARCO_MM, ancho_caj, alto_caj)
    _cajetin(caj, filas, "Sensor remolcado tipo 2 · DBF 2026-27 (UPB) · plano de optimización")
    pared = " + ".join(f"{cp.material} {cp.t / MM:.1f}" for cp in cfg.pared["cuerpo"])
    notas = [
        "NOTAS",
        "1. Cotas en mm; x desde la punta de la nariz. Escala válida al imprimir el PDF en A3 al 100 %.",
        f"2. Pared: {pared} mm (afuera → adentro). Plomo macizo ({cfg.lastre.rho_b:.0f} kg/m³), rayado.",
        f"3. Aletas: {g.params.n if g else '—'} × {g.params.material if g else ''}, espesor "
        f"{(g.params.t / MM if g else math.nan):.1f} mm, raíz al ras del extremo del tubo de cola.",
        "4. CG y CP del modelo propio (Barrowman); validación con OpenRocket en el cajetín.",
        "5. Diseño sin validación aerodinámica independiente: requiere CFD o ensayo antes de fabricarse.",
    ]
    for i, t in enumerate(notas):
        marco.text(MARCO_MM + 6.0, MARCO_MM + alto_caj - 4.0 - i * 5.0, t, fontsize=6.8 if i else 7.5,
                   weight="bold" if i == 0 else "normal", va="top")
    ruta_base.parent.mkdir(parents=True, exist_ok=True)
    # sin with_suffix: el cand_id tiene puntos (r1.2) que Path tomaría como extensión
    rutas = [ruta_base.parent / f"{ruta_base.name}.png", ruta_base.parent / f"{ruta_base.name}.pdf"]
    fig.savefig(rutas[0], dpi=dpi)
    fig.savefig(rutas[1])
    plt.close(fig)
    return Plano(rutas=rutas, cotas=cot.d, infactible=infactible)


# --------------------------------------------------------------------------- comparativo y perfiles


def comparativo(cfg: ConfigOpt, filas: list[tuple[str, pd.Series]], ruta: Path, dpi: int = 200) -> Path | None:
    """Siluetas laterales superpuestas a la misma escala, con masa y D aparente en la leyenda."""
    if not filas:
        return None
    fig, ax = plt.subplots(figsize=(14, 4.2), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for i, (etq, f) in enumerate(filas):
        c, a = specs_de_fila(f)
        det = detalle(cfg, c, a)
        cav = det.cu.cav
        col = cmap(i / max(len(filas) - 1, 1))
        x, r = cav.x / MM, cav.r_e / MM
        lab = f"{etq} · {f['m_total_g'] / 1000:.2f} kg · D_ap {f['D_ap_mm']:.1f} mm"
        ax.plot(np.r_[x, x[::-1]], np.r_[r, -r[::-1]], color=col, lw=1.2, label=lab)
        if det.g is not None:
            P = det.g.poligono / MM
            for s in (1, -1):
                ax.plot(np.r_[P[:, 0], P[0, 0]], s * np.r_[P[:, 1], P[0, 1]], color=col, lw=0.9)
    ax.axhline(0, color=TINTA2, lw=0.4, ls=(0, (10, 2, 2, 2)))
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm] (desde la punta)")
    ax.set_ylabel("r [mm]")
    ax.set_title("Siluetas de los mejores candidatos (misma escala; aletas en verdadera magnitud)", loc="left",
                 fontsize=10)
    ax.grid(color=REJILLA, lw=0.5)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ruta.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta, dpi=dpi)
    plt.close(fig)
    return ruta
