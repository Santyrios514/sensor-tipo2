"""Llenado de lastre de plomo hasta m_max con SM ≥ SM_min (el criterio de dbf-sensor, sin cambios).

Con la electrónica detrás del tapón delantero (x_e = x_b0 + ℓ + holgura):

    x_CG(ℓ) = (M_0 + ρ_b Φ1(ℓ) + m_e x̄_e(ℓ)) / (m_0 + ρ_b ∀(ℓ) + m_e),     SM = (x_CP − x_CG) / D

Orden de llenado (el de `sensor_lastre._analizar_relleno`):

1. tapón delantero desde x_b0 hasta ℓ_SM,max = sup{ℓ ∈ [0, ℓ_geo] : SM(ℓ) ≥ SM_min};
2. si llegó a ℓ_geo, tapón trasero detrás de la electrónica (transición y tubo de cola) hasta
   SM = SM_min o cavidad llena (SM(ℓ₂) es monótona: `brentq`);
3. si la masa total supera m_max, el lastre se recorta a m_max − m_vacío y se ubica igual
   (delantero primero); si así el SM no llega a SM_min → `SM_inalcanzable_con_masa_max`;
4. si la tolerancia del amarre queda por debajo de tol_min, se retrocede por el camino de
   llenado (primero el trasero) hasta el punto de mayor masa que cumple ambas condiciones:

    tol_amarre = α_max q S_ref C_Nα (x_CP − x_CG) / (m g)

El SM_max no entra al llenado: si el SM final lo supera, el candidato es infactible
(`SM_sobre_max`); lo resuelven aletas más pequeñas.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from . import G0
from .config import ConfigOpt, Electronica
from .geometria import Cavidad, Cuerpo

MM, G = 1e-3, 1e-3

INFACTIBLES = ("inviable_geo", "SM_inalcanzable", "SM_inalcanzable_por_geometria",
               "SM_inalcanzable_con_masa_max", "tol_amarre_inalcanzable", "trim_no_lineal")
RESTRICCION = {"masa_max_sensor": "masa", "SM": "SM_min", "geometria": "volumen", "tol_amarre": "tol_amarre"}


@dataclass
class ModeloLastre:
    cav: Cavidad
    m0: float  # vacío sin electrónica ni lastre (casco, aletas, masas puntuales)
    M0: float
    el: Electronica
    x_b0: float
    rho_b: float
    fll: float
    margen_popa: float

    @property
    def ell_cavidad(self) -> float:
        return max(self.cav.x_fin - self.x_b0, 0.0)

    def V_b(self, ell):
        return self.fll * self.cav.vol(self.x_b0, self.x_b0 + np.asarray(ell, dtype=float))

    def Phi1(self, ell):
        return self.fll * self.cav.mom(self.x_b0, self.x_b0 + np.asarray(ell, dtype=float))

    def m_b(self, ell):
        return self.rho_b * self.V_b(ell)

    def xbar_b(self, ell):
        return self.Phi1(ell) / self.V_b(ell)

    def x_e(self, ell):
        return self.x_b0 + np.asarray(ell, dtype=float) + self.el.holgura

    def m(self, ell):
        return self.m0 + self.m_b(ell) + self.el.me

    def x_CG(self, ell):
        num = self.M0 + self.rho_b * self.Phi1(ell) + self.el.me * self.el.x_cg(self.x_e(ell))
        return num / self.m(ell)

    # --- tapón trasero: desde x_r0(ℓ), detrás de la electrónica, hacia popa
    def x_r0(self, ell):
        return self.x_e(ell) + self.el.Le + self.el.holgura

    def ell2_max(self, ell) -> float:
        return max(self.cav.x_fin - self.margen_popa - float(self.x_r0(ell)), 0.0)

    def V_tras(self, ell, ell2):
        a = self.x_r0(ell)
        return self.fll * self.cav.vol(a, a + np.asarray(ell2, dtype=float))

    def m_con_trasero(self, ell, ell2):
        return self.m(ell) + self.rho_b * self.V_tras(ell, ell2)

    def x_CG_con_trasero(self, ell, ell2):
        a = self.x_r0(ell)
        M_tras = self.rho_b * self.fll * self.cav.mom(a, a + np.asarray(ell2, dtype=float))
        return (self.x_CG(ell) * self.m(ell) + M_tras) / self.m_con_trasero(ell, ell2)

    def ell2_por_SM(self, ell, x_CP: float, D: float, SM_min: float, tol: float) -> float:
        l2max = self.ell2_max(ell)

        def g(l2):
            return (x_CP - float(self.x_CG_con_trasero(ell, l2))) / D - SM_min

        if l2max <= 0 or g(0.0) < 0:
            return 0.0
        if g(l2max) >= 0:
            return l2max
        return brentq(g, 0.0, l2max, xtol=tol)

    def ell_de_masa(self, m_obj: float, tol: float) -> float:
        """ρ_b ∀(ℓ) = m_obj en [0, ℓ_cavidad]; NaN si no cabe."""
        ell_c = self.ell_cavidad
        if m_obj <= 0:
            return 0.0
        if self.m_b(ell_c) < m_obj:
            return math.nan
        return brentq(lambda l: float(self.m_b(l)) - m_obj, 0.0, ell_c, xtol=tol)


def tolerancia_amarre(m: float, x_CG: float, x_CP: float, q: float, S_ref: float, CNa: float,
                      alpha_max: float) -> float:
    brazo = x_CP - x_CG
    return alpha_max * q * S_ref * CNa * brazo / (m * G0) if brazo > 0 else math.nan


def _ell_SM_max(mod: ModeloLastre, ell_geo: float, x_max: float, n: int, tol: float) -> float:
    """sup{ℓ ∈ [0, ℓ_geo] : x_CG(ℓ) ≤ x_max}; NaN si el conjunto es vacío."""
    ells = np.linspace(0.0, ell_geo, n)
    h = mod.x_CG(ells) - x_max
    ok = np.nonzero(h <= 0)[0]
    if not ok.size:
        return math.nan
    i = int(ok[-1])
    if i == n - 1:
        return float(ells[-1])
    return brentq(lambda l: float(mod.x_CG(l)) - x_max, ells[i], ells[i + 1], xtol=tol)


def _n_valles(y: np.ndarray, tol: float) -> int:
    d = np.diff(y)
    s = np.sign(np.where(np.abs(d) <= tol, 0.0, d))
    s = s[s != 0]
    if s.size == 0:
        return 1
    return int(np.sum((s[:-1] < 0) & (s[1:] > 0))) + int(s[-1] < 0)


def _ubicar_masa(mod: ModeloLastre, m_obj: float, ell_geo: float, trasero: bool, tol: float):
    """(ℓ, ℓ₂) para m_obj de lastre: delantero desde x_b0 y, si no alcanza con ℓ_geo, el resto
    detrás de la electrónica. ℓ = NaN si no cabe."""
    ell, ell2 = mod.ell_de_masa(m_obj, tol), 0.0
    m_del_max = float(mod.m_b(ell_geo)) if ell_geo > 0 else 0.0
    if trasero and ell_geo > 0 and m_obj > m_del_max:
        resto = m_obj - m_del_max
        l2max = mod.ell2_max(ell_geo)
        if mod.rho_b * float(mod.V_tras(ell_geo, l2max)) >= resto:
            ell = ell_geo
            ell2 = brentq(lambda l2: mod.rho_b * float(mod.V_tras(ell_geo, l2)) - resto, 0.0, l2max, xtol=tol)
        else:
            ell = math.nan
    return ell, ell2


def _recortar_por_amarre(mod, ell1, ell2, x_CP, D, SM_min, q, S_ref, CNa, alpha_max, tol_min, n_ell, tol):
    """Punto de mayor masa sobre el camino de llenado (delantero 0 → ℓ1, luego trasero 0 → ℓ2)
    con tol_amarre ≥ tol_min y SM ≥ SM_min. (NaN, NaN) si ninguno cumple."""
    def punto(s):
        return (min(s, ell1), max(s - ell1, 0.0))

    def holgura(s):
        """min(holgura relativa de tol_amarre, holgura de SM) en el punto s del camino (vectorial)."""
        s = np.asarray(s, dtype=float)
        a, b = np.minimum(s, ell1), np.maximum(s - ell1, 0.0)
        m = mod.m_con_trasero(a, b)
        xcg = mod.x_CG_con_trasero(a, b)
        brazo = x_CP - xcg
        t = alpha_max * q * S_ref * CNa * brazo / (m * G0)
        h_tol = np.where(brazo > 0, (t - tol_min) / tol_min, -1.0)
        return np.minimum(h_tol, brazo / D - SM_min + 1e-6)

    s_tot = ell1 + ell2
    if holgura(s_tot) >= 0:
        return ell1, ell2
    n = max(n_ell // 10, 40)
    ss = np.linspace(0.0, s_tot, n)
    h = holgura(ss)
    ok = np.nonzero(h >= 0)[0]
    if not ok.size:
        return math.nan, math.nan
    i = int(ok[-1])
    s = ss[i]
    if i + 1 < n:
        s = brentq(lambda v: float(holgura(v)), ss[i], ss[i + 1], xtol=tol)
        if holgura(s) < 0:
            s = max(s - tol, ss[i])
    return punto(s)


@dataclass
class Llenado:
    modelo: ModeloLastre
    ell: float = math.nan  # tapón delantero [m]
    ell2: float = 0.0  # tapón trasero [m]
    limitante: str = ""
    banderas: list[str] = field(default_factory=list)
    fila: dict = field(default_factory=dict)

    @property
    def motivos(self) -> list[str]:
        return [b for b in self.banderas if b in INFACTIBLES or b == "SM_sobre_max"]

    @property
    def ok(self) -> bool:
        return math.isfinite(self.ell) and not self.motivos

    @property
    def restriccion_activa(self) -> str:
        return RESTRICCION.get(self.limitante, "ninguna")


def llenar(cfg: ConfigOpt, cu: Cuerpo, m_aletas: float, x_aletas: float, x_CP: float, CNa: float) -> Llenado:
    """Llenado de máxima masa con SM ≥ SM_min (pasos 1–4 del módulo)."""
    rest, nu, la, vu = cfg.restricciones, cfg.numerico, cfg.lastre, cfg.vuelo
    D = cu.spec.D
    S_ref = math.pi * D**2 / 4
    SM_min = rest.SM_min
    mod = ModeloLastre(cav=cu.cav, m0=cu.m_casco + m_aletas + cu.m_puntuales,
                       M0=cu.M_casco + m_aletas * x_aletas + cu.M_puntuales, el=cfg.electronica,
                       x_b0=cu.x_b0, rho_b=la.rho_b, fll=la.fll, margen_popa=la.margen_popa)
    Ll = Llenado(modelo=mod)
    ell_geo = cu.lim.ell_geo
    if not ell_geo > 0:
        Ll.banderas.append("inviable_geo")
        return Ll

    def SM(xcg):
        return (x_CP - xcg) / D

    x_max = x_CP - SM_min * D
    ells = np.linspace(0.0, ell_geo, nu.n_ell)
    xcg = mod.x_CG(ells)
    i = int(np.argmin(xcg))
    lo, hi = ells[max(i - 1, 0)], ells[min(i + 1, nu.n_ell - 1)]
    ell_star = float(ells[i])
    if hi > lo:
        r = minimize_scalar(lambda l: float(mod.x_CG(l)), bounds=(lo, hi), method="bounded",
                            options={"xatol": nu.tol})
        if r.fun <= xcg[i]:
            ell_star = float(r.x)
    SM_inf = (x_CP - float(mod.xbar_b(ell_star))) / D if ell_star > 0 else math.nan
    Ll.fila.update({"SM_max_alcanzable_cal": SM(float(mod.x_CG(ell_star))), "SM_inf_cal": SM_inf})
    if _n_valles(xcg, nu.tol * 1e-3) > 1:
        Ll.banderas.append("no_unimodal")
    if math.isfinite(SM_inf) and SM_inf < SM_min:
        Ll.banderas += ["SM_inalcanzable_por_geometria", "SM_inalcanzable"]
        return Ll
    if math.isfinite(SM_inf) and SM_inf - SM_min < cfg.umbral_margen_bajo:
        Ll.banderas.append("margen_SM_bajo")

    ell_u = _ell_SM_max(mod, ell_geo, x_max, nu.n_ell, nu.tol)
    if not math.isfinite(ell_u):
        Ll.banderas.append("SM_inalcanzable")
        return Ll
    # 1–2: delantero; trasero solo con el delantero lleno
    ell2 = 0.0
    lleno = ell_u >= ell_geo - nu.tol
    if la.trasero and lleno:
        ell2 = mod.ell2_por_SM(ell_u, x_CP, D, SM_min, nu.tol)
    if not lleno:
        limitante = "SM"
    elif la.trasero:
        limitante = "SM" if ell2 < mod.ell2_max(ell_u) - nu.tol else "geometria"
    else:
        limitante = "geometria"
    # 3: tope de masa del sensor
    if float(mod.m_con_trasero(ell_u, ell2)) > cfg.m_max:
        ell, ell2 = _ubicar_masa(mod, cfg.m_max - float(mod.m(0.0)), ell_geo, la.trasero, nu.tol)
        limitante = "masa_max_sensor"
        if not (math.isfinite(ell) and SM(float(mod.x_CG_con_trasero(ell, ell2))) >= SM_min - 1e-6):
            Ll.banderas.append("SM_inalcanzable_con_masa_max")
            ell = math.nan
        ell_u = ell
    # 4: tolerancia mínima del amarre
    if cfg.remolque.tol_min is not None and math.isfinite(ell_u):
        l1, l2 = _recortar_por_amarre(mod, ell_u, ell2, x_CP, D, SM_min, vu.q, S_ref, CNa,
                                      cfg.remolque.alpha_max, cfg.remolque.tol_min, nu.n_ell, nu.tol)
        if not math.isfinite(l1):
            Ll.banderas.append("tol_amarre_inalcanzable")
        elif (l1, l2) != (ell_u, ell2):
            limitante = "tol_amarre"
        ell_u, ell2 = l1, l2
    if not math.isfinite(ell_u):
        return Ll

    V1, V2 = float(mod.V_b(ell_u)), float(mod.V_tras(ell_u, ell2))
    m_u = float(mod.m_con_trasero(ell_u, ell2))
    xcg_u = float(mod.x_CG_con_trasero(ell_u, ell2))
    SM_u = SM(xcg_u)
    # amarre en el CG: α_trim = 0; NaN (no lineal) si el CP no está detrás del amarre
    alpha = 0.0 if x_CP > xcg_u else math.nan
    if not abs(alpha) <= cfg.remolque.alpha_lineal:
        Ll.banderas.append("trim_no_lineal")
    if SM_u > rest.SM_max + 1e-9:
        Ll.banderas.append("SM_sobre_max")
    Ll.ell, Ll.ell2, Ll.limitante = ell_u, ell2, limitante
    Ll.fila.update({
        "m_total_g": m_u / G, "m_lastre_g": la.rho_b * (V1 + V2) / G,
        "m_lastre_delantero_g": la.rho_b * V1 / G, "m_lastre_trasero_g": la.rho_b * V2 / G,
        "ell_mm": ell_u / MM, "ell_trasero_mm": ell2 / MM, "x_CG_mm": xcg_u / MM, "SM_cal": SM_u,
        "alpha_trim_deg": math.degrees(alpha), "x_T_mm": xcg_u / MM,
        "tol_amarre_mm": tolerancia_amarre(m_u, xcg_u, x_CP, vu.q, S_ref, CNa, cfg.remolque.alpha_max) / MM,
        "costo_lastre_usd": la.rho_b * (V1 + V2) * la.costo_usd_kg,
        "x_electronica_mm": float(mod.x_e(ell_u)) / MM,
    })
    return Ll
