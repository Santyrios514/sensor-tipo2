"""Aletas freeform sobre el tubo de cola: geometría, masa, Barrowman por franjas y flutter.

Parametrización (x global desde la punta; la raíz termina en el extremo del tubo, como el
`offset bottom = 0` del .ork):

    c_r = μ L_tc,   x_LE = L − c_r,   c_t = γ c_r,   x_s = σ (c_r − c_t)
    r_tip = (r_tip/R) R,   h = r_tip − r_tc

Puntos locales de OpenRocket (origen en el borde de ataque de la raíz, y desde la superficie
del tubo): P0 = (0, 0), P1 = (x_s, h), P2 = (x_s + c_t, h), P3 = (c_r, 0). Con σ ≤ 1 la punta no
pasa del extremo del tubo, así que el largo total es el del cuerpo.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .atmosfera import isa
from .config import AletaSpec, ConfigOpt, ParamsAleta
from .perfiles import Perfil


@dataclass(frozen=True)
class GeomAleta:
    params: ParamsAleta
    x_LE: float
    r_raiz: float  # radio del tubo de cola
    c_r: float
    c_t: float
    x_s: float
    h: float

    @property
    def r_tip(self) -> float:
        return self.r_raiz + self.h

    @property
    def x_TE(self) -> float:
        return self.x_LE + self.c_r

    @property
    def puntos(self) -> np.ndarray:
        return np.array([(0.0, 0.0), (self.x_s, self.h), (self.x_s + self.c_t, self.h), (self.c_r, 0.0)])

    @property
    def poligono(self) -> np.ndarray:
        """Polígono global (x, r); el cierre P3 → P0 va por la superficie del tubo."""
        return self.puntos + np.array([self.x_LE, self.r_raiz])

    @property
    def area(self) -> float:
        return area_centroide(self.poligono)[0]

    @property
    def x_cg(self) -> float:
        return area_centroide(self.poligono)[1]

    @property
    def masa(self) -> float:
        p = self.params
        return p.n * p.rho * p.phi * p.t * self.area


def area_centroide(P: np.ndarray) -> tuple[float, float]:
    """Área y x̄ de un polígono simple (fórmula del polígono)."""
    x, y = P[:, 0], P[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cruz = x * y1 - x1 * y
    A = 0.5 * cruz.sum()
    return abs(A), float(((x + x1) * cruz).sum() / (6.0 * A))


def construir(cfg: ConfigOpt, perfil: Perfil, a: AletaSpec) -> tuple[GeomAleta, list[str]]:
    """Geometría de la aleta del candidato y motivos de descarte (lista vacía si es válida)."""
    rest = cfg.restricciones
    c_r = a.mu_cr * perfil.L_tc
    c_t = a.gamma_ct * c_r
    g = GeomAleta(params=cfg.aleta, x_LE=perfil.L - c_r, r_raiz=perfil.r_tc, c_r=c_r, c_t=c_t,
                  x_s=a.sigma_flecha * (c_r - c_t), h=a.r_tip_rel_R * perfil.R - perfil.r_tc)
    motivos = []
    if not g.h >= rest.h_min - 1e-12:
        motivos.append("h_menor_minimo")
    if not c_r >= rest.c_min - 1e-12:
        motivos.append("c_r_menor_minimo")
    if not c_t >= rest.c_min - 1e-12:
        motivos.append("c_t_menor_minimo")
    if rest.D_ap_max is not None and 2 * max(perfil.R, g.r_tip) > rest.D_ap_max + 1e-12:
        motivos.append("D_ap_mayor_maximo")
    if rest.D_acostado_max_rel is not None and \
            max(envolvente(g.r_tip, perfil.R, cfg.aleta.n, cfg.rot_guardado)) > rest.D_acostado_max_rel * perfil.D + 1e-9:
        motivos.append("D_acostado_mayor_maximo")
    return g, motivos


def envolvente(r_tip: float, R: float, n: int, rotacion: float) -> tuple[float, float]:
    """(alto, ancho) de la sección: cuerpo de radio R y n aletas hasta r_tip (φ desde la vertical)."""
    phi = rotacion + 2 * np.pi * np.arange(n) / n
    c, s = r_tip * np.cos(phi), r_tip * np.sin(phi)
    return float(max(R, c.max()) + max(R, -c.min())), float(max(R, s.max()) + max(R, -s.min()))


# --------------------------------------------------------------------------- Barrowman (Niskanen §3.2.2)


@dataclass(frozen=True)
class CPAleta:
    CNa: float  # del conjunto, referido a A_ref
    x_CP: float  # global
    span: float
    mac: float
    cos_gamma: float


def barrowman(g: GeomAleta, A_ref: float, mach: float, n_franjas: int) -> CPAleta:
    """CP y C_Nα del conjunto de aletas por franjas en y, como FinSetCalc de OpenRocket 24.12:

        C_Nα,1 = 2π s²/A_ref / (1 + sqrt(1 + (β s² / (A_f cos Γ))²))
        C_Nα = C_Nα,1 · (n/2) · K_TB,   K_TB = 1 + r_t / (s + r_t)
        x_CP = x_LE + x_MAC,LE + MAC/4

    con s el span, A_f el área de una aleta, Γ la flecha de la línea de cuerdas medias y r_t el
    radio del tubo de cola. Las cuerdas se cortan contra la polilínea P0 → P3 (sin el cierre).
    """
    P = g.puntos
    span = g.h
    ys = np.linspace(0.0, span, n_franjas)
    lead = np.full(n_franjas, np.inf)
    trail = np.full(n_franjas, -np.inf)
    for (x1, y1), (x2, y2) in zip(P[:-1], P[1:]):
        if abs(y2 - y1) < 1e-12:
            continue
        lo, hi = min(y1, y2), max(y1, y2)
        m = (ys >= lo - 1e-12) & (ys <= hi + 1e-12)
        x = x1 + (ys[m] - y1) * (x2 - x1) / (y2 - y1)
        lead[m] = np.minimum(lead[m], x)
        trail[m] = np.maximum(trail[m], x)
    ok = np.isfinite(lead) & np.isfinite(trail) & (trail > lead)
    c = np.where(ok, trail - lead, 0.0)
    lead = np.where(ok, lead, 0.0)
    w = np.full(n_franjas, span / (n_franjas - 1))
    w[[0, -1]] /= 2
    A_fr = float((c * w).sum())
    mac = float((c**2 * w).sum() / A_fr)
    x_mac_le = float((lead * c * w).sum() / A_fr)
    tan_g = np.gradient((lead + c / 2)[ok], ys[ok])
    cos_g = float(np.mean(np.cos(np.arctan(tan_g))))
    beta = math.sqrt(max(1 - mach**2, 1e-6))
    CNa1 = 2 * math.pi * span**2 / A_ref / (1 + math.sqrt(1 + (beta * span**2 / (g.area * cos_g)) ** 2))
    n = g.params.n
    K_TB = 1 + g.r_raiz / (span + g.r_raiz)
    return CPAleta(CNa=CNa1 * (n / 2 if n >= 3 else 1.0) * K_TB, x_CP=g.x_LE + x_mac_le + 0.25 * mac,
                   span=span, mac=mac, cos_gamma=cos_g)


# --------------------------------------------------------------------------- flutter


def velocidad_flutter(g: GeomAleta, altitud: float) -> float:
    """Martin (NACA TN 4197) sobre la aleta trapezoidal [m/s]:

        V_f = a sqrt(G / (1.337 AR³ P (λ + 1) / (2 (AR + 2) (t/c_r)³))),
        AR = h²/A_f, λ = c_t/c_r, G = E / (2 (1 + ν))
    """
    p = g.params
    A_f = g.area
    if not (g.h > 0 and A_f > 0 and g.c_r > 0 and p.t > 0 and p.E > 0):
        return math.nan
    atm = isa(altitud)
    AR = g.h**2 / A_f
    den = 1.337 * AR**3 * atm.p * (g.c_t / g.c_r + 1) / (2 * (AR + 2) * (p.t / g.c_r) ** 3)
    return atm.a * math.sqrt(p.E / (2 * (1 + p.nu)) / den)
