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
    def AR(self) -> float:
        """Alargamiento de una aleta, h² / A_f."""
        return self.h**2 / self.area

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
    if not c_r >= max(rest.c_min, rest.c_r_min) - 1e-12:
        motivos.append("c_r_menor_minimo")
    if not c_t >= rest.c_min - 1e-12:
        motivos.append("c_t_menor_minimo")
    if rest.r_tip_rel_R_max is not None and a.r_tip_rel_R > rest.r_tip_rel_R_max + 1e-12:
        motivos.append("r_tip_sobre_tope")
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


def extremos_seccion(R: float, r_tc: float, r_tip: float, t: float, n: int, phi0) -> np.ndarray:
    """(x_min, x_max, y_min, y_max) exactos de la sección transversal: el cuerpo (círculo de radio R)
    y n aletas rectangulares de espesor t entre r_tc y r_tip, con la aleta 0 a phi0 de la vertical.
    Con phi0 vector (m,), devuelve un arreglo (m, 4)."""
    phi0 = np.asarray(phi0, dtype=float)
    phi = phi0[..., None] + 2 * np.pi * np.arange(n) / n
    s_, c_ = np.sin(phi), np.cos(phi)
    xs, ys = [], []
    for r in (r_tc, r_tip):
        for sg in (1.0, -1.0):
            # r u + sg t/2 w, con u = (sin φ, cos φ) radial y w = (cos φ, −sin φ) normal a la aleta
            xs.append(r * s_ + sg * t / 2 * c_)
            ys.append(r * c_ - sg * t / 2 * s_)
    x, y = np.concatenate(xs, axis=-1), np.concatenate(ys, axis=-1)
    return np.stack([np.minimum(-R, x.min(-1)), np.maximum(R, x.max(-1)),
                     np.minimum(-R, y.min(-1)), np.maximum(R, y.max(-1))], axis=-1)


@dataclass(frozen=True)
class CajaMinima:
    """Caja de menor altura que contiene al sensor acostado (eje horizontal), girándolo sobre su eje."""

    alto: float  # altura aparente H_ap [m]
    ancho: float  # ancho de la caja con ese giro [m]
    phi: float  # giro de la aleta 0 respecto de la vertical que da H_ap [rad]
    extremos: tuple[float, float, float, float]  # (x_min, x_max, y_min, y_max) con ese giro [m]


def altura_aparente(R: float, r_tc: float, r_tip: float, t: float, n: int, phi_vuelo: float,
                    n_muestras: int) -> CajaMinima:
    """H_ap = min_φ [y_max(φ) − y_min(φ)]: la altura mínima de una caja que contiene al sensor acostado.

    La sección tiene periodo 2π/n en φ. Se muestrea un periodo con n_muestras puntos y se refina el
    mínimo con una búsqueda acotada. Si el giro de vuelo (`phi_vuelo`) ya da la altura mínima, se usa
    ese giro (con 4 aletas a 45° y r_tip ≤ √2 R − t/2 la altura es la del cuerpo, 2R)."""
    from scipy.optimize import minimize_scalar

    def alto(phi):
        e = extremos_seccion(R, r_tc, r_tip, t, n, phi)
        return float(e[3] - e[2])

    e = extremos_seccion(R, r_tc, r_tip, t, n, phi_vuelo)
    if e[3] - e[2] <= 2 * R + 1e-15:  # H_ap ≥ 2R siempre: el giro de vuelo ya es el mínimo
        return CajaMinima(alto=float(e[3] - e[2]), ancho=float(e[1] - e[0]), phi=float(phi_vuelo),
                          extremos=tuple(float(v) for v in e))
    T = 2 * np.pi / n
    phis = phi_vuelo + np.linspace(0.0, T, n_muestras, endpoint=False)
    E = extremos_seccion(R, r_tc, r_tip, t, n, phis)
    H = E[:, 3] - E[:, 2]
    i = int(np.argmin(H))
    paso = T / n_muestras
    r = minimize_scalar(alto, bounds=(phis[i] - paso, phis[i] + paso), method="bounded",
                        options={"xatol": 1e-10})
    phi, H_min = (float(r.x), float(r.fun)) if r.fun < H[i] else (float(phis[i]), float(H[i]))
    if alto(phi_vuelo) <= H_min + 1e-12:
        phi, H_min = phi_vuelo, alto(phi_vuelo)
    e = extremos_seccion(R, r_tc, r_tip, t, n, phi)
    return CajaMinima(alto=float(H_min), ancho=float(e[1] - e[0]), phi=phi, extremos=tuple(float(v) for v in e))


# --------------------------------------------------------------------------- Barrowman (Niskanen §3.2.2)

# Interferencia aleta-aleta de OpenRocket 24.12 (FinSetCalc.calculateNonaxialForces, leído del código
# fuente de release-24.12): factor sobre C_Nα según cuántas aletas comparten la misma estación.
INTERFERENCIA_ALETAS = {5: 0.948, 6: 0.913, 7: 0.854, 8: 0.81}
INTERFERENCIA_MAS_DE_8 = 0.75


def factor_interferencia(n: int) -> float:
    if n <= 4:
        return 1.0
    return INTERFERENCIA_ALETAS.get(n, INTERFERENCIA_MAS_DE_8)


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
        C_Nα = C_Nα,1 · (n/2) · f_n · K_TB,   K_TB = 1 + r_t / (s + r_t)
        x_CP = x_LE + x_MAC,LE + MAC/4

    con s el span, A_f el área de una aleta, Γ la flecha de la línea de cuerdas medias y r_t el
    radio del tubo de cola y f_n la interferencia aleta-aleta de OpenRocket (1 hasta 4 aletas, 0.913
    con 6, 0.81 con 8). OpenRocket suma sin²(θ − φ_i) sobre las aletas, que vale n/2 para n ≥ 3
    aletas equiespaciadas. Las cuerdas se cortan contra la polilínea P0 → P3 (sin el cierre).
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
    return CPAleta(CNa=CNa1 * (n / 2 if n >= 3 else 1.0) * factor_interferencia(n) * K_TB,
                   x_CP=g.x_LE + x_mac_le + 0.25 * mac,
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
