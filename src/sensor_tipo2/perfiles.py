"""Perfil exterior analítico r_e(x) del sensor tipo 2: nariz, cuerpo, transición y tubo de cola.

Las funciones de forma replican `Transition.Shape.getRadius(x, radius, length, param)` de
OpenRocket 24.12, con x medido desde el extremo de menor radio:

- nariz:                 r(x) = shape(x, R, L_n, p)
- transición:            r(u) = r_tc + shape(L_t − u, R − r_tc, L_t, p), u desde la unión con el cuerpo
- transición recortada:  r(u) = shape(c + L_t − u, R, c + L_t, p), con shape(c, R, c + L_t, p) = r_tc
  (`clipped`, el defecto de OpenRocket para elipsoide, potencia y Haack)

x global se mide desde la punta de la nariz, positivo hacia popa. Todo en metros.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy.optimize import brentq

FORMAS = ("conica", "elipsoide", "ogiva", "potencia", "parabolica", "haack")
RECORTABLES = ("elipsoide", "potencia", "haack")  # Transition.Shape.isClippable()
PARAMETRO_DEFECTO = {"ogiva": 1.0, "parabolica": 1.0, "haack": 0.0}
FORMA_OR = {"conica": "CONICAL", "ogiva": "OGIVE", "elipsoide": "ELLIPSOID", "potencia": "POWER",
            "parabolica": "PARABOLIC", "haack": "HAACK"}
ESTACIONES = ("nariz", "cuerpo", "cola", "tubo_cola")


def forma_or(x, forma: str, radius: float, length: float, param: float | None) -> np.ndarray:
    """Radio de la forma `forma` a la distancia x del extremo fino (convención de OpenRocket)."""
    x = np.clip(np.asarray(x, dtype=float), 0.0, length)
    xi = x / length
    if forma == "conica":
        return radius * xi
    if forma == "ogiva":
        p = 1.0 if param is None else float(param)
        if p < 0.001:
            return radius * xi
        Rc = np.sqrt((length**2 + radius**2) * (((2 - p) * length) ** 2 + (p * radius) ** 2)
                     / (4 * (p * radius) ** 2))
        Lc = length / p
        y0 = np.sqrt(max(Rc**2 - Lc**2, 0.0))
        return np.sqrt(np.clip(Rc**2 - (Lc - x) ** 2, 0.0, None)) - y0
    if forma == "elipsoide":
        return radius * np.sqrt(np.clip(2 * xi - xi**2, 0.0, None))
    if forma == "potencia":
        return radius * xi ** float(param)
    if forma == "parabolica":
        K = 1.0 if param is None else float(param)
        return radius * (2 * xi - K * xi**2) / (2 - K)
    if forma == "haack":
        C = 0.0 if param is None else float(param)
        th = np.arccos(1.0 - 2.0 * xi)
        return radius * np.sqrt(np.clip((th - np.sin(2 * th) / 2 + C * np.sin(th) ** 3) / np.pi, 0.0, None))
    raise ValueError(f"forma desconocida: {forma}")


def longitud_recorte(forma: str, R: float, Ra: float, Lt: float, parametro: float | None) -> float:
    """c ≥ 0 tal que la forma de radio R y longitud c + L_t vale R_a en c (Transition.calculateClip)."""
    if Ra <= 0:
        return 0.0

    def f(c):
        return float(forma_or(np.array([c]), forma, R, c + Lt, parametro)[0]) - Ra

    hi = Lt
    while f(hi) < 0:
        hi *= 2
        if hi > 1e6 * Lt:
            raise ValueError(f"no se pudo recortar la forma '{forma}' (R = {R}, R_a = {Ra})")
    return brentq(f, 0.0, hi, xtol=1e-15, rtol=1e-14)


@dataclass(frozen=True)
class Perfil:
    """Geometría exterior resuelta (SI). L es el largo total: nariz + cuerpo + transición + tubo."""

    L: float
    D: float
    Ln: float
    forma_n: str
    param_n: float | None
    Lt: float
    forma_t: str
    param_t: float | None
    recortada: bool
    d_tc: float
    L_tc: float

    @property
    def R(self) -> float:
        return self.D / 2

    @property
    def r_tc(self) -> float:
        return self.d_tc / 2

    @property
    def k(self) -> float:
        return self.d_tc / self.D

    @property
    def x_t0(self) -> float:
        """Inicio de la transición."""
        return self.L - self.L_tc - self.Lt

    @property
    def x_tc0(self) -> float:
        """Inicio del tubo de cola."""
        return self.L - self.L_tc

    @property
    def L_c(self) -> float:
        """Largo del cuerpo cilíndrico."""
        return self.x_t0 - self.Ln

    @property
    def uniones(self) -> tuple[float, float, float]:
        return (self.Ln, self.x_t0, self.x_tc0)

    @cached_property
    def c_recorte(self) -> float:
        if self.recortada and self.forma_t in RECORTABLES and self.r_tc < self.R:
            return longitud_recorte(self.forma_t, self.R, self.r_tc, self.Lt, self.param_t)
        return 0.0

    def radio_transicion(self, u) -> np.ndarray:
        """Radio de la transición con u = x − x_t0 ∈ [0, L_t] (se estrecha hacia popa)."""
        u = np.clip(np.asarray(u, dtype=float), 0.0, self.Lt)
        if self.Lt == 0 or self.R == self.r_tc:
            return np.full_like(u, self.R)
        if self.recortada and self.forma_t in RECORTABLES:
            c = self.c_recorte
            return forma_or(c + self.Lt - u, self.forma_t, self.R, c + self.Lt, self.param_t)
        return self.r_tc + forma_or(self.Lt - u, self.forma_t, self.R - self.r_tc, self.Lt, self.param_t)

    def radio(self, x) -> np.ndarray:
        """r_e(x) en [0, L]; 0 fuera de ese intervalo."""
        x = np.asarray(x, dtype=float)
        r = np.full_like(x, self.R)
        m_n = x < self.Ln
        r[m_n] = forma_or(x[m_n], self.forma_n, self.R, self.Ln, self.param_n)
        m_t = (x > self.x_t0) & (x < self.x_tc0)
        r[m_t] = self.radio_transicion(x[m_t] - self.x_t0)
        r[x >= self.x_tc0] = self.r_tc
        r[(x < 0) | (x > self.L)] = 0.0
        return r

    def estacion(self, x) -> np.ndarray:
        """0 nariz, 1 cuerpo, 2 transición, 3 tubo de cola."""
        x = np.asarray(x, dtype=float)
        return np.searchsorted(np.array(self.uniones), x, side="right")
