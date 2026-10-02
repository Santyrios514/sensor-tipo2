"""CP sustituto (Barrowman interno) con guarda numérica y calibración lineal contra OpenRocket.

    x_CP = Σ C_Nα,i x_i / Σ C_Nα,i,     Σ C_Nα = 2 + 2 (k² − 1) + C_Nα,f = 2 k² + C_Nα,f

La singularidad aparece cuando 2k² + C_Nα,f → 0; por eso el candidato se descarta si
Σ C_Nα < ε_CN, antes de dividir.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .aletas import GeomAleta, barrowman
from .geometria import Contribucion, Cuerpo

MOTIVO_CN = "CN_total_pequeno"


@dataclass(frozen=True)
class ResultadoCP:
    x_CP: float
    CNa: float
    partes: tuple[Contribucion, ...]

    def con_x(self, x: float) -> "ResultadoCP":
        return ResultadoCP(x_CP=x, CNa=self.CNa, partes=self.partes)


def combinar(partes: tuple[Contribucion, ...], eps_CN: float) -> tuple[ResultadoCP | None, str | None]:
    """Suma de contribuciones con la guarda Σ C_Nα ≥ ε_CN. Nunca devuelve inf ni NaN."""
    CNa = sum(p.CNa for p in partes)
    if not (math.isfinite(CNa) and CNa >= eps_CN):
        return None, MOTIVO_CN
    x = sum(p.CNa * p.x for p in partes) / CNa
    if not math.isfinite(x):
        return None, MOTIVO_CN
    return ResultadoCP(x_CP=x, CNa=CNa, partes=partes), None


def cp_sustituto(cu: Cuerpo, g: GeomAleta, mach: float, n_franjas: int,
                 eps_CN: float) -> tuple[ResultadoCP | None, str | None]:
    """Nariz y transición del caché del cuerpo + aletas por franjas."""
    f = barrowman(g, math.pi * cu.spec.D**2 / 4, mach, n_franjas)
    return combinar(cu.partes + (Contribucion("aletas", f.CNa, f.x_CP),), eps_CN)


@dataclass(frozen=True)
class Calibracion:
    """x_CP,cal = α + β x_CP,sust (m)."""

    alpha: float = 0.0
    beta: float = 1.0
    R2: float = math.nan
    res_max: float = math.nan  # m
    n: int = 0
    iteracion: int = 0

    def aplicar(self, x):
        return self.alpha + self.beta * np.asarray(x, dtype=float)

    @property
    def identidad(self) -> bool:
        return self.alpha == 0.0 and self.beta == 1.0

    def a_dict(self) -> dict:
        return {"alpha_mm": self.alpha * 1e3, "beta": self.beta, "R2": self.R2,
                "residuo_max_mm": self.res_max * 1e3, "n": self.n, "iteracion": self.iteracion}


IDENTIDAD = Calibracion()


def ajustar(x_sust, x_or, iteracion: int = 0) -> Calibracion:
    """Mínimos cuadrados de x_OR ≈ α + β x_sust."""
    x = np.asarray(x_sust, dtype=float)
    y = np.asarray(x_or, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 2 or np.ptp(x) == 0:
        raise ValueError("se necesitan al menos dos x_CP sustitutos distintos para calibrar")
    A = np.column_stack([np.ones_like(x), x])
    (alpha, beta), *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - (alpha + beta * x)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    R2 = 1.0 - float((r**2).sum()) / ss_tot if ss_tot > 0 else 1.0
    return Calibracion(alpha=float(alpha), beta=float(beta), R2=R2, res_max=float(np.abs(r).max()),
                       n=int(x.size), iteracion=iteracion)
