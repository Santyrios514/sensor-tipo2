"""CP del modelo propio (Barrowman interno) con guarda numérica. OpenRocket solo valida a los
ganadores (verificacion.py); no hay calibración.

    x_CP = Σ C_Nα,i x_i / Σ C_Nα,i,     Σ C_Nα = 2 + 2 (k² − 1) + C_Nα,f = 2 k² + C_Nα,f

La singularidad aparece cuando 2k² + C_Nα,f → 0; por eso el candidato se descarta si
Σ C_Nα < ε_CN, antes de dividir.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

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
