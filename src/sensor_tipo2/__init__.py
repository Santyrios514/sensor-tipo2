"""Optimización del sensor remolcado tipo 2 (DBF 2026-27, UPB).

Nariz + cuerpo cilíndrico + transición + tubo de cola con aletas. Internamente todo está en SI
(m, kg, s, rad); la E/S (YAML, CSV) usa mm, g y grados.
"""

__version__ = "0.1.0"

G0 = 9.80665  # m/s^2, gravedad estándar
MM = 1e-3
G = 1e-3
