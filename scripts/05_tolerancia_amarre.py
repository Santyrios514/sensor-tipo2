#!/usr/bin/env python
"""Script 5: tolerancia de amarre de un .ork (OpenRocket 24.12, necesita Java).

    python scripts/05_tolerancia_amarre.py archivo.ork [--config config/optimizacion.yaml] [--mach M]

Con el amarre en el CG (trim estático 0), la tolerancia es el error con que hay que ubicar el
amarre respecto del CG para que el ángulo de trim no pase de α_max:

    tol = α_max q S_ref C_Nα (x_CP − x_CG) / (m g)

Masa y CG salen de OpenRocket (estructura + Mass components: plomo, electrónica, herraje); CP y C_Nα,
de Barrowman a α = 0. El herraje (Mass component "herraje…") se toma siempre en el CG, esté donde
esté en el .ork: cuenta su masa pero no mueve el CG. V, altitud (ISA), Mach, α_max (remolque.alpha_trim_max_deg) y la tolerancia
mínima (remolque.tol_amarre_min_mm) salen de la configuración.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from sensor_tipo2.config import ConfigError, cargar  # noqa: E402

MM, G = 1e-3, 1e-3


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("ork", help="archivo .ork")
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--mach", type=float, help="por defecto, el de condiciones_vuelo")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    if not Path(a.ork).exists():
        print(f"No existe {a.ork}", file=sys.stderr)
        return 2
    try:
        from sensor_tipo2.verificacion import tolerancia_amarre_ork
        r = tolerancia_amarre_ork(cfg, a.ork, a.mach)
    except Exception as e:  # sin Java/orlab o .ork ilegible
        print(f"No se pudo evaluar {a.ork} con OpenRocket: {e}", file=sys.stderr)
        return 2
    vu, tol_min = cfg.vuelo, cfg.remolque.tol_min
    print(f"{Path(a.ork).name}")
    print(f"  masa {r.m / G:.1f} g · largo {r.L / MM:.1f} mm · x_CG {r.x_CG / MM:.2f} mm · x_CP {r.x_CP / MM:.2f} mm · "
          f"C_Nα {r.CNa:.3f} · D_ref {r.D_ref / MM:.1f} mm · SM {r.SM:.3f} cal")
    print(f"  V {vu.V:g} m/s a {vu.altitud:g} m (ISA) · Mach {r.mach:.4f} · q {r.q:.1f} Pa · "
          f"α_max {math.degrees(r.alpha_max):g}°")
    if not math.isfinite(r.tol):
        print("  Tolerancia de amarre: no definida (el CP no está detrás del CG: inestable)")
        return 1
    tol_mm = round(r.tol / MM, 3)  # el veredicto se da sobre el valor que se imprime
    linea = f"  Tolerancia de amarre: {tol_mm:.3f} mm"
    if tol_min is not None:
        linea += f"  ({'cumple' if tol_mm >= tol_min / MM else 'NO cumple'} el mínimo de {tol_min / MM:g} mm)"
    print(linea)
    for xh in r.x_herraje:
        if abs(xh - r.x_CG) > 0.5 * MM:
            print(f"  herraje: en el .ork está en x = {xh / MM:.2f} mm; se tomó en el CG (x = {r.x_CG / MM:.2f} mm), "
                  "que es donde va el amarre. Conviene moverlo también en el .ork.")
    if not r.x_herraje:
        print("  aviso: el .ork no tiene un Mass component llamado 'herraje…': no se incluye su masa")
    if r.L > cfg.restricciones.L_max + 1e-9:
        print(f"  aviso: el largo total ({r.L / MM:.1f} mm) supera restricciones.L_max_mm = {cfg.restricciones.L_max / MM:g} mm")
    for w in r.advertencias:
        print(f"  aviso de OpenRocket: {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
