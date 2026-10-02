#!/usr/bin/env python
"""Script 1: barrido en malla de la geometría del sensor tipo 2 con el CP sustituto (sin JVM).

    python scripts/01_optimizar_malla.py [--config config/optimizacion.yaml] [--forzar]
        [--sin-refinamiento] [--procesos N] [--sin-figuras]

Salidas en `salida.dir` (ranking.csv, pareto.csv, tolerancias_implicitas.json) y figuras en
`salida.dir_figuras`. El script 2 verifica con OpenRocket y reescribe el ranking.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from sensor_tipo2 import barrido, exportar  # noqa: E402
from sensor_tipo2.config import ConfigError, cargar  # noqa: E402
from sensor_tipo2.objetivo import tolerancias  # noqa: E402


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--forzar", action="store_true", help="corre aunque se supere max_evaluaciones")
    p.add_argument("--sin-refinamiento", action="store_true")
    p.add_argument("--procesos", type=int, help="sobrescribe ejecucion.procesos")
    p.add_argument("--sin-figuras", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    n_c, n_a = len(cfg.cuerpos()), len(cfg.aletas())
    n = n_c * n_a
    procesos = a.procesos or cfg.procesos
    lim = int(cfg.ejecucion.get("max_evaluaciones", 2_000_000))
    print(f"Malla: {n_c:,} cuerpos × {n_a} aletas = {n:,} evaluaciones (límite {lim:,}); {procesos} procesos.")
    print(f"Estimado: ~{n * 1.2e-3 / procesos / 60:.1f} min a ~1.2 ms por candidato (sin refinamiento).")
    if n > lim and not a.forzar:
        largos = sorted(((len(v), k) for k, v in cfg.malla.items() if isinstance(v, list)), reverse=True)[:3]
        print("La malla supera max_evaluaciones. Reduzca las listas más largas ("
              + ", ".join(f"{k} ({m} valores)" for m, k in largos) + ") o use --forzar.", file=sys.stderr)
        return 3

    t0 = time.time()

    def progreso(i, total):
        if i == total or i % max(1, total // 20) == 0:
            print(f"  {i}/{total} cuerpos · {time.time() - t0:.0f} s", flush=True)

    df = barrido.optimizar(cfg, refinar=not a.sin_refinamiento, procesos=procesos, progreso=progreso)
    dir_d, dir_f = cfg.dir_salida(), cfg.dir_figuras()
    exportar.exportar_ranking(cfg, df, dir_d)
    fac = df[df["factible"]]
    print(f"\n{len(df):,} candidatos ({int(df['refinamiento'].sum()):,} del refinamiento), {len(fac):,} factibles, "
          f"{int(df['en_pareto'].sum()):,} en el frente de Pareto · {time.time() - t0:.0f} s")
    print("\nTolerancias implícitas de los pesos:")
    for t in tolerancias(cfg)["lectura"]:
        print("  - " + t)
    if fac.empty:
        print("\nNingún candidato factible.")
    else:
        g = fac.iloc[0]
        print(f"\nMejor (CP sustituto, sin verificar): {g['cand_id']}\n"
              f"  m_total {g['m_total_g']:.0f} g (plomo {g['m_lastre_g']:.0f} g) · D {g['D_mm']:.1f} mm · "
              f"d_tc {g['d_tc_mm']:.1f} mm · D_ap {g['D_ap_mm']:.1f} mm (r_tip/R = {g['r_tip_rel_R']:g}) · "
              f"SM {g['SM_cal']:.2f} · tol. amarre {g['tol_amarre_mm']:.2f} mm · restricción activa: "
              f"{g['restriccion_activa']}")
        r1 = df[(df["r_tip_rel_R"] <= 1.0 + 1e-9)]
        print(f"  Con aletas dentro del diámetro (r_tip ≤ R): {int(r1['factible'].sum())} factibles de {len(r1):,}.")
    if not a.sin_figuras:
        rutas = exportar.figuras_barrido(cfg, df, dir_f)
        print(f"\n{len(rutas)} figuras → {dir_f}")
    print(f"Datos → {dir_d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
