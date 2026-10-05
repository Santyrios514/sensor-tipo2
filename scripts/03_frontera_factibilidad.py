#!/usr/bin/env python
"""Script 3 (spec v2 §6): estudio de factibilidad. Se corre siempre.

    python scripts/03_frontera_factibilidad.py [--config config/optimizacion.yaml] [--completa]
        [--procesos N]

Para cada tope de r_tip/R ∈ {1.0, 1.1, 1.2, 1.3, 1.4, 1.6} y n ∈ {4, 6, 8} aletas: número de
candidatos factibles y el de mayor masa (D_ap, SM, tolerancia de amarre, cand_id). Las filas con
tope > restricciones.r_tip_rel_R_max o n ≠ geometria_fija.aletas.n son solo diagnóstico.
Por defecto usa toda la malla de la configuración (reutiliza la caché de cuerpos para las tres n);
--gruesa usa un valor de cada dos en D, tubo y L_tc, para una corrida rápida.

Salidas: frontera_factibilidad.csv, casi_factibles.csv (si no hay factibles con las restricciones
reales) y masa_vs_tope.png.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from sensor_tipo2 import frontera  # noqa: E402
from sensor_tipo2.config import ConfigError, cargar  # noqa: E402
from sensor_tipo2.exportar import escribir_csv  # noqa: E402


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--gruesa", action="store_true", help="un valor de cada dos en D, tubo y L_tc")
    p.add_argument("--procesos", type=int)
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    t0 = time.time()

    def progreso(i, total):
        if i == total or i % max(1, total // 10) == 0:
            print(f"  {i}/{total} cuerpos · {time.time() - t0:.0f} s", flush=True)

    m = frontera.malla_frontera(cfg, gruesa=a.gruesa)
    print(f"Malla {'gruesa' if a.gruesa else 'completa'}: D {m['D_mm']}, {cfg.var_tubo} {m[cfg.var_tubo]}, "
          f"L_tc {m['L_tc_mm']}, r_tip/R {m['r_tip_rel_R']}; n = {list(frontera.NS)}")
    tab, cf, n_eval = frontera.evaluar(cfg, gruesa=a.gruesa, procesos=a.procesos, progreso=progreso)
    dir_d, dir_f = cfg.dir_salida(), cfg.dir_figuras()
    escribir_csv(tab, dir_d / "frontera_factibilidad.csv")
    frontera.fig_masa_vs_tope(tab, cfg, dir_f / "masa_vs_tope.png")
    cols = ["n_aletas", "tope_r_tip_rel_R", "n_factibles", "m_total_max_g", "D_ap_mm", "SM_cal", "tol_amarre_mm",
            "solo_diagnostico"]
    print(f"\n{n_eval:,} evaluaciones · {time.time() - t0:.0f} s\n")
    print(tab.reindex(columns=cols).to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    real = tab[~tab["solo_diagnostico"]]
    n_real = int(real["n_factibles"].max()) if len(real) else 0
    if n_real == 0:
        escribir_csv(cf, dir_d / "casi_factibles.csv")
        print(f"\nSin factibles con las restricciones reales → casi_factibles.csv ({len(cf)} candidatos)")
    else:
        print(f"\nHay {n_real:,} factibles con las restricciones reales (n = {cfg.aleta.n}, "
              f"r_tip/R ≤ {cfg.restricciones.r_tip_rel_R_max:g}) en esta malla.")
    print(f"Datos → {dir_d} · Figuras → {dir_f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
