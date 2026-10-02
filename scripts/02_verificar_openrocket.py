#!/usr/bin/env python
"""Script 2: verificación y calibración con OpenRocket 24.12 del ranking del script 1.

    python scripts/02_verificar_openrocket.py [--config config/optimizacion.yaml]
        [--ranking data_opt/ranking.csv] [--procesos N] [--sin-figuras]

Verifica los mejores y una muestra estratificada, calibra el CP sustituto, recalcula la malla si
hace falta, elige el ganador con el CP de OpenRocket y reescribe ranking.csv. Escribe además
verificacion_or.csv, calibracion.json y los `.ork` de los mejores (lastre, electrónica y herraje
como Mass components), y comprueba que el `.ork` del ganador reproduce su CP al reabrirlo.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from sensor_tipo2 import exportar  # noqa: E402
from sensor_tipo2.barrido import RANKING  # noqa: E402
from sensor_tipo2.config import ConfigError, cargar  # noqa: E402
from sensor_tipo2.verificacion import VERIFICACION, PuenteTipo2, verificar_candidato, verificar_y_calibrar  # noqa: E402

MM = 1e-3
TOL_REAPERTURA_MM = 0.5


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--ranking", help="por defecto, {salida.dir}/ranking.csv")
    p.add_argument("--procesos", type=int)
    p.add_argument("--sin-figuras", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    dir_d, dir_f = cfg.dir_salida(), cfg.dir_figuras()
    ruta_rank = Path(a.ranking) if a.ranking else dir_d / "ranking.csv"
    if not ruta_rank.exists():
        print(f"No existe {ruta_rank}: corra primero scripts/01_optimizar_malla.py", file=sys.stderr)
        return 2
    ranking = exportar.leer_ranking(ruta_rank)
    orc = cfg.openrocket
    ork = cfg.resolver(orc.get("ork_base", "modelos/analisis_tipo2.ork"))

    import orlab
    with orlab.OpenRocketInstance(jar_path=orc.get("jar"), log_level="ERROR") as inst:
        pr = PuenteTipo2(inst, ork)
        print(f"OpenRocket {pr.version_or} · {len(ranking):,} candidatos en {ruta_rank}")
        res = verificar_y_calibrar(cfg, pr, ranking, a.procesos)
        df, ver, cal = res.ranking, res.verificacion, res.calibracion
        df["ganador"] = df["cand_id"] == res.ganador
        n_ork = int(cfg.salida.get("exportar_ork_top", 5))
        mejores = (ver[ver["factible_or"].fillna(False).astype(bool)].sort_values(["J_or", "cand_id"]).head(n_ork)
                   if "factible_or" in ver else ver.iloc[0:0])
        filas = df.set_index("cand_id")
        reapertura = math.nan
        for i, cid in enumerate(mejores["cand_id"], start=1):
            ruta = dir_d / "ork" / f"rank{i:02d}_{cid}.ork"
            fila = filas.loc[cid].copy()
            fila["cand_id"] = cid
            r = verificar_candidato(cfg, pr, fila, cal, guardar=ruta)
            if cid == res.ganador:
                pr.cargar(ruta)
                reapertura = abs(pr.aero(cfg.vuelo.mach).x_CP / MM - r["x_CP_or_mm"])
            print(f"  .ork → {ruta.name}")

    exportar.escribir_csv(df, dir_d / "ranking.csv", RANKING + ["ganador"])
    exportar.escribir_csv(df[df["en_pareto"].astype(bool)], dir_d / "pareto.csv", RANKING)
    exportar.escribir_csv(ver, dir_d / "verificacion_or.csv", VERIFICACION)
    exportar.escribir_json({**cal.a_dict(), "historial": res.historial, "ganador": res.ganador,
                            "reapertura_ganador_dif_CP_mm": reapertura}, dir_d / "calibracion.json")
    if res.ganador is None:
        print("\nNingún candidato verificado es factible con el CP de OpenRocket.")
        return 1
    g = df[df["ganador"]].iloc[0]
    print(f"\nGanador (verificado con OpenRocket): {g['cand_id']}\n"
          f"  m_total {g['m_total_or_g']:.0f} g · D {g['D_mm']:.1f} mm · d_tc {g['d_tc_mm']:.1f} mm · "
          f"D acostado {g['D_acostado_mm']:.1f} mm · D aparente {g['D_ap_mm']:.1f} mm · SM_OR {g['SM_or_cal']:.2f} (sustituto {g['SM_cal']:.2f}) · "
          f"x_CP OR {g['x_CP_or_mm']:.1f} mm · C_D OR {g['CD_or']:.3f} · dif. masa {g['dif_masa_or_pct']:+.3f} %")
    print(f"  .ork del ganador reabierto: |ΔCP| = {reapertura:.3f} mm (tolerancia {TOL_REAPERTURA_MM} mm)")
    if not a.sin_figuras:
        dir_f.mkdir(parents=True, exist_ok=True)
        exportar.fig_calibracion(ver, cal, dir_f / "calibracion.png")
        exportar.fig_pareto(df.sort_values(["ganador"], ascending=False), dir_f / "pareto.png", cfg.objetivo.col_diametro)
        exportar.fig_factibilidad(df, dir_f / "factibilidad.png")
        v = ver.set_index("cand_id").loc[res.ganador]
        exportar.fig_ganador(cfg, g, dir_f, cal, x_CP=g["x_CP_or_mm"] * MM, CNa=float(v["CNa_or"]))
        exportar.fig_sensibilidad(cfg, g, dir_f / "sensibilidad.png", cal)
        print(f"Figuras → {dir_f}")
    print(f"Datos → {dir_d}")
    return 0 if reapertura <= TOL_REAPERTURA_MM else 4


if __name__ == "__main__":
    sys.exit(main())
