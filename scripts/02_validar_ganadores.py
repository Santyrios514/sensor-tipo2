#!/usr/bin/env python
"""Script 2 (spec v2 rev. 2 §7, opcional): validación de los ganadores con OpenRocket 24.12.

    python scripts/02_validar_ganadores.py [--config config/optimizacion.yaml] [--ranking data_opt/ranking.csv]
        [--cand ID ...]

Arma en OpenRocket (orlab + JPype, plantilla analisis_tipo2.ork) los `ejecucion.validacion_or.N_ganadores`
mejores factibles del ranking, con el plomo, la electrónica y el herraje como Mass components, y
compara CP, C_Nα, masa y CG con el modelo propio. Escribe validacion_or.csv y marca en ranking.csv
`validado_or` y `ganador` (el primero validado). Cada --cand se valida además (rol "extra"), sin
entrar en la elección del ganador. No guarda ningún `.ork`.

Códigos de salida: 0 (con o sin ganador), 2 (configuración, ranking o OpenRocket no disponibles),
5 (algún candidato supera tol_cp_mm: el modelo propio dejó de coincidir con OpenRocket; se detiene
y se reporta, sin calibrar).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from sensor_tipo2 import exportar  # noqa: E402
from sensor_tipo2.barrido import RANKING  # noqa: E402
from sensor_tipo2.config import ConfigError, cargar  # noqa: E402


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--ranking", help="por defecto, {salida.dir}/ranking.csv")
    p.add_argument("--cand", action="append", default=[], help="valida además este cand_id (se puede repetir)")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    dir_d = cfg.dir_salida()
    ruta_rank = Path(a.ranking) if a.ranking else dir_d / "ranking.csv"
    if not ruta_rank.exists():
        print(f"No existe {ruta_rank}: corra primero scripts/01_optimizar_malla.py", file=sys.stderr)
        return 2
    rk = exportar.leer_ranking(ruta_rank)
    try:
        from sensor_tipo2.verificacion import VALIDACION, abrir_openrocket, validar_candidato, validar_ganadores
        ctx = abrir_openrocket(cfg)
        pr = ctx.__enter__()
    except Exception as e:  # sin orlab, JPype, Java o jar
        print(f"OpenRocket no disponible ({e}). La validación es opcional: los planos quedan 'sin validar'.",
              file=sys.stderr)
        return 2
    try:
        vo = cfg.ejecucion.get("validacion_or") or {}
        print(f"OpenRocket {pr.version_or} · validando los {vo.get('N_ganadores', 5)} mejores factibles de {ruta_rank.name} "
              f"(|Δx_CP| ≤ {vo.get('tol_cp_mm', 2.0)} mm, |Δm| ≤ {vo.get('tol_masa_pct', 0.5)} %, SM_OR ∈ "
              f"[{cfg.restricciones.SM_min:g}, {cfg.restricciones.SM_max:g}])")
        res = validar_ganadores(cfg, pr, rk)
        res.tabla["rol"] = "ganadores"
        extra = []
        filas = rk.set_index("cand_id")
        for cid in a.cand:
            if cid in set(res.tabla["cand_id"]):
                continue
            if cid not in filas.index:
                print(f"--cand {cid}: no está en {ruta_rank}", file=sys.stderr)
                continue
            f = filas.loc[cid].copy()
            f["cand_id"] = cid
            r = validar_candidato(cfg, pr, f, float(vo.get("tol_cp_mm", 2.0)) * 1e-3, float(vo.get("tol_masa_pct", 0.5)))
            r["rol"] = "extra"
            print(f"    extra {cid}  Δx_CP {r.get('dx_cp_mm', float('nan')):+.3f} mm  SM_OR {r.get('SM_or_cal', float('nan')):.3f}  "
                  + ("validado" if r["validado"] else f"NO ({r.get('motivo', '')})"))
            extra.append(r)
        if extra:
            import pandas as pd
            res.tabla = pd.concat([res.tabla, pd.DataFrame(extra)], ignore_index=True)
    finally:
        ctx.__exit__(None, None, None)
    exportar.escribir_csv(res.tabla, dir_d / "validacion_or.csv", VALIDACION + ["rol"])
    val = dict(zip(res.tabla["cand_id"], res.tabla["validado"].astype(bool))) if len(res.tabla) else {}
    rk["validado_or"] = rk["cand_id"].map(val)
    rk["ganador"] = rk["cand_id"] == res.ganador
    exportar.escribir_csv(rk, ruta_rank, RANKING + ["validado_or", "ganador"])
    print(f"Datos → {dir_d / 'validacion_or.csv'}")
    if res.cp_fuera_de_tolerancia:
        print("\nDETENIDO: el CP del modelo propio difiere de OpenRocket más de tol_cp_mm en "
              f"{', '.join(res.cp_fuera_de_tolerancia)}. El modelo dejó de coincidir para esa familia de formas; "
              "revíselo antes de usar el ranking (no se calibra en silencio). Sin ganador.", file=sys.stderr)
        return 5
    if res.ganador is None:
        print("\nSin ganador: ninguno de los candidatos validados cumple el criterio (ver validacion_or.csv).")
        return 0
    v = res.tabla.set_index("cand_id").loc[res.ganador]
    print(f"\nGanador (validado con OpenRocket): {res.ganador}\n"
          f"  masa {v['m_modelo_g']:.0f} g (OR {v['m_or_g']:.0f} g, {v['dif_masa_or_pct']:+.3f} %) · "
          f"x_CP {v['x_CP_modelo_mm']:.2f} / {v['x_CP_or_mm']:.2f} mm · SM {v['SM_modelo_cal']:.3f} / {v['SM_or_cal']:.3f} · "
          f"C_D (OR) {v['CD_or']:.3f}")
    print("Planos con la línea de validación: python scripts/04_planos.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
