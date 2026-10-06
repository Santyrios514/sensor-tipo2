#!/usr/bin/env python
"""Script 4 (spec v2 rev. 2 §8): planos de los mejores candidatos, sin JVM.

    python scripts/04_planos.py [--config config/optimizacion.yaml] [--ranking data_opt/ranking.csv]
        [--cand ID ...] [--n N]

Para los `salida.N_planos` mejores factibles del ranking (y para cada --cand del ranking) dibuja
`{salida.dir_planos}/rankNN_<cand_id>.png` (300 dpi) y `.pdf`, escribe `{salida.dir}/rankNN_perfil.csv`
(x, r_e, r_i y el polígono de la aleta, para CAD/CFD) y el comparativo de siluetas
`comparativo_top.png`. Si el script 2 ya corrió, el cajetín muestra la validación con
OpenRocket; si no, "sin validar". Sin factibles, dibuja los 3 primeros de casi_factibles.csv con el
rótulo INFACTIBLE. No escribe ningún `.ork`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import pandas as pd  # noqa: E402

from sensor_tipo2 import exportar, frontera, planos  # noqa: E402
from sensor_tipo2.barrido import evaluar_cuerpo, specs_de_fila  # noqa: E402
from sensor_tipo2.config import ConfigError, cargar  # noqa: E402

N_CASI_PLANOS = 3


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--ranking", help="por defecto, {salida.dir}/ranking.csv")
    p.add_argument("--cand", action="append", default=[], help="cand_id del ranking (se puede repetir)")
    p.add_argument("--n", type=int, help="sobrescribe salida.N_planos")
    return p.parse_args(argv)


def _validacion(dir_d: Path) -> dict[str, dict]:
    ruta = dir_d / "validacion_or.csv"
    if not ruta.exists():
        return {}
    v = pd.read_csv(ruta)
    if "validado" in v:
        v["validado"] = v["validado"].astype(str).str.lower().isin(["true", "1"])
    return {r["cand_id"]: r.to_dict() for _, r in v.iterrows() if pd.notna(r.get("x_CP_or_mm"))}


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    dir_d = cfg.dir_salida()
    dir_p = cfg.resolver(cfg.salida.get("dir_planos", "planos"))
    dpi = int(cfg.salida.get("dpi_planos", 300))
    ruta_rank = Path(a.ranking) if a.ranking else dir_d / "ranking.csv"
    if not ruta_rank.exists():
        print(f"No existe {ruta_rank}: corra primero scripts/01_optimizar_malla.py", file=sys.stderr)
        return 2
    rk = exportar.leer_ranking(ruta_rank)
    val = _validacion(dir_d)
    N = a.n if a.n is not None else int(cfg.salida.get("N_planos", 5))
    fac = rk[rk["factible"]]
    hechos = []

    def dibujar(fila, prefijo, puesto=None, falla=""):
        cid = fila["cand_id"]
        pl = planos.plano(cfg, fila, dir_p / f"{prefijo}_{cid}", puesto=puesto, validacion=val.get(cid), falla=falla,
                          dpi=dpi)
        exportar.escribir_csv(exportar.perfil_candidato(cfg, fila), dir_d / f"{prefijo}_perfil.csv")
        hechos.append((prefijo, fila))
        print(f"  {pl.rutas[0].name} (+ .pdf){' · INFACTIBLE' if pl.infactible else ''}")
        return pl

    dir_p.mkdir(parents=True, exist_ok=True)
    if len(fac):
        print(f"Planos de los {min(N, len(fac))} mejores factibles → {dir_p}")
        for _, f in fac.head(N).iterrows():
            dibujar(f, f"rank{int(f['puesto']):02d}", f["puesto"])
    else:
        ruta_cf = dir_d / "casi_factibles.csv"
        if ruta_cf.exists():
            cf = pd.read_csv(ruta_cf)
        else:
            cf = frontera.casi_factibles(cfg, rk.assign(n_aletas=cfg.aleta.n))
        print(f"Sin factibles: planos de los {min(N_CASI_PLANOS, len(cf))} primeros casi factibles → {dir_p}")
        filas = rk.set_index("cand_id")
        for i, (_, c) in enumerate(cf.head(N_CASI_PLANOS).iterrows(), start=1):
            if c["cand_id"] in filas.index:
                f = filas.loc[c["cand_id"]].copy()
                f["cand_id"] = c["cand_id"]
            else:  # de la frontera (r_tip fuera de la malla del ranking): se reevalúa
                cu, al = specs_de_fila(c)
                f = pd.Series(evaluar_cuerpo(cfg, cu, [al])[0])
            falla = f"{c['restriccion_que_falla']} (SM máx. {c['SM_max_alcanzable_cal']:.2f} cal, faltan " \
                    f"{c['deficit_SM_cal']:.2f} cal)"
            dibujar(f, f"casi{i:02d}", None, falla)
    filas = rk.set_index("cand_id")
    for cid in a.cand:
        if cid not in filas.index:
            print(f"--cand {cid}: no está en {ruta_rank}", file=sys.stderr)
            return 2
        f = filas.loc[cid].copy()
        f["cand_id"] = cid
        puesto = f["puesto"] if pd.notna(f["puesto"]) else None
        prefijo = f"rank{int(puesto):02d}" if puesto is not None else "cand"
        dibujar(f, prefijo, puesto, "" if f["factible"] else frontera._falla(f["motivos"]))
    comp = [(p, f) for p, f in hechos if bool(f["factible"])]
    if comp:
        r = planos.comparativo(cfg, comp, dir_p / "comparativo_top.png")
        print(f"  {r.name}")
    print(f"Perfiles → {dir_d}/*_perfil.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
