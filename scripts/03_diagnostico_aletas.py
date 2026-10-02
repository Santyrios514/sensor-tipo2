#!/usr/bin/env python
"""Script 3: ¿qué tan grandes tienen que ser las aletas? Barrido de r_tip/R sobre un cuerpo fijo.

    python scripts/03_diagnostico_aletas.py [--config config/optimizacion.yaml]
        [--cand CAND_ID --ranking data_opt/ranking.csv] [--r-max 3.0] [--paso 0.05]

Por defecto usa el cuerpo y la forma de aleta del `.ork` tipo 2 (D = 70, L_n = L_t = 65,
d_tc = 30, L_tc = 50, aleta 50/35/15 mm). Con --cand usa el cuerpo y la aleta de ese candidato.
Para cada r_tip/R reporta el C_Nα por componente, el CP, el SM máximo alcanzable con el tapón
delantero de plomo (SM_max_alc), el techo por geometría SM_∞ y el llenado completo (masa, SM,
tolerancia de amarre). También el C_Nα de cola + aletas por teoría de cuerpos esbeltos:

    C_Nα,cola+aletas = 2 [s_m² − r_tc² + r_tc⁴/s_m² − R²] / R²

que es ≤ 0 mientras s_m = r_tip ≤ R: con las aletas dentro del diámetro del cuerpo, la
transición y las aletas juntas no estabilizan.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from sensor_tipo2.aletas import barrowman, construir  # noqa: E402
from sensor_tipo2.barrido import specs_de_fila  # noqa: E402
from sensor_tipo2.config import AletaSpec, ConfigError, CuerpoSpec, cargar  # noqa: E402
from sensor_tipo2.exportar import AZUL, NARANJA, AQUA, TINTA2, REJILLA, escribir_csv, leer_ranking  # noqa: E402
from sensor_tipo2.geometria import construir_cuerpo  # noqa: E402
from sensor_tipo2.lastre import llenar  # noqa: E402

MM = 1e-3
BASE_ORK = (CuerpoSpec(400, 70, 65 / 70, "elipsoide", None, 65 / 70, 30 / 70, 50), AletaSpec(1.0, 0.7, 1.0, 1.0))


def _args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(RAIZ / "config" / "optimizacion.yaml"))
    p.add_argument("--cand", help="cand_id del ranking (por defecto, la geometría del .ork tipo 2)")
    p.add_argument("--ranking", help="por defecto, {salida.dir}/ranking.csv")
    p.add_argument("--r-max", type=float, default=3.0)
    p.add_argument("--paso", type=float, default=0.05)
    return p.parse_args(argv)


def barrer(cfg, c: CuerpoSpec, a: AletaSpec, r_rel: np.ndarray) -> pd.DataFrame:
    cu = construir_cuerpo(cfg, c)
    if cu.perfil is None:
        raise ValueError(f"cuerpo inválido: {cu.motivos}")
    p = cu.perfil
    A_ref = math.pi * p.D**2 / 4
    CN = {q.nombre: q for q in cu.partes}
    filas = []
    for r in r_rel:
        g, mot = construir(cfg, p, AletaSpec(a.mu_cr, a.gamma_ct, a.sigma_flecha, float(r)))
        if "h_menor_minimo" in mot:
            continue
        f = barrowman(g, A_ref, cfg.vuelo.mach, cfg.numerico.n_franjas)
        CNa = CN["nariz"].CNa + CN["cola"].CNa + f.CNa
        x_CP = (CN["nariz"].CNa * CN["nariz"].x + CN["cola"].CNa * CN["cola"].x + f.CNa * f.x_CP) / CNa
        s = g.r_tip
        esbelto = 2 * (s**2 - p.r_tc**2 + p.r_tc**4 / s**2 - p.R**2) / p.R**2
        Ll = llenar(cfg, cu, g.masa, g.x_cg, x_CP, CNa)
        filas.append({"r_tip_rel_R": float(r), "r_tip_mm": s / MM, "D_ap_mm": 2 * max(p.R, s) / MM,
                      "CNa_nariz": CN["nariz"].CNa, "CNa_cola": CN["cola"].CNa, "CNa_aletas": f.CNa,
                      "CNa_cola_mas_aletas": CN["cola"].CNa + f.CNa, "CNa_cola_aletas_esbelto": esbelto,
                      "CNa_total": CNa, "x_CP_mm": x_CP / MM, "SM_max_alc_cal": Ll.fila.get("SM_max_alcanzable_cal"),
                      "SM_inf_cal": Ll.fila.get("SM_inf_cal"), "factible": Ll.ok, "motivos": ";".join(Ll.motivos),
                      "m_total_g": Ll.fila.get("m_total_g"), "SM_cal": Ll.fila.get("SM_cal"),
                      "tol_amarre_mm": Ll.fila.get("tol_amarre_mm"), "restriccion_activa": Ll.restriccion_activa})
    return pd.DataFrame(filas)


def figura(df: pd.DataFrame, SM_min: float, titulo: str, ruta: Path):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    a1.plot(df["r_tip_rel_R"], df["CNa_aletas"], color=AZUL, label="aletas (Barrowman)")
    a1.plot(df["r_tip_rel_R"], df["CNa_cola_mas_aletas"], color=NARANJA, label="transición + aletas (Barrowman)")
    a1.plot(df["r_tip_rel_R"], df["CNa_cola_aletas_esbelto"], color=NARANJA, ls="--",
            label="transición + aletas (cuerpos esbeltos)")
    a1.axhline(0, color=TINTA2, lw=0.8)
    a1.axvline(1.0, color=TINTA2, lw=0.8, ls=":")
    a1.set_xlabel("r_tip / R")
    a1.set_ylabel("C_Nα [1/rad] (ref. π D²/4)")
    a1.set_title("Sustentación de la cola", loc="left", fontsize=10)
    a1.legend(frameon=False, fontsize=8)
    sm = df["SM_max_alc_cal"].astype(float).clip(lower=-6)
    a2.plot(df["r_tip_rel_R"], sm, color=AZUL, label="SM máx. con el tapón delantero de plomo")
    ok = df["factible"].astype(bool)
    a2.scatter(df.loc[ok, "r_tip_rel_R"], df.loc[ok, "SM_cal"], s=14, color=AQUA, zorder=3,
               label="llenado factible (SM final)")
    a2.axhline(SM_min, color=TINTA2, lw=1.0, ls=":")
    a2.axvline(1.0, color=TINTA2, lw=0.8, ls=":")
    a2.set_xlabel("r_tip / R")
    a2.set_ylabel("SM [cal]")
    a2.set_title("Margen estático", loc="left", fontsize=10)
    a2.legend(frameon=False, fontsize=8, loc="lower right")
    for ax in (a1, a2):
        ax.grid(color=REJILLA, lw=0.6)
    fig.suptitle(titulo, x=0.01, ha="left", fontsize=10)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta, dpi=140)
    plt.close(fig)


def main(argv=None) -> int:
    a = _args(argv)
    try:
        cfg = cargar(a.config)
    except ConfigError as e:
        print(e, file=sys.stderr)
        return 2
    if a.cand:
        df_r = leer_ranking(Path(a.ranking) if a.ranking else cfg.dir_salida() / "ranking.csv")
        c, al = specs_de_fila(df_r[df_r["cand_id"] == a.cand].iloc[0])
        nombre = a.cand
    else:
        c, al = BASE_ORK
        nombre = "geometría del .ork tipo 2"
    r_rel = np.round(np.arange(1.0, a.r_max + 1e-9, a.paso), 6)
    df = barrer(cfg, c, al, r_rel)
    SM_min = cfg.restricciones.SM_min
    print(f"{nombre}: {c.id}, aleta {al.id}")
    print(df[["r_tip_rel_R", "D_ap_mm", "CNa_cola_mas_aletas", "CNa_total", "x_CP_mm", "SM_max_alc_cal",
              "factible", "m_total_g", "SM_cal"]].iloc[::4].to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    est = df[df["SM_max_alc_cal"].astype(float) >= SM_min]
    fac = df[df["factible"].astype(bool)]
    for etq, d in (("SM ≥ SM_min con el tapón delantero", est), ("llenado factible (todos los criterios)", fac)):
        print(f"r_tip/R mínimo para {etq}: " + (f"{d['r_tip_rel_R'].min():.2f} (D_ap = {d['D_ap_mm'].min():.1f} mm)"
                                              if len(d) else f"ninguno hasta {a.r_max:g}"))
    base = df.iloc[0]
    print(f"Con r_tip = R: C_Nα transición + aletas = {base['CNa_cola_mas_aletas']:.3f} (Barrowman), "
          f"{base['CNa_cola_aletas_esbelto']:.3f} (cuerpos esbeltos); x_CP = {base['x_CP_mm']:.1f} mm.")
    escribir_csv(df, cfg.dir_salida() / "diagnostico_aletas.csv")
    figura(df, SM_min, f"{nombre} · {c.id} · aleta μ={al.mu_cr:g}, γ={al.gamma_ct:g}, σ={al.sigma_flecha:g}",
           cfg.dir_figuras() / "diagnostico_aletas.png")
    print(f"→ {cfg.dir_salida() / 'diagnostico_aletas.csv'}, {cfg.dir_figuras() / 'diagnostico_aletas.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
