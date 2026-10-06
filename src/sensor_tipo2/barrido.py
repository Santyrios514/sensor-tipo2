"""Barrido en malla con caché por cuerpo, en paralelo; ranking, Pareto y refinamiento.

Unidad de trabajo: un cuerpo con su lista de aletas. El cuerpo (perfil, cavidad, masas del
casco, x_b0, ℓ_geo, CP de nariz y transición) se construye una vez; cada aleta agrega su
polígono, su masa, su C_Nα y el llenado. Las filas vuelven en el orden pedido, así que el
ranking es reproducible sin importar el número de procesos.
"""

from __future__ import annotations

import itertools
import math
import multiprocessing as mp
from collections import OrderedDict

import numpy as np
import pandas as pd

from .aletas import construir, envolvente, velocidad_flutter
from .config import AletaSpec, ConfigOpt, CuerpoSpec, cand_id
from .geometria import construir_cuerpo, motivos_previos
from .lastre import llenar
from .objetivo import J, col_criterio, columnas, f_valores, pareto
from .sustituto import cp_sustituto

MM, G = 1e-3, 1e-3

RANKING = [
    "puesto", "cand_id", "factible", "motivos", "J", "f1_masa", "f2_diametro", "f3_SM", "f4_taper",
    # cuerpo
    "L_mm", "D_mm", "L_n_mm", "f_cil", "L_c_mm", "cola_forma", "cola_parametro", "L_t_mm", "k", "d_tc_mm",
    "L_tc_mm", "theta_eq_deg", "fineza_cola", "f_base_roma", "k_efectivo", "esbeltez_tubo", "tubo_esbelto",
    "aletas_en_estela",
    # aletas
    "x_LE_mm", "c_r_mm", "c_t_mm", "x_s_mm", "h_mm", "r_tip_mm", "D_ap_mm", "D_acostado_mm",
    "A_aleta_mm2", "AR_aleta", "frac_h_fuera_sombra", "m_aletas_g", "V_flutter_m_s", "flutter_margen_bajo",
    # masa, lastre y estabilidad
    "m_total_g", "m_lastre_g", "m_lastre_delantero_g", "m_lastre_trasero_g", "ell_mm", "ell_trasero_mm",
    "restriccion_activa", "x_CG_mm", "x_CP_sust_mm", "CN_alpha_total", "CN_alpha_aletas",
    "CN_alpha_cola", "SM_cal", "SM_inf_cal", "SM_max_alcanzable_cal", "tol_amarre_mm", "alpha_trim_deg",
    "m_en_SM_max_g", "tol_en_SM_max_mm", "banderas_lastre", "en_pareto",
    # malla y otros
    "cuerpo_id", "refinamiento", "L_n_rel_D", "mu_cr", "gamma_ct", "sigma_flecha", "r_tip_rel_R",
    "x_b0_mm", "ell_geo_mm", "m_casco_g", "x_electronica_mm", "costo_lastre_usd",
]


# --------------------------------------------------------------------------- evaluación


def _fila_base(c: CuerpoSpec, a: AletaSpec) -> dict:
    return {"cand_id": cand_id(c, a), "cuerpo_id": c.id, "factible": False, "motivos": "",
            "L_mm": c.L_mm, "D_mm": c.D_mm, "L_n_rel_D": c.L_n_rel_D, "L_n_mm": c.Ln / MM, "cola_forma": c.forma,
            "cola_parametro": c.parametro, "f_cil": c.f_cil, "L_t_mm": c.Lt / MM, "k": c.k,
            "d_tc_mm": c.d_tc_mm, "L_tc_mm": c.L_tc_mm, "L_c_mm": c.L_c / MM, "esbeltez_tubo": c.L_tc_mm / c.d_tc_mm,
            "mu_cr": a.mu_cr, "gamma_ct": a.gamma_ct, "sigma_flecha": a.sigma_flecha, "r_tip_rel_R": a.r_tip_rel_R,
            "en_pareto": False}


def evaluar_cuerpo(cfg: ConfigOpt, c: CuerpoSpec, aletas: list[AletaSpec], cu=None) -> list[dict]:
    """Filas del ranking para todas las aletas de un cuerpo (factibles o no, con motivos). `cu`
    reutiliza un cuerpo ya construido (la geometría del cuerpo no depende de las aletas)."""
    cu = cu if cu is not None else construir_cuerpo(cfg, c)
    rest = cfg.restricciones
    comunes = {"theta_eq_deg": math.degrees(cu.theta_eq), "fineza_cola": cu.fineza, "f_base_roma": cu.f_base_roma,
               "k_efectivo": cu.k_ef,
               "tubo_esbelto": rest.esbeltez_tubo_aviso is not None and c.L_tc_mm / c.d_tc_mm > rest.esbeltez_tubo_aviso,
               "aletas_en_estela": rest.angulo_estela_aviso is not None and cu.theta_eq > rest.angulo_estela_aviso}
    if not cu.ok:
        return [{**_fila_base(c, a), **comunes, "motivos": ";".join(cu.motivos)} for a in aletas]
    vu = cfg.vuelo
    R = c.D / 2
    comunes.update({"x_b0_mm": cu.x_b0 / MM, "ell_geo_mm": cu.lim.ell_geo / MM, "m_casco_g": cu.m_casco / G,
                    "CN_alpha_cola": next((p.CNa for p in cu.partes if p.nombre == "cola"), 0.0)})
    filas = []
    for a in aletas:
        f = {**_fila_base(c, a), **comunes}
        g, motivos = construir(cfg, cu.perfil, a)
        f.update({"x_LE_mm": g.x_LE / MM, "c_r_mm": g.c_r / MM, "c_t_mm": g.c_t / MM, "x_s_mm": g.x_s / MM,
                  "h_mm": g.h / MM, "r_tip_mm": g.r_tip / MM, "D_ap_mm": 2 * max(R, g.r_tip) / MM,
                  "frac_h_fuera_sombra": max(0.0, g.r_tip - R) / g.h if g.h > 0 else math.nan,
                  "D_acostado_mm": max(envolvente(g.r_tip, R, cfg.aleta.n, cfg.rot_guardado)) / MM})
        if motivos:
            f["motivos"] = ";".join(motivos)
            filas.append(f)
            continue
        V_f = velocidad_flutter(g, vu.altitud)
        f.update({"A_aleta_mm2": g.area / MM**2, "AR_aleta": g.AR, "m_aletas_g": g.masa / G, "V_flutter_m_s": V_f,
                  "flutter_margen_bajo": bool(V_f < cfg.aleta.factor_flutter * vu.V)})
        cp, motivo = cp_sustituto(cu, g, vu.mach, cfg.numerico.n_franjas, rest.eps_CN)
        if cp is None:
            f["motivos"] = motivo
            filas.append(f)
            continue
        Ll = llenar(cfg, cu, g.masa, g.x_cg, cp.x_CP, cp.CNa)
        f.update({"x_CP_sust_mm": cp.x_CP / MM, "CN_alpha_total": cp.CNa,
                  "CN_alpha_aletas": cp.partes[-1].CNa, "restriccion_activa": Ll.restriccion_activa,
                  "banderas_lastre": ";".join(Ll.banderas), **Ll.fila})
        f["motivos"] = ";".join(Ll.motivos) if Ll.motivos else ("" if Ll.ok else "sin_llenado")
        f["factible"] = Ll.ok
        filas.append(f)
    return filas


_CFG: ConfigOpt | None = None


def _init(cfg: ConfigOpt):
    global _CFG
    _CFG = cfg


def _tarea(args):
    c, aletas = args
    return evaluar_cuerpo(_CFG, c, aletas)


Grupos = "OrderedDict[CuerpoSpec, list[AletaSpec]]"


def evaluar(cfg: ConfigOpt, grupos: Grupos, procesos: int | None = None, progreso=None) -> pd.DataFrame:
    """Evalúa {cuerpo: [aletas]} y devuelve las filas en el mismo orden, con f, J y el ranking."""
    procesos = procesos or cfg.procesos
    tareas = list(grupos.items())
    filas: list[dict] = []
    if procesos <= 1 or len(tareas) <= 1:
        _init(cfg)
        for i, t in enumerate(tareas):
            filas += _tarea(t)
            if progreso:
                progreso(i + 1, len(tareas))
    else:
        with mp.get_context("fork").Pool(procesos, initializer=_init, initargs=(cfg,)) as pool:
            for i, r in enumerate(pool.imap(_tarea, tareas, chunksize=4)):
                filas += r
                if progreso:
                    progreso(i + 1, len(tareas))
    return completar(cfg, pd.DataFrame(filas))


def completar(cfg: ConfigOpt, df: pd.DataFrame) -> pd.DataFrame:
    """Objetivos, J, puesto y frente de Pareto; columnas en el orden del esquema."""
    df = df.reindex(columns=list(dict.fromkeys(RANKING + list(df.columns))))
    df["factible"] = df["factible"].fillna(False).astype(bool)
    m_total = df["m_total_g"].astype(float) * G
    SM = df["SM_cal"].astype(float)
    F = f_valores(cfg, m_total, df[cfg.objetivo.col_diametro].astype(float) * MM, df["k_efectivo"].astype(float), SM)
    valido = np.isfinite(m_total) & np.isfinite(SM)
    cols = columnas(cfg)
    for j, col in enumerate(cols):
        df[col] = np.where(valido, F[:, j], np.nan)
    df["J"] = np.where(valido, J(cfg, np.nan_to_num(F)), np.nan)
    df["en_pareto"] = False
    df = ordenar(df)
    ok = df["factible"]
    if ok.any():
        fac = df[ok]
        tol = cfg.numerico.tol_masa_max
        alcanza = fac["m_total_g"] * G >= cfg.m_max - tol
        fm = col_criterio(cfg, "masa")
        base = fac[alcanza] if alcanza.any() else fac[fac[fm] <= fac[fm].min() + tol / cfg.m_max]
        df.loc[base.index[pareto(base[[c for c in cols if c != fm]].to_numpy())], "en_pareto"] = True
    return df


def ordenar(df: pd.DataFrame) -> pd.DataFrame:
    """Factibles por J (desempate por cand_id), luego los infactibles; puesto = 1..n factibles."""
    df = df.assign(_inf=~df["factible"]).sort_values(["_inf", "J", "cand_id"], na_position="last",
                                                     kind="stable").drop(columns="_inf").reset_index(drop=True)
    n_ok = int(df["factible"].sum())
    df["puesto"] = np.where(np.arange(len(df)) < n_ok, np.arange(1, len(df) + 1), np.nan)
    return df


# --------------------------------------------------------------------------- malla y refinamiento


def aleta_imposible(cfg: ConfigOpt, a: AletaSpec, L_tc: float) -> bool:
    """Aleta descartable sin construir el cuerpo: cuerda de raíz bajo el mínimo o r_tip sobre el tope."""
    rest = cfg.restricciones
    return (a.mu_cr * L_tc < max(rest.c_min, rest.c_r_min) - 1e-12
            or (rest.r_tip_rel_R_max is not None and a.r_tip_rel_R > rest.r_tip_rel_R_max + 1e-12))


def grupos_malla(cfg: ConfigOpt, malla: dict | None = None, prefiltrar: bool = True) -> Grupos:
    """{cuerpo: [aletas]}. Con `prefiltrar`, quita los cuerpos geométricamente imposibles
    (motivos_previos) y las aletas imposibles por su cuerda o su r_tip, antes de evaluar."""
    aletas = cfg.aletas(malla)
    out: OrderedDict = OrderedDict()
    for c in cfg.cuerpos(malla):
        if prefiltrar and motivos_previos(cfg, c):
            continue
        al = [a for a in aletas if not (prefiltrar and aleta_imposible(cfg, a, c.L_tc))]
        if al:
            out[c] = al
    return out


def n_candidatos(grupos: Grupos) -> int:
    return sum(len(a) for a in grupos.values())


def specs_de_fila(r) -> tuple[CuerpoSpec, AletaSpec]:
    get = r.get if hasattr(r, "get") else (lambda k: getattr(r, k))
    par = get("cola_parametro")
    par = None if par is None or (isinstance(par, float) and math.isnan(par)) else float(par)
    c = CuerpoSpec(float(get("L_mm")), float(get("D_mm")), float(get("L_n_rel_D")), str(get("cola_forma")), par,
                   float(get("f_cil")), float(get("d_tc_mm")), float(get("L_tc_mm")))
    a = AletaSpec(float(get("mu_cr")), float(get("gamma_ct")), float(get("sigma_flecha")), float(get("r_tip_rel_R")))
    return c, a


def _vecinos(v: float, lista, lo: float, hi: float) -> list[float]:
    """v y los puntos a medio paso hacia sus vecinos de la malla original, dentro de [lo, hi]."""
    xs = sorted(set(float(x) for x in lista))
    out = {v}
    menores = [x for x in xs if x < v - 1e-12]
    mayores = [x for x in xs if x > v + 1e-12]
    if menores:
        out.add(v - (v - menores[-1]) / 2)
    if mayores:
        out.add(v + (mayores[0] - v) / 2)
    return sorted(round(x, 10) for x in out if lo - 1e-12 <= x <= hi + 1e-12)


# Spec v2 rev. 2 §3.1/§5: el refinamiento solo toca los grupos A (masa), B (diámetro aparente) y
# C (estabilidad); las variables de ajuste fino (μ, γ, σ, forma de la transición) y L quedan fijas.
REFINABLES_CUERPO = ("D_mm", "L_n_rel_D", "f_cil", "tubo", "L_tc_mm")  # tubo = k o d_tc_mm
REFINABLES_ALETA = ("r_tip_rel_R",)


def grupos_refinamiento(cfg: ConfigOpt, ranking: pd.DataFrame, top_K: int) -> Grupos:
    """Malla a medio paso alrededor de los top_K factibles, sin repetir candidatos, en las variables
    de los grupos A, B y C (REFINABLES_*). El tubo se refina en la variable de la malla (k o d_tc_mm)."""
    m, rest = cfg.malla, cfg.restricciones
    vt = cfg.var_tubo
    vars_c = tuple(vt if v == "tubo" else v for v in REFINABLES_CUERPO)
    refinables = vars_c + REFINABLES_ALETA
    fac = ranking[ranking["factible"]].head(top_K)
    vistos = set(ranking["cand_id"])
    limites = {k: (min(float(v) for v in m[k]), max(float(v) for v in m[k])) for k in refinables}
    if vt == "k":
        limites["k"] = (max(limites["k"][0], rest.k_min), limites["k"][1])
    else:
        limites["d_tc_mm"] = (max(limites["d_tc_mm"][0], rest.d_tc_min / MM), limites["d_tc_mm"][1])
    if rest.r_tip_rel_R_max is not None:
        limites["r_tip_rel_R"] = (limites["r_tip_rel_R"][0], min(limites["r_tip_rel_R"][1], rest.r_tip_rel_R_max))
    out: OrderedDict = OrderedDict()
    for _, r in fac.iterrows():
        c0, a0 = specs_de_fila(r)
        vals = {k: _vecinos(float(r[k]), m[k], *limites[k]) for k in refinables}
        for D, n, fc, tubo, tc in itertools.product(*(vals[v] for v in vars_c)):
            d_tc = round(tubo * D, 9) if vt == "k" else tubo
            c = CuerpoSpec(c0.L_mm, D, n, c0.forma, c0.parametro, fc, d_tc, tc)
            if motivos_previos(cfg, c):
                continue
            lista = out.setdefault(c, [])
            ya = {a.id for a in lista}
            for rt in vals["r_tip_rel_R"]:
                a = AletaSpec(a0.mu_cr, a0.gamma_ct, a0.sigma_flecha, rt)
                if a.id in ya or cand_id(c, a) in vistos or aleta_imposible(cfg, a, c.L_tc):
                    continue
                ya.add(a.id)
                lista.append(a)
    return OrderedDict((c, a) for c, a in out.items() if a)


def optimizar(cfg: ConfigOpt, refinar: bool = True, procesos: int | None = None, progreso=None,
              grupos: Grupos | None = None) -> pd.DataFrame:
    """Barrido completo (+ refinamiento opcional) → ranking con los candidatos evaluados (los
    descartados por el prefiltro no aparecen)."""
    df = evaluar(cfg, grupos if grupos is not None else grupos_malla(cfg), procesos, progreso)
    df = df.assign(refinamiento=False)
    ref = cfg.ejecucion.get("refinamiento") or {}
    if refinar and ref.get("activar", True):
        gr = grupos_refinamiento(cfg, df, int(ref.get("top_K", 10)))
        if gr:
            nuevo = evaluar(cfg, gr, procesos, progreso).assign(refinamiento=True)
            df = completar(cfg, pd.concat([df, nuevo], ignore_index=True))
    return df
