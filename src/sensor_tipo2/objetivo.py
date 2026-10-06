"""Función objetivo ponderada, tolerancias implícitas y frente de Pareto.

    J = Σ w_j f_j,  f_j ∈ [0, 1] normalizados con cotas a priori de la malla configurada, en el
    orden de `objetivo.orden` (spec v2 rev. 2: masa, diámetro, SM, taper; el de dbf-sensor es
    masa, diámetro, taper, SM):

    masa      (m_max − m_total)/m_max
    diametro  (D_2 − D_lo)/(D_hi − D_lo),      D_2 = D aparente 2 max(R, r_tip) o D acostado
                                               max(D, √2 r_tip) (4 aletas a 45°), según objetivo.diametro
    SM        (SM_max − SM)/(SM_max − SM_min)  (modo 'max': se prefiere el mayor)
              |SM − SM_c| / ((SM_max − SM_min)/2)   (modo 'centro')
    taper     (k_ef − k_lo)/(1 − k_lo),        k_ef = sqrt(k² + f_b (1 − k²))

Tolerancia implícita: un criterio inferior solo compensa en J a lo sumo ε_j = Σ_{i>j} w_i / w_j
del rango del superior.
"""

from __future__ import annotations

import numpy as np

from .config import ConfigOpt

UNIDADES = {"masa": "g", "diametro": "mm", "SM": "cal", "taper": "(k efectivo)"}


def columnas(cfg: ConfigOpt) -> list[str]:
    """Nombres de las columnas f del ranking en el orden del objetivo: f1_masa, f2_diametro, …"""
    return [f"f{j + 1}_{c}" for j, c in enumerate(cfg.objetivo.orden)]


def col_criterio(cfg: ConfigOpt, criterio: str) -> str:
    return columnas(cfg)[cfg.objetivo.orden.index(criterio)]


def _escala_SM(cfg: ConfigOpt) -> float:
    r = cfg.restricciones
    return (r.SM_max - r.SM_min) / 2 if cfg.objetivo.f_SM_modo == "centro" else (r.SM_max - r.SM_min)


def f_valores(cfg: ConfigOpt, m_total, D_2, k, SM) -> np.ndarray:
    """Matriz (n, 4) de objetivos normalizados en el orden de `objetivo.orden`; cada columna
    recortada a [0, 1]. D_2 es el diámetro del criterio 2 (columna `cfg.objetivo.col_diametro`)."""
    c, r, ob = cfg.cotas, cfg.restricciones, cfg.objetivo
    m_total, D_2, k, SM = (np.asarray(v, dtype=float) for v in (m_total, D_2, k, SM))
    f = {"masa": (cfg.m_max - m_total) / cfg.m_max,
         "diametro": (D_2 - c["D_lo"]) / (c["D_hi"] - c["D_lo"]),
         "taper": (k - c["k_lo"]) / (c["k_hi"] - c["k_lo"]) if c["k_hi"] > c["k_lo"] else np.zeros_like(k),
         "SM": (np.abs(SM - ob.SM_centro) if ob.f_SM_modo == "centro" else r.SM_max - SM) / _escala_SM(cfg)}
    return np.clip(np.column_stack([f[n] for n in ob.orden]), 0.0, 1.0)


def J(cfg: ConfigOpt, F: np.ndarray) -> np.ndarray:
    return np.asarray(F, dtype=float) @ np.asarray(cfg.objetivo.pesos, dtype=float)


def tolerancias(cfg: ConfigOpt) -> dict:
    """ε_j adimensionales y en unidades físicas (g, mm, cal, k) del criterio en la posición j."""
    w = np.asarray(cfg.objetivo.pesos, dtype=float)
    orden = cfg.objetivo.orden
    eps = [float(w[j + 1:].sum() / w[j]) for j in range(4)]
    c = cfg.cotas
    rango = {"masa": cfg.m_max * 1e3, "diametro": (c["D_hi"] - c["D_lo"]) * 1e3, "SM": _escala_SM(cfg),
             "taper": c["k_hi"] - c["k_lo"]}
    fis = {n: eps[j] * rango[n] for j, n in enumerate(orden)}
    nombre = {"masa": "masa", "diametro": f"diámetro {cfg.objetivo.diametro}",
              "SM": "SM (" + ("más alto" if cfg.objetivo.f_SM_modo == "max" else "más cerca del objetivo") + ")",
              "taper": "taper"}
    textos = []
    for j, n in enumerate(orden[:-1]):
        iguales = " y ".join(nombre[o] for o in orden[:j])
        resto = ", ".join(nombre[o] for o in orden[j + 1:])
        fmt = {"masa": ".1f", "diametro": ".2f", "SM": ".3f", "taper": ".4f"}[n]
        textos.append((f"Con igual {iguales}, u" if j else "U")
                      + f"n candidato con hasta {fis[n]:{fmt}} {UNIDADES[n]} peor en {nombre[n]} puede ganar "
                        f"si es mejor en {resto}.")
    textos.append(f"{nombre[orden[-1]].capitalize()} es el último criterio: no compensa nada.")
    return {"orden": list(orden), "pesos": w.tolist(),
            "eps": {f"f{j + 1}_{n}": eps[j] for j, n in enumerate(orden)},
            "fisicas": {n: fis[n] for n in orden},
            "cotas": {"diametro": cfg.objetivo.diametro, "D_lo_mm": c["D_lo"] * 1e3, "D_hi_mm": c["D_hi"] * 1e3,
                      "k_lo": c["k_lo"], "k_hi": c["k_hi"], "m_max_g": cfg.m_max * 1e3, "SM_modo": cfg.objetivo.f_SM_modo},
            "lectura": textos}


def pareto(F: np.ndarray) -> np.ndarray:
    """Máscara de los no dominados (minimización en todas las columnas)."""
    F = np.asarray(F, dtype=float)
    n = len(F)
    if n == 0:
        return np.zeros(0, dtype=bool)
    nd = np.zeros(n, dtype=bool)
    frente = np.empty((0, F.shape[1]))
    for i in np.lexsort(F.T[::-1]):
        f = F[i]
        if frente.size and np.any(np.all(frente <= f, axis=1) & np.any(frente < f, axis=1)):
            continue
        nd[i] = True
        if not (frente.size and np.any(np.all(frente == f, axis=1))):
            frente = np.vstack([frente, f])
    return nd
