"""Cuerpo del candidato: perfil, cavidad interior, masas del casco, CP de nariz y transición y
región de lastre. Todo lo que no depende de las aletas se calcula aquí una vez por cuerpo.

**Cavidad.** La erosión morfológica del meridiano {|y| ≤ r_e(x)} con un disco de radio T da la
superficie interior a espesor normal constante T:

    r_T(x) = max(0, min_{|ζ|≤T} [ r_e(x+ζ) − sqrt(T² − ζ²) ])

Es exacta para superficies de revolución y compone (erosionar por t1 y luego por t2 equivale a
erosionar por t1 + t2). Volumen y primer momento acumulados por Simpson:

    ∀(x) = ∫_0^x π r_i² dx,     Φ1(x) = ∫_0^x π r_i² x dx

**CP del cuerpo** (forma general de Barrowman, SymmetricComponentCalc de OpenRocket):

    C_Nα = 2 (A_a − A_f) / A_ref,     x_CP = x_0 + (ℓ A_a − ∀) / (A_a − A_f)

La nariz aporta +2 y la transición 2 (k² − 1) < 0; el cuerpo y el tubo de cola, 0.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import cumulative_simpson, simpson

from .config import ConfigOpt, CuerpoSpec, espesor_total
from .perfiles import ESTACIONES, Perfil

MM = 1e-3


# --------------------------------------------------------------------------- erosión y cavidad


def malla(L: float, dx: float) -> np.ndarray:
    return np.linspace(0.0, L, int(round(L / dx)) + 1)


def erosionar(r: np.ndarray, T: float, dx: float) -> np.ndarray:
    """Erosión del perfil r (malla uniforme de paso dx) con un disco de radio T.

    A proa r = 0 fuera de la malla (la nariz termina en punta); a popa el tubo de cola es abierto
    y se prolonga el último radio para que no aparezca pared de fondo.
    """
    if T <= 0:
        return r.copy()
    m = int(np.floor(T / dx + 1e-9))
    rp = np.concatenate([np.zeros(m), r, np.full(m, r[-1])])
    n = r.size
    out = np.full(n, np.inf)
    for j in range(-m, m + 1):
        z = j * dx
        np.minimum(out, rp[m + j: m + j + n] - math.sqrt(max(T * T - z * z, 0.0)), out=out)
    return np.maximum(out, 0.0)


@dataclass
class Cavidad:
    x: np.ndarray
    dx: float
    r_e: np.ndarray
    r_i: np.ndarray
    estacion: np.ndarray  # 0..3 por nodo
    fronteras: dict[str, list[np.ndarray]]  # estación -> [r_0 = r_e, r_1, ..., r_N = r_i]
    C0: np.ndarray = field(init=False)
    C1: np.ndarray = field(init=False)

    def __post_init__(self):
        A = np.pi * self.r_i**2
        self.C0 = cumulative_simpson(A, x=self.x, initial=0.0)
        self.C1 = cumulative_simpson(A * self.x, x=self.x, initial=0.0)

    def vol(self, a, b):
        """∀ entre a y b (np.interp satura fuera de [0, L])."""
        return np.interp(b, self.x, self.C0) - np.interp(a, self.x, self.C0)

    def mom(self, a, b):
        return np.interp(b, self.x, self.C1) - np.interp(a, self.x, self.C1)

    def area(self, x):
        return np.pi * np.interp(np.asarray(x, dtype=float), self.x, self.r_i) ** 2

    @property
    def V_ext(self) -> float:
        return float(simpson(np.pi * self.r_e**2, x=self.x))

    @property
    def V_int(self) -> float:
        return float(self.C0[-1])

    @property
    def x_fin(self) -> float:
        idx = np.nonzero(self.r_i > 0)[0]
        return float(self.x[idx[-1]]) if idx.size else 0.0


def construir_cavidad(perfil: Perfil, pared: dict, dx: float) -> Cavidad:
    """Erosión por estaciones con sus propias paredes. Cada nodo toma la erosión de su estación;
    en una ventana ±T_max alrededor de cada unión, la cavidad es el mínimo de ambas estaciones."""
    x = malla(perfil.L, dx)
    dx = float(x[1] - x[0])
    r_e = perfil.radio(x)
    est = perfil.estacion(x)
    E: dict[str, list[np.ndarray]] = {}
    cache: dict[float, np.ndarray] = {}
    for s in ESTACIONES:
        T = np.cumsum([c.t for c in pared[s]])
        E[s] = [r_e] + [cache.setdefault(round(float(Tk), 12), erosionar(r_e, float(Tk), dx)) for Tk in T]
    r_i = np.empty_like(r_e)
    for i, s in enumerate(ESTACIONES):
        r_i[est == i] = E[s][-1][est == i]
    T_max = max(espesor_total(pared[s]) for s in ESTACIONES)
    for i, xu in enumerate(perfil.uniones):
        s1, s2 = ESTACIONES[i], ESTACIONES[i + 1]
        w = np.abs(x - xu) <= T_max
        r_i[w] = np.minimum(E[s1][-1][w], E[s2][-1][w])
    fronteras = {s: [r_e] + [np.maximum(Ek, r_i) for Ek in E[s][1:-1]] + [r_i] for s in ESTACIONES}
    return Cavidad(x=x, dx=dx, r_e=r_e, r_i=r_i, estacion=est, fronteras=fronteras)


@dataclass(frozen=True)
class MasaCapa:
    estacion: str
    capa_idx: int  # 1 = exterior
    material: str
    m: float
    M: float  # primer momento m x̄


def masas_pared(cav: Cavidad, pared: dict) -> list[MasaCapa]:
    """m_k = ρ_k φ_k π ∫ (r_{k−1}² − r_k²) dx por estación, sin aproximación de pared delgada."""
    out = []
    for i, s in enumerate(ESTACIONES):
        mask = cav.estacion == i
        if not mask.any():
            continue
        fr = cav.fronteras[s]
        for k, capa in enumerate(pared[s], start=1):
            dA = np.pi * (fr[k - 1] ** 2 - fr[k] ** 2) * mask
            f = capa.rho * capa.phi
            out.append(MasaCapa(s, k, capa.material, f * float(simpson(dA, x=cav.x)),
                                f * float(simpson(dA * cav.x, x=cav.x))))
    return out


# --------------------------------------------------------------------------- Barrowman del cuerpo


@dataclass(frozen=True)
class Contribucion:
    nombre: str
    CNa: float
    x: float


def cuerpo_revolucion(x: np.ndarray, r: np.ndarray, x0: float, x1: float, A_ref: float,
                      nombre: str) -> Contribucion | None:
    m = (x >= x0 - 1e-12) & (x <= x1 + 1e-12)
    xs, rs = x[m], r[m]
    if xs.size < 3:
        return None
    A_f, A_a = np.pi * rs[0] ** 2, np.pi * rs[-1] ** 2
    if np.isclose(A_a, A_f, rtol=1e-12, atol=0):
        return None
    V = float(simpson(np.pi * rs**2, x=xs))
    ell = xs[-1] - xs[0]
    return Contribucion(nombre, 2.0 * (A_a - A_f) / A_ref, xs[0] + (ell * A_a - V) / (A_a - A_f))


# --------------------------------------------------------------------------- arrastre de la transición


def theta_eq(R: float, k: float, L_t: float) -> float:
    """Ángulo equivalente de la transición, arctan(R (1 − k) / L_t) [rad] (informativo)."""
    return math.atan2(R * (1.0 - k), L_t)


def fineza_cola(D: float, k: float, L_t: float) -> float:
    """L_t / ΔD con ΔD = D (1 − k): la 'fineness' que usa OpenRocket para el arrastre de base."""
    dD = D * (1.0 - k)
    return L_t / dD if dD > 0 else math.inf


def fraccion_base_roma(fineza: float, sin_arrastre: float = 3.0, base_roma: float = 1.0) -> float:
    """Fracción f_b del arrastre de base que suma la transición (OpenRocket 24.12, Hoerner):
    0 si L_t/ΔD ≥ 3, 1 si ≤ 1 y lineal entre ambos."""
    return min(max((sin_arrastre - fineza) / (sin_arrastre - base_roma), 0.0), 1.0)


def k_efectivo(k: float, f_b: float) -> float:
    """k_ef = sqrt(k² + f_b (1 − k²)): razón de popa de una base con el mismo arrastre."""
    return math.sqrt(k * k + f_b * (1.0 - k * k))


# --------------------------------------------------------------------------- cuerpo (caché)


@dataclass(frozen=True)
class LimiteGeo:
    ell_geo: float
    limitante: str
    candidatos: dict[str, float]
    x_e_a: float = -math.inf  # tramo admisible de la electrónica [x_a, x_b] (con electronica.r_min_mm)
    x_e_b: float = math.inf


def tramo_electronica(cav: Cavidad, r_min: float) -> tuple[float, float] | None:
    """[x_a, x_b]: el tramo contiguo más largo donde r_i ≥ r_min (un intervalo si el perfil es
    unimodal). None si no hay ninguno."""
    ok = cav.r_i >= r_min
    if not ok.any():
        return None
    d = np.diff(np.r_[0, ok.astype(np.int8), 0])
    ini, fin = np.nonzero(d == 1)[0], np.nonzero(d == -1)[0] - 1
    i = int(np.argmax(cav.x[fin] - cav.x[ini]))
    return float(cav.x[ini[i]]), float(cav.x[fin[i]])


def x_inicio_lastre(cfg: ConfigOpt, cav: Cavidad) -> float:
    """x_b0: primer x donde r_i ≥ r_min,util."""
    idx = np.nonzero(cav.r_i >= cfg.lastre.r_min_util)[0]
    return float(cav.x[idx[0]]) if idx.size else math.nan


def limite_geometrico(cfg: ConfigOpt, perfil: Perfil, cav: Cavidad, x_b0: float) -> LimiteGeo:
    """ℓ_geo del tapón delantero: fin de la cavidad, inicio de la transición (menos margen),
    electrónica detrás del lastre y fracción máxima de L."""
    if not math.isfinite(x_b0):
        return LimiteGeo(0.0, "sin_cavidad_util", {})
    la, el = cfg.lastre, cfg.electronica
    if el.r_min is not None:
        # la electrónica va donde r_i ≥ r_min (nariz, cuerpo o transición); el tapón delantero la
        # precede y puede entrar en la transición: ℓ_geo = x_b − L_e − 2 holgura − x_b0
        tramo = tramo_electronica(cav, el.r_min)
        if tramo is None or tramo[1] - tramo[0] < el.Le:
            return LimiteGeo(0.0, "electronica_no_cabe", {})
        x_a, x_b = tramo
        cand = {"fin_cavidad": cav.x_fin - x_b0, "electronica_no_cabe": x_b - el.Le - 2 * el.holgura - x_b0}
        if la.fmax is not None:
            cand["fraccion_max_L"] = la.fmax * perfil.L - x_b0
        lim = min(cand, key=cand.get)
        return LimiteGeo(cand[lim], lim, cand, x_a, x_b)
    x_cola = perfil.x_t0 - la.margen_cola
    cand = {"fin_cavidad": cav.x_fin - x_b0, "cola": max(x_cola - x_b0, 0.0),
            "electronica_detras": x_cola - el.Le - el.holgura - x_b0}
    if la.fmax is not None:
        cand["fraccion_max_L"] = la.fmax * perfil.L - x_b0
    lim = min(cand, key=cand.get)
    return LimiteGeo(cand[lim], lim, cand)


@dataclass
class Cuerpo:
    """Caché por cuerpo. `motivos` no vacío = cuerpo descartado."""

    spec: CuerpoSpec
    motivos: list[str] = field(default_factory=list)
    perfil: Perfil | None = None
    cav: Cavidad | None = None
    capas: list[MasaCapa] = field(default_factory=list)
    m_puntuales: float = 0.0
    M_puntuales: float = 0.0
    x_b0: float = math.nan
    lim: LimiteGeo | None = None
    partes: tuple[Contribucion, ...] = ()  # nariz y transición
    theta_eq: float = math.nan
    fineza: float = math.nan
    f_base_roma: float = math.nan
    k_ef: float = math.nan

    @property
    def ok(self) -> bool:
        return not self.motivos

    @property
    def m_casco(self) -> float:
        return sum(c.m for c in self.capas)

    @property
    def M_casco(self) -> float:
        return sum(c.M for c in self.capas)

    def masas_por_estacion(self) -> dict[str, tuple[float, float]]:
        out = {}
        for s in ESTACIONES:
            cs = [c for c in self.capas if c.estacion == s]
            m = sum(c.m for c in cs)
            out[s] = (m, sum(c.M for c in cs) / m if m > 0 else math.nan)
        return out


def perfil_de(cfg: ConfigOpt, spec: CuerpoSpec) -> Perfil:
    return Perfil(L=spec.L, D=spec.D, Ln=spec.Ln, forma_n=cfg.nariz[0], param_n=cfg.nariz[1], Lt=spec.Lt,
                  forma_t=spec.forma, param_t=spec.parametro, recortada=cfg.recortada, d_tc=spec.d_tc,
                  L_tc=spec.L_tc)


def motivos_previos(cfg: ConfigOpt, spec: CuerpoSpec) -> list[str]:
    """Motivos de descarte que no necesitan la cavidad (sirven de prefiltro de la malla)."""
    cu = Cuerpo(spec=spec)
    _chequeos_previos(cfg, cu)
    return cu.motivos


def _chequeos_previos(cfg: ConfigOpt, cu: Cuerpo):
    spec, rest, el = cu.spec, cfg.restricciones, cfg.electronica
    R = spec.D / 2
    Lt = max(spec.Lt, 0.0)
    cu.theta_eq = theta_eq(R, spec.k, Lt)
    cu.fineza = fineza_cola(spec.D, spec.k, Lt)
    cu.f_base_roma = fraccion_base_roma(cu.fineza, rest.fineza_sin_arrastre, rest.fineza_base_roma)
    cu.k_ef = k_efectivo(spec.k, cu.f_base_roma)
    if spec.L > rest.L_max + 1e-12:
        cu.motivos.append("L_mayor_maximo")
    if spec.k < rest.k_min - 1e-12:
        cu.motivos.append("k_bajo_minimo")
    if not spec.d_tc < spec.D:
        cu.motivos.append("d_tubo_cola_no_menor_que_D")
    if spec.d_tc < rest.d_tc_min - 1e-12:
        cu.motivos.append("d_tubo_cola_bajo_minimo")
    if rest.angulo_cola_max is not None and cu.theta_eq > rest.angulo_cola_max + 1e-12:
        cu.motivos.append("angulo_cola")
    if rest.base_roma_infactible and cu.fineza <= rest.fineza_base_roma + 1e-12:
        cu.motivos.append("cola_base_roma")
    T_tc = espesor_total(cfg.pared["tubo_cola"])
    if not T_tc < spec.d_tc / 2:
        cu.motivos.append("pared_mayor_que_tubo_cola")
    if spec.L_n_rel_D < rest.L_n_rel_D_min - 1e-12:
        cu.motivos.append("nariz_bajo_minimo")
    if not spec.L_disp > 0:
        cu.motivos.append("L_disp_no_positivo")  # la nariz y el tubo de cola no dejan lugar a la transición
    elif rest.L_c_min_rel_cola is not None and spec.L_c < rest.L_c_min_rel_cola * (spec.Lt + spec.L_tc) - 1e-12:
        cu.motivos.append("cuerpo_central_corto")
    if el.r_min is None and spec.L_disp > 0 and \
            spec.Ln + spec.L_c - cfg.lastre.margen_cola - el.Le - el.holgura <= 0:
        cu.motivos.append("electronica_no_cabe")  # ℓ_geo ≤ −x_b0 ≤ 0 (lo mismo que daría limite_geometrico)


def construir_cuerpo(cfg: ConfigOpt, spec: CuerpoSpec) -> Cuerpo:
    cu = Cuerpo(spec=spec)
    _chequeos_previos(cfg, cu)
    if cu.motivos:
        return cu
    cu.perfil = perfil_de(cfg, spec)
    cu.cav = construir_cavidad(cu.perfil, cfg.pared, cfg.numerico.dx)
    cu.capas = masas_pared(cu.cav, cfg.pared)
    cu.m_puntuales = sum(p.m for p in cfg.puntuales)
    cu.M_puntuales = sum(p.m * p.x for p in cfg.puntuales)
    cu.x_b0 = x_inicio_lastre(cfg, cu.cav)
    cu.lim = limite_geometrico(cfg, cu.perfil, cu.cav, cu.x_b0)
    if not cu.lim.ell_geo > 0:
        cu.motivos.append("electronica_no_cabe" if cu.lim.limitante in ("electronica_detras", "electronica_no_cabe")
                          else "sin_region_lastre")
    A_ref = math.pi * spec.D**2 / 4
    p = cu.perfil
    partes = (cuerpo_revolucion(cu.cav.x, cu.cav.r_e, 0.0, p.Ln, A_ref, "nariz"),
              cuerpo_revolucion(cu.cav.x, cu.cav.r_e, p.x_t0, p.x_tc0, A_ref, "cola"))
    cu.partes = tuple(c for c in partes if c is not None)
    return cu
