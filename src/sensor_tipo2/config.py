"""Carga y validación de `config/optimizacion.yaml` y expansión de la malla.

La malla se separa en dos niveles, porque todo lo que no depende de las aletas (perfil, erosión,
volúmenes, masas del casco, CP de nariz y transición, región de lastre) se calcula una vez por
cuerpo:

- cuerpo: (L, D, L_n/D, forma de la transición, L_t/D, diámetro d_tc y largo L_tc del tubo de cola);
- aleta:  (μ, γ, σ, r_tip/R).

Todo lo que sale de aquí está en SI (m, kg, rad); los mm, g y grados del YAML se convierten aquí.
"""

from __future__ import annotations

import copy
import itertools
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .atmosfera import isa
from .perfiles import ESTACIONES, FORMAS, PARAMETRO_DEFECTO

MM = 1e-3
G = 1e-3

VARS_CUERPO = ("L_total_mm", "D_mm", "L_n_rel_D", "cola_forma", "L_t_rel_D", "d_tc_mm", "L_tc_mm")
VARS_ALETA = ("mu_cr", "gamma_ct", "sigma_flecha", "r_tip_rel_R")


class ConfigError(ValueError):
    """Error de configuración con todos los mensajes."""

    def __init__(self, errores: list[str]):
        self.errores = errores
        super().__init__("Configuración inválida:\n  - " + "\n  - ".join(errores))


def _num(v: float) -> str:
    return f"{v:g}"


# --------------------------------------------------------------------------- especificaciones de la malla


@dataclass(frozen=True)
class CuerpoSpec:
    L_mm: float
    D_mm: float
    L_n_rel_D: float
    forma: str
    parametro: float | None
    L_t_rel_D: float
    d_tc_mm: float  # diámetro del tubo de cola (absoluto: lo fijan el anclaje de aletas y la pared)
    L_tc_mm: float

    @property
    def id(self) -> str:
        p = "" if self.parametro is None else _num(self.parametro)
        return (f"L{_num(self.L_mm)}_D{_num(self.D_mm)}_n{_num(self.L_n_rel_D)}_{self.forma}{p}"
                f"_Lt{_num(self.L_t_rel_D)}_dtc{_num(self.d_tc_mm)}_tc{_num(self.L_tc_mm)}")

    @property
    def L(self) -> float:
        return self.L_mm * MM

    @property
    def D(self) -> float:
        return self.D_mm * MM

    @property
    def Ln(self) -> float:
        return self.L_n_rel_D * self.D

    @property
    def Lt(self) -> float:
        return self.L_t_rel_D * self.D

    @property
    def d_tc(self) -> float:
        return self.d_tc_mm * MM

    @property
    def k(self) -> float:
        """Razón de popa k = d_tc / D."""
        return self.d_tc_mm / self.D_mm

    @property
    def L_tc(self) -> float:
        return self.L_tc_mm * MM


@dataclass(frozen=True)
class AletaSpec:
    mu_cr: float
    gamma_ct: float
    sigma_flecha: float
    r_tip_rel_R: float

    @property
    def id(self) -> str:
        return (f"m{_num(self.mu_cr)}_g{_num(self.gamma_ct)}_s{_num(self.sigma_flecha)}"
                f"_r{_num(self.r_tip_rel_R)}")


def cand_id(c: CuerpoSpec, a: AletaSpec) -> str:
    return f"{c.id}_{a.id}"


# --------------------------------------------------------------------------- parámetros resueltos (SI)


@dataclass(frozen=True)
class Capa:
    material: str
    rho: float
    t: float
    phi: float = 1.0


def espesor_total(capas) -> float:
    return sum(c.t for c in capas)


def densidad_equivalente(capas) -> float:
    """ρ_eq = Σ ρ_k φ_k t_k / Σ t_k (una sola capa equivalente para OpenRocket)."""
    return sum(c.rho * c.phi * c.t for c in capas) / espesor_total(capas)


@dataclass(frozen=True)
class ParamsAleta:
    n: int
    rotacion: float  # rad
    material: str
    rho: float  # kg/m³ (sin fracción sólida)
    t: float
    phi: float
    E: float  # Pa
    nu: float
    factor_flutter: float


@dataclass(frozen=True)
class Electronica:
    Le: float
    me: float
    holgura: float
    x_cg_rel: float | None

    def x_cg(self, x_e):
        return x_e + (self.x_cg_rel if self.x_cg_rel is not None else self.Le / 2.0)


@dataclass(frozen=True)
class MasaPuntual:
    nombre: str
    m: float
    x: float


@dataclass(frozen=True)
class Lastre:
    rho_b: float
    fll: float
    r_min_util: float
    fmax: float | None
    margen_cola: float
    trasero: bool
    margen_popa: float
    costo_usd_kg: float = math.nan


@dataclass(frozen=True)
class Vuelo:
    V: float
    rho: float
    a: float
    altitud: float
    mach: float

    @property
    def q(self) -> float:
        return 0.5 * self.rho * self.V**2


@dataclass(frozen=True)
class Remolque:
    alpha_max: float  # rad
    alpha_lineal: float  # rad
    tol_min: float | None  # m


@dataclass(frozen=True)
class Restricciones:
    L_max: float
    SM_min: float
    SM_max: float
    k_min: float
    d_tc_min: float
    D_ap_max: float | None
    D_acostado_max_rel: float | None  # D_acostado ≤ este factor · D (1 = aletas dentro del D acostado)
    L_c_min_rel_cola: float | None  # convención del .ork: L_c ≥ este factor · (L_t + L_tc)
    h_min: float
    c_min: float
    eps_CN: float
    angulo_cola_max: float | None  # rad
    fineza_sin_arrastre: float = 3.0  # L_t/ΔD ≥ esto: la transición no suma arrastre de base (≈ 9.5°)
    fineza_base_roma: float = 1.0  # L_t/ΔD ≤ esto: equivale a una base roma (≈ 26.6°)
    base_roma_infactible: bool = True


@dataclass(frozen=True)
class Objetivo:
    pesos: tuple[float, float, float, float]
    f4_modo: str
    SM_centro: float
    diametro: str = "acostado"  # criterio 2: acostado | aparente

    @property
    def col_diametro(self) -> str:
        """Columna del ranking que usa el criterio 2."""
        return "D_acostado_mm" if self.diametro == "acostado" else "D_ap_mm"


@dataclass(frozen=True)
class Numerico:
    dx: float
    n_ell: int
    tol: float
    n_franjas: int
    tol_masa_max: float  # kg


@dataclass
class ConfigOpt:
    raw: dict[str, Any]
    ruta: Path | None
    m_max: float  # kg
    pared: dict[str, tuple[Capa, ...]]
    electronica: Electronica
    puntuales: tuple[MasaPuntual, ...]
    lastre: Lastre
    umbral_margen_bajo: float
    vuelo: Vuelo
    remolque: Remolque
    rot_guardado: float
    nariz: tuple[str, float | None]
    recortada: bool
    aleta: ParamsAleta
    malla: dict[str, list]
    restricciones: Restricciones
    objetivo: Objetivo
    numerico: Numerico
    ejecucion: dict[str, Any]
    openrocket: dict[str, Any]
    salida: dict[str, Any]
    raiz: Path = field(default_factory=Path.cwd)

    # ------------------------------------------------------------------ malla
    def cuerpos(self, malla: dict | None = None) -> list[CuerpoSpec]:
        m = malla or self.malla
        return [CuerpoSpec(float(L), float(D), float(n), f["forma"],
                           None if f.get("parametro") is None else float(f["parametro"]), float(lt),
                           float(dtc), float(tc))
                for L, D, n, f, lt, dtc, tc in itertools.product(*(m[v] for v in VARS_CUERPO))]

    def aletas(self, malla: dict | None = None) -> list[AletaSpec]:
        m = malla or self.malla
        return [AletaSpec(*(float(v) for v in combo)) for combo in itertools.product(*(m[k] for k in VARS_ALETA))]

    @property
    def n_evaluaciones(self) -> int:
        return len(self.cuerpos()) * len(self.aletas())

    @property
    def procesos(self) -> int:
        p = self.ejecucion.get("procesos", "auto")
        return max(1, os.cpu_count() or 1) if p in (None, "auto") else max(1, int(p))

    @property
    def cotas(self) -> dict[str, float]:
        """Cotas de normalización del objetivo, tomadas de la malla configurada (no de resultados).
        D_lo = min D; D_hi = el diámetro del criterio 2 con el mayor D y el mayor r_tip/R."""
        from .aletas import envolvente
        D = [float(d) * MM for d in self.malla["D_mm"]]
        R, r_tip = max(D) / 2, max(float(v) for v in self.malla["r_tip_rel_R"]) * max(D) / 2
        if self.objetivo.diametro == "acostado":
            D_hi = max(envolvente(r_tip, R, self.aleta.n, self.rot_guardado))
        else:
            D_hi = 2 * max(R, r_tip)
        k_lo = min(float(v) for v in self.malla["d_tc_mm"]) * MM / max(D)
        return {"D_lo": min(D), "D_hi": D_hi if D_hi > min(D) else min(D) + MM, "k_lo": k_lo, "k_hi": 1.0}

    def resolver(self, v: str | Path) -> Path:
        p = Path(v)
        return p if p.is_absolute() else self.raiz / p

    def dir_salida(self) -> Path:
        return self.resolver(self.salida.get("dir", "data_opt"))

    def dir_figuras(self) -> Path:
        return self.resolver(self.salida.get("dir_figuras", "figs_opt"))


# --------------------------------------------------------------------------- carga


def _material(nombre: str, mats: dict, donde: str, errores: list[str]) -> float:
    if nombre not in mats:
        errores.append(f"{donde}: material '{nombre}' no existe en 'materiales'")
        return math.nan
    return float(mats[nombre])


def _capas(lista, mats: dict, donde: str, errores: list[str]) -> tuple[Capa, ...]:
    if not lista:
        errores.append(f"{donde}: la pared necesita al menos una capa")
        return ()
    out = []
    for i, c in enumerate(lista):
        d = f"{donde}[{i}]"
        t = float(c.get("espesor_mm", 0.0)) * MM
        phi = float(c.get("fraccion_solida", 1.0))
        if not t > 0:
            errores.append(f"{d}: espesor_mm debe ser > 0")
        if not 0 < phi <= 1:
            errores.append(f"{d}: fraccion_solida debe estar en (0, 1]")
        out.append(Capa(c.get("material", "?"), _material(c.get("material", ""), mats, d, errores), t, phi))
    return tuple(out)


def _fineza(angulo_deg: float) -> float:
    """L_t/ΔD de una transición cónica equivalente con ese semiángulo: 1 / (2 tan θ)."""
    return 1.0 / (2.0 * math.tan(math.radians(angulo_deg)))


def _base_roma(b: dict, errores: list[str]) -> dict:
    a0 = float(b.get("angulo_inicio_deg", math.degrees(math.atan(1 / 6))))
    a1 = float(b.get("angulo_base_roma_deg", math.degrees(math.atan(0.5))))
    if not 0 < a0 < a1 < 90:
        errores.append("restricciones.cola_base_roma: se requiere 0 < angulo_inicio_deg < angulo_base_roma_deg < 90")
        return {}
    return {"fineza_sin_arrastre": _fineza(a0), "fineza_base_roma": _fineza(a1),
            "base_roma_infactible": bool(b.get("infactible", True))}


def _parametro(forma: str, p, donde: str, errores: list[str]) -> float | None:
    if forma not in FORMAS:
        errores.append(f"{donde}: forma '{forma}' no válida {FORMAS}")
        return None
    if p is None:
        if forma == "potencia":
            errores.append(f"{donde}: la forma 'potencia' necesita 'parametro'")
        return PARAMETRO_DEFECTO.get(forma)
    p = float(p)
    if forma in ("parabolica", "ogiva") and not 0 <= p <= 1:
        errores.append(f"{donde}: el parámetro de '{forma}' debe estar en [0, 1]")
    if forma == "potencia" and not p > 0:
        errores.append(f"{donde}: el exponente de 'potencia' debe ser > 0")
    return p


def deep_merge(base: dict, override: dict) -> dict:
    """Fusión en profundidad: los dict se combinan, lo demás (listas incluidas) se reemplaza."""
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _leer(path: Path) -> dict:
    """Lee un YAML; si tiene `hereda: otro.yaml` (relativo a su carpeta), lo fusiona encima de ese."""
    with open(path, encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}
    base = raw.pop("hereda", None)
    if base is None:
        return raw
    bp = Path(base)
    return deep_merge(_leer(bp if bp.is_absolute() else path.parent / bp), raw)


REQUERIDAS = ("materiales", "masa", "pared", "electronica", "lastre", "condiciones_vuelo", "remolque",
              "geometria_fija", "malla", "restricciones", "objetivo")


def cargar(ruta: str | Path | dict, raiz: Path | None = None) -> ConfigOpt:
    """Carga y valida optimizacion.yaml (o un dict). Lanza ConfigError con todos los problemas."""
    if isinstance(ruta, dict):
        raw, path = copy.deepcopy(ruta), None
        raiz = raiz or Path.cwd()
    else:
        path = Path(ruta).resolve()
        raw = _leer(path)
        raiz = raiz or path.parents[1]
    faltan = [s for s in REQUERIDAS if s not in raw]
    if faltan:
        raise ConfigError([f"falta la sección '{s}'" for s in faltan])
    errores: list[str] = []
    mats = raw["materiales"]

    ms = raw["masa"]
    m_max = float(ms["m_max_g"]) * G
    if not m_max > 0:
        errores.append("masa.m_max_g debe ser > 0")
    fll = float(ms.get("factor_llenado", 1.0))
    if not 0 < fll <= 1:
        errores.append("masa.factor_llenado debe estar en (0, 1]")

    p = raw["pared"]
    defecto = _capas(p.get("por_defecto"), mats, "pared.por_defecto", errores)
    pared = {e: (_capas(p[e], mats, f"pared.{e}", errores) if p.get(e) else defecto) for e in ESTACIONES}

    el = raw["electronica"]
    xr = el.get("x_cg_relativo_mm")
    electronica = Electronica(Le=float(el["longitud_mm"]) * MM, me=float(el["masa_g"]) * G,
                              holgura=float(el.get("holgura_mm", 0.0)) * MM,
                              x_cg_rel=None if xr is None else float(xr) * MM)
    if not electronica.Le > 0 or electronica.me < 0:
        errores.append("electronica: longitud_mm debe ser > 0 y masa_g ≥ 0")
    puntuales = tuple(MasaPuntual(m.get("nombre", f"p{i}"), float(m["masa_g"]) * G, float(m["x_mm"]) * MM)
                      for i, m in enumerate(raw.get("masas_puntuales") or []))

    la = raw["lastre"]
    fmax = la.get("fraccion_max_L")
    lastre = Lastre(rho_b=_material("plomo", mats, "lastre (plomo)", errores), fll=fll,
                    r_min_util=float(la.get("r_min_util_mm", 0.0)) * MM,
                    fmax=None if fmax is None else float(fmax),
                    margen_cola=float(la.get("margen_cola_mm", 0.0)) * MM,
                    trasero=bool(la.get("lastre_trasero", True)),
                    margen_popa=float(la.get("margen_popa_mm", 0.0)) * MM,
                    costo_usd_kg=float((raw.get("costos_usd_kg") or {}).get("plomo", math.nan)))

    cv = raw["condiciones_vuelo"]
    h = float(cv.get("altitud_m", 0.0))
    est = isa(h)
    V = float(cv["V_m_s"])
    mach = cv.get("mach", "auto")
    vuelo = Vuelo(V=V, rho=est.rho, a=est.a, altitud=h,
                  mach=V / est.a if mach in (None, "auto") else float(mach))

    rq = raw["remolque"]
    a_max = math.radians(float(rq.get("alpha_trim_max_deg", 5.0)))
    a_lin = math.radians(float(rq.get("alpha_lineal_max_deg", 15.0)))
    if not 0 < a_max <= a_lin:
        errores.append("remolque: se requiere 0 < alpha_trim_max_deg ≤ alpha_lineal_max_deg")
    tmin = rq.get("tol_amarre_min_mm")
    remolque = Remolque(alpha_max=a_max, alpha_lineal=a_lin, tol_min=None if tmin is None else float(tmin) * MM)

    gf = raw["geometria_fija"]
    nz = gf.get("nariz") or {}
    nariz = (nz.get("forma", "elipsoide"), _parametro(nz.get("forma", "elipsoide"), nz.get("parametro"),
                                                      "geometria_fija.nariz", errores))
    a = gf["aletas"]
    aleta = ParamsAleta(
        n=int(a.get("n", 4)), rotacion=math.radians(float(a.get("rotacion_deg", 0.0))),
        material=a.get("material", "?"), rho=_material(a.get("material", ""), mats, "geometria_fija.aletas", errores),
        t=float(a["espesor_mm"]) * MM, phi=float(a.get("fraccion_solida", 1.0)),
        E=float(a.get("E_GPa", math.nan)) * 1e9, nu=float(a.get("nu", math.nan)),
        factor_flutter=float(a.get("factor_flutter", 3.0)))
    if not aleta.t > 0:
        errores.append("geometria_fija.aletas.espesor_mm debe ser > 0")
    if not 0 < aleta.phi <= 1:
        errores.append("geometria_fija.aletas.fraccion_solida debe estar en (0, 1]")
    if aleta.n < 1:
        errores.append("geometria_fija.aletas.n debe ser ≥ 1")

    r = raw["restricciones"]
    dmax = r.get("D_ap_max_mm")
    rest = Restricciones(
        L_max=float(r.get("L_max_mm", 400.0)) * MM, SM_min=float(r["SM_min_cal"]), SM_max=float(r["SM_max_cal"]),
        k_min=float(r.get("k_min", 0.0)), d_tc_min=float(r.get("d_tc_min_mm", 0.0)) * MM,
        D_ap_max=None if dmax is None else float(dmax) * MM,
        D_acostado_max_rel=None if r.get("D_acostado_max_rel_D") is None else float(r["D_acostado_max_rel_D"]),
        L_c_min_rel_cola=None if r.get("L_c_min_rel_cola") is None else float(r["L_c_min_rel_cola"]),
        h_min=float(r.get("h_min_mm", 5.0)) * MM, c_min=float(r.get("c_min_mm", 5.0)) * MM,
        eps_CN=float(r.get("eps_CN", 0.5)),
        angulo_cola_max=None if r.get("angulo_cola_max_deg") is None else math.radians(float(r["angulo_cola_max_deg"])),
        **_base_roma(r.get("cola_base_roma") or {}, errores))
    if not rest.SM_max > rest.SM_min:
        errores.append("restricciones: SM_max_cal debe ser > SM_min_cal")
    if rest.D_acostado_max_rel is not None and not rest.D_acostado_max_rel >= 1.0:
        errores.append("restricciones.D_acostado_max_rel_D debe ser ≥ 1 (el D acostado nunca es menor que D)")

    m = raw["malla"]
    for k in VARS_CUERPO + VARS_ALETA:
        if not m.get(k):
            errores.append(f"malla.{k} está vacía o falta")
    if all(m.get(k) for k in VARS_CUERPO + VARS_ALETA):
        for f in m["cola_forma"]:
            _parametro(f.get("forma"), f.get("parametro"), "malla.cola_forma", errores)
        for v in m["L_total_mm"]:
            if not 0 < float(v) * MM <= rest.L_max + 1e-12:
                errores.append(f"malla.L_total_mm = {v} debe estar en (0, L_max = {rest.L_max / MM:g}]")
        for key in ("D_mm", "L_n_rel_D", "L_t_rel_D", "d_tc_mm", "L_tc_mm", "r_tip_rel_R"):
            for v in m[key]:
                if not float(v) > 0:
                    errores.append(f"malla.{key} = {v} debe ser > 0")
        for key in ("mu_cr", "gamma_ct"):
            for v in m[key]:
                if not 0 < float(v) <= 1:
                    errores.append(f"malla.{key} = {v} fuera de (0, 1]")
        for v in m["sigma_flecha"]:
            if not 0 <= float(v) <= 1:
                errores.append(f"malla.sigma_flecha = {v} fuera de [0, 1]")

    ob = raw["objetivo"]
    pesos = tuple(float(w) for w in ob.get("pesos", (1e6, 1e4, 1e2, 1.0)))
    if len(pesos) != 4 or not all(w > 0 for w in pesos):
        errores.append("objetivo.pesos: se esperan 4 pesos positivos")
    modo = ob.get("f4_modo", "centro")
    if modo not in ("centro", "max"):
        errores.append(f"objetivo.f4_modo '{modo}' no válido (centro | max)")
    diam = ob.get("diametro", "acostado")
    if diam not in ("acostado", "aparente"):
        errores.append(f"objetivo.diametro '{diam}' no válido (acostado | aparente)")

    nu = raw.get("numerico") or {}
    numerico = Numerico(dx=float(nu.get("dx_mm", 0.1)) * MM, n_ell=int(nu.get("n_barrido_ell", 600)),
                        tol=float(nu.get("tol_raiz_mm", 0.01)) * MM, n_franjas=int(nu.get("n_franjas_aleta", 48)),
                        tol_masa_max=float(nu.get("tol_masa_max_g", 0.5)) * G)
    if not numerico.dx > 0:
        errores.append("numerico.dx_mm debe ser > 0")
    if errores:
        raise ConfigError(errores)

    return ConfigOpt(
        raw=raw, ruta=path, m_max=m_max, pared=pared, electronica=electronica, puntuales=puntuales,
        lastre=lastre, umbral_margen_bajo=float((raw.get("estabilidad") or {}).get("umbral_margen_bajo_cal", 0.3)),
        vuelo=vuelo, remolque=remolque,
        rot_guardado=math.radians(float((raw.get("envolvente") or {}).get("rotacion_guardado_deg", 45.0))),
        nariz=nariz, recortada=bool(gf.get("cola_recortada", True)), aleta=aleta, malla=m, restricciones=rest,
        objetivo=Objetivo(pesos=pesos, f4_modo=modo, SM_centro=float(ob.get("SM_centro_cal", 1.5)), diametro=diam),
        numerico=numerico, ejecucion=raw.get("ejecucion") or {}, openrocket=raw.get("openrocket") or {},
        salida=raw.get("salida") or {}, raiz=raiz)
