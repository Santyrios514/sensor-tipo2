"""Validación de los ganadores con OpenRocket 24.12 (orlab + JPype). Único módulo con Java.

`PuenteTipo2` abre `modelos/analisis_tipo2.ork` como plantilla, ubica sus cinco componentes por
tipo y orden (los dos tubos se llaman igual, "Body Tube"), fija la geometría del candidato (nariz,
cuerpo, transición, tubo de cola, paredes con ρ_eq y aletas freeform sobre el tubo con offset
*bottom* 0) y agrega el lastre, la electrónica y las masas puntuales como *Mass components*.

Flujo (spec v2 rev. 2 §7, `validar_ganadores`): el ranking sale completo del modelo propio; aquí
solo se arman en OpenRocket los N_ganadores mejores factibles y se comparan CP, C_Nα, masa y CG.
Criterio por candidato: |Δx_CP| ≤ tol_cp_mm, |Δm| ≤ tol_masa_pct y SM_OR ∈ [SM_min, SM_max]. Si
alguno supera tol_cp_mm, el modelo propio dejó de coincidir con OpenRocket para esa familia de
formas: se reporta y no se declara ganador (no hay calibración). No se guarda ningún `.ork`.
"""

from __future__ import annotations

import contextlib
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import simpson

from .barrido import specs_de_fila
from .config import ConfigOpt, densidad_equivalente, espesor_total
from .exportar import Detalle, detalle, masas_del_llenado
from .perfiles import ESTACIONES, FORMA_OR

MM, G = 1e-3, 1e-3
L_MASA = 1e-4  # m, largo de los Mass components (masas casi puntuales en su CG)
COMPONENTES = ("nariz", "cuerpo", "cola", "tubo_cola", "aletas")
TIPOS = {"nariz": "NoseCone", "cuerpo": "BodyTube", "cola": "Transition", "tubo_cola": "BodyTube",
         "aletas": "FreeformFinSet"}

VALIDACION = [
    "puesto", "cand_id", "validado", "motivo", "x_CP_modelo_mm", "x_CP_or_mm", "dx_cp_mm", "CNa_modelo", "CNa_or",
    "m_modelo_g", "m_or_g", "dif_masa_or_pct", "x_CG_modelo_mm", "x_CG_or_mm", "dx_cg_mm", "SM_modelo_cal",
    "SM_or_cal", "CD_or", "advertencias",
]


def densidad_pared(cav, capas, estacion: str) -> float:
    """Densidad de la capa única equivalente de OpenRocket para la pared multicapa de una estación:
    masa de las capas / volumen de la pared (ponderada por área). La ponderada por espesor
    (`densidad_equivalente`) solo vale para paredes planas: en un tubo de Ø 16 mm la capa exterior
    (fibra, más densa) tiene más área que la interior y la masa quedaba 0.7 % corta."""
    i = ESTACIONES.index(estacion)
    m = cav.estacion == i
    if m.sum() < 3:
        return densidad_equivalente(capas)
    fr = cav.fronteras[estacion]
    x = cav.x[m]
    vol = [float(simpson(np.pi * (fr[k - 1][m] ** 2 - fr[k][m] ** 2), x=x)) for k in range(1, len(capas) + 1)]
    V = sum(vol)
    return sum(c.rho * c.phi * v for c, v in zip(capas, vol)) / V if V > 0 else densidad_equivalente(capas)


@contextlib.contextmanager
def abrir_openrocket(cfg: ConfigOpt):
    """OpenRocket 24.12 vía orlab con el jar de `openrocket.jar` (o ORLAB_JAR / caché de orlab)."""
    import orlab
    with orlab.OpenRocketInstance(jar_path=cfg.openrocket.get("jar"), log_level="ERROR") as inst:
        yield PuenteTipo2(inst, cfg.resolver(cfg.openrocket.get("ork_base", "modelos/analisis_tipo2.ork")))


class ErrorAleta(RuntimeError):
    """OpenRocket no aceptó los puntos de la aleta tal como se pidieron."""


@dataclass
class ResultadoAero:
    x_CP: float
    CNa: float
    CD: float
    por_componente: dict[str, tuple[float, float]]  # (CNa, x_CP global)
    warnings: list[str]


@dataclass
class ResultadoMasas:
    m_total: float
    x_CG: float
    por_componente: dict[str, tuple[float, float]] = field(default_factory=dict)


class PuenteTipo2:
    """Envoltorio de una instancia de OpenRocket abierta con orlab:

        with orlab.OpenRocketInstance(jar_path=...) as inst:
            pr = PuenteTipo2(inst, "modelos/analisis_tipo2.ork")
    """

    def __init__(self, inst, ork_base: str | Path):
        import jpype
        import orlab

        self._jpype = jpype
        self.orl = orlab.Helper(inst)
        self.core = inst.openrocket
        self.ork_base = str(ork_base)
        self.doc = self.rocket = None
        self.comp: dict[str, object] = {}

    @property
    def version_or(self) -> str:
        return str(self.core.util.BuildProperties.getVersion())

    def cargar(self, ruta: str | Path | None = None):
        ruta = str(ruta or self.ork_base)
        self.doc = self.orl.load_doc(ruta)
        self.rocket = self.doc.getRocket()
        por_tipo: dict[str, list] = {}
        it = self.rocket.iterator(True)
        while it.hasNext():
            c = it.next()
            por_tipo.setdefault(str(c.getClass().getSimpleName()), []).append(c)
        try:
            tubos = por_tipo["BodyTube"]
            self.comp = {"nariz": por_tipo["NoseCone"][0], "cuerpo": tubos[0], "cola": por_tipo["Transition"][0],
                         "tubo_cola": tubos[1], "aletas": por_tipo["FreeformFinSet"][0]}
        except (KeyError, IndexError) as e:
            raise KeyError(f"{ruta}: se esperan nariz, tubo, transición, tubo y aletas freeform; hay "
                           f"{ {k: len(v) for k, v in por_tipo.items()} }") from e
        if self.comp["aletas"].getParent() != self.comp["tubo_cola"]:
            raise ValueError(f"{ruta}: las aletas deben colgar del tubo de cola")
        return self

    @property
    def configuracion(self):
        return self.rocket.getSelectedConfiguration()

    def x_abs(self, clave: str) -> float:
        return float(self.comp[clave].getComponentLocations()[0].x)

    def _material(self, nombre: str, rho: float):
        M = self.core.material.Material
        return M.newMaterial(M.Type.BULK, nombre, float(rho), True)

    def _forma(self, forma: str):
        return self.core.rocketcomponent.Transition.Shape.valueOf(FORMA_OR[forma])

    # ------------------------------------------------------------------ escritura
    def aplicar(self, cfg: ConfigOpt, det: Detalle, eps_puntos: float = 1e-7) -> dict:
        """Recarga la base y fija el cuerpo, las paredes y las aletas del candidato."""
        self.cargar()
        p, g = det.cu.perfil, det.g
        n, c, t, tc, f = (self.comp[k] for k in COMPONENTES)
        n.setShapeType(self._forma(p.forma_n))
        if p.param_n is not None and n.getShapeType().usesParameter():
            n.setShapeParameter(float(p.param_n))
        n.setLength(float(p.Ln))
        n.setBaseRadius(float(p.R))
        c.setLength(float(p.L_c))  # 0 en un cuerpo abombado: OpenRocket 24.12 lo acepta
        c.setOuterRadiusAutomatic(True)
        t.setShapeType(self._forma(p.forma_t))
        if p.param_t is not None and t.getShapeType().usesParameter():
            t.setShapeParameter(float(p.param_t))
        if t.getShapeType().isClippable():
            t.setClipped(bool(p.recortada))
        t.setLength(float(p.Lt))
        t.setForeRadiusAutomatic(True)
        t.setAftRadius(float(p.r_tc))
        tc.setLength(float(p.L_tc))
        tc.setOuterRadius(float(p.r_tc))
        for k in ("nariz", "cuerpo", "cola", "tubo_cola"):
            capas = cfg.pared[k]
            self.comp[k].setThickness(float(espesor_total(capas)))
            self.comp[k].setMaterial(self._material(f"PARED_EQ_{k}", densidad_pared(det.cu.cav, capas, k)))
        pa = g.params
        pts = [(float(x), float(y)) for x, y in g.puntos]
        C = self.core.util.Coordinate
        f.setPoints(self._jpype.JArray(C)([C(x, y) for x, y in pts]))
        f.setFinCount(int(pa.n))
        f.setBaseRotation(float(pa.rotacion))
        f.setThickness(float(pa.t))
        f.setMaterial(self._material(pa.material.upper(), pa.rho * pa.phi))
        f.setAxialMethod(self.core.rocketcomponent.position.AxialMethod.BOTTOM)
        f.setAxialOffset(0.0)
        aceptados = [(float(q.x), float(q.y)) for q in f.getFinPoints()]
        if len(aceptados) != len(pts) or np.max(np.abs(np.array(aceptados) - np.array(pts))) > eps_puntos:
            raise ErrorAleta(f"OpenRocket modificó los puntos de la aleta: {pts} → {aceptados}")
        return {"x_LE": self.x_abs("aletas"), "puntos": aceptados}

    def agregar_masas(self, masas: list[tuple[str, float, float]]):
        MC = self.core.rocketcomponent.MassComponent
        AM = self.core.rocketcomponent.position.AxialMethod
        padre = self.comp["cuerpo"]
        for nombre, m, x in masas:
            mc = MC(L_MASA, 0.001, float(m))
            mc.setName(nombre)
            padre.addChild(mc)
            mc.setAxialMethod(AM.ABSOLUTE)
            mc.setAxialOffset(float(x - L_MASA / 2))

    # ------------------------------------------------------------------ cálculos
    def aero(self, mach: float, aoa: float = 0.0) -> ResultadoAero:
        fc = self.configuracion
        calc = self.core.aerodynamics.BarrowmanCalculator()
        cond = self.core.aerodynamics.FlightConditions(fc)
        cond.setMach(float(mach))
        cond.setAOA(float(aoa))
        ws = self.core.logging.WarningSet()
        cp = calc.getCP(fc, cond, ws)
        CD = float(calc.getAerodynamicForces(fc, cond, self.core.logging.WarningSet()).getCD())
        ids = {str(v.getID()): k for k, v in self.comp.items()}
        por = {}
        for comp, fz in calc.getForceAnalysis(fc, cond, self.core.logging.WarningSet()).items():
            k = ids.get(str(comp.getID()))
            if k is not None:
                c = fz.getCP()
                por[k] = (float(c.weight), float(c.x) if float(c.weight) != 0 else math.nan)
        return ResultadoAero(x_CP=float(cp.x), CNa=float(cp.weight), CD=CD, por_componente=por,
                             warnings=[str(w) for w in ws])

    def masas(self) -> ResultadoMasas:
        rb = self.core.masscalc.MassCalculator.calculateStructure(self.configuracion)
        out = ResultadoMasas(m_total=float(rb.getMass()), x_CG=float(rb.getCM().x))
        for k in COMPONENTES:
            c = self.comp[k]
            out.por_componente[k] = (float(c.getComponentMass()), self.x_abs(k) + float(c.getComponentCG().x))
        return out

    def perfil(self, n: int) -> list[tuple[str, float, float]]:
        filas = []
        for k in ("nariz", "cuerpo", "cola", "tubo_cola"):
            c = self.comp[k]
            x0, L = self.x_abs(k), float(c.getLength())
            for xl in np.linspace(0.0, L, n):
                filas.append((k, x0 + xl, float(c.getRadius(float(xl)))))
        return filas


# --------------------------------------------------------------------------- validación


def validar_candidato(cfg: ConfigOpt, pr: PuenteTipo2, fila: pd.Series, tol_cp: float, tol_masa_pct: float) -> dict:
    """Arma el candidato en OpenRocket con el lastre del modelo propio y compara CP, C_Nα, masa y CG."""
    c, a = specs_de_fila(fila)
    det = detalle(cfg, c, a)
    out = {"puesto": fila.get("puesto"), "cand_id": fila["cand_id"], "validado": False}
    Ll = det.llenado
    if det.g is None or det.cp is None or Ll is None or not math.isfinite(Ll.ell):
        out["motivo"] = "modelo_sin_llenado:" + ";".join(det.motivos)
        return out
    try:
        pr.aplicar(cfg, det)
    except ErrorAleta as e:
        out.update({"motivo": "aleta_rechazada", "advertencias": str(e)})
        return out
    aero = pr.aero(cfg.vuelo.mach)
    pr.agregar_masas(masas_del_llenado(cfg, det))
    mas = pr.masas()
    f, D = Ll.fila, c.D
    out.update({
        "x_CP_modelo_mm": det.cp.x_CP / MM, "x_CP_or_mm": aero.x_CP / MM, "dx_cp_mm": (aero.x_CP - det.cp.x_CP) / MM,
        "CNa_modelo": det.cp.CNa, "CNa_or": aero.CNa, "m_modelo_g": f["m_total_g"], "m_or_g": mas.m_total / G,
        "dif_masa_or_pct": 100 * (mas.m_total / G / f["m_total_g"] - 1), "x_CG_modelo_mm": f["x_CG_mm"],
        "x_CG_or_mm": mas.x_CG / MM, "dx_cg_mm": mas.x_CG / MM - f["x_CG_mm"], "SM_modelo_cal": f["SM_cal"],
        "SM_or_cal": (aero.x_CP - mas.x_CG) / D, "CD_or": aero.CD, "advertencias": ";".join(aero.warnings)})
    for comp, (cna, x) in aero.por_componente.items():
        out[f"CNa_or_{comp}"], out[f"x_CP_or_{comp}_mm"] = cna, x / MM
    for comp, (m, _) in mas.por_componente.items():
        out[f"m_or_{comp}_g"] = m / G
    rest = cfg.restricciones
    fallas = []
    if abs(out["dx_cp_mm"]) > tol_cp / MM:
        fallas.append("dx_cp")
    if abs(out["dif_masa_or_pct"]) > tol_masa_pct:
        fallas.append("dif_masa")
    if not rest.SM_min - 1e-9 <= out["SM_or_cal"] <= rest.SM_max + 1e-9:
        fallas.append("SM_or_fuera_de_rango")
    out["validado"] = not fallas
    out["motivo"] = ";".join(fallas)
    return out


@dataclass
class ResultadoValidacion:
    tabla: pd.DataFrame
    ganador: str | None
    cp_fuera_de_tolerancia: list[str] = field(default_factory=list)


def validar_ganadores(cfg: ConfigOpt, pr: PuenteTipo2, ranking: pd.DataFrame, log=print) -> ResultadoValidacion:
    """Valida los N_ganadores mejores factibles. El ganador es el primero validado en el orden del
    ranking, salvo que algún candidato supere tol_cp_mm: entonces no hay ganador y se reporta."""
    vo = cfg.ejecucion.get("validacion_or") or {}
    N = int(vo.get("N_ganadores", 5))
    tol_cp = float(vo.get("tol_cp_mm", 2.0)) * MM
    tol_m = float(vo.get("tol_masa_pct", 0.5))
    filas = []
    for _, fila in ranking[ranking["factible"].astype(bool)].head(N).iterrows():
        r = validar_candidato(cfg, pr, fila, tol_cp, tol_m)
        filas.append(r)
        log(f"    {r['cand_id']}  x_CP {r.get('x_CP_modelo_mm', math.nan):7.2f} / {r.get('x_CP_or_mm', math.nan):7.2f} mm  "
            f"Δm {r.get('dif_masa_or_pct', math.nan):+.3f} %  SM_OR {r.get('SM_or_cal', math.nan):.3f}  "
            + ("validado" if r["validado"] else f"NO ({r.get('motivo', '')})"))
    tabla = pd.DataFrame(filas).reindex(columns=list(dict.fromkeys(VALIDACION + [k for f in filas for k in f])))
    fuera = [r["cand_id"] for r in filas if "dx_cp" in str(r.get("motivo", "")).split(";")]
    ok = tabla[tabla["validado"].fillna(False).astype(bool)] if len(tabla) else tabla
    ganador = None if fuera or ok.empty else str(ok.iloc[0]["cand_id"])
    return ResultadoValidacion(tabla, ganador, fuera)


# --------------------------------------------------------------------------- tolerancia de amarre de un .ork


@dataclass
class AmarreOrk:
    """Tolerancia de amarre de un .ork cualquiera, con el amarre en el CG (trim estático 0)."""

    m: float  # kg (estructura + Mass components)
    x_CG: float  # m desde la punta
    x_CP: float  # m desde la punta
    CNa: float
    S_ref: float  # m² (referencia de OpenRocket: la del mayor diámetro)
    D_ref: float  # m
    mach: float
    q: float  # Pa
    alpha_max: float  # rad
    tol: float  # m
    advertencias: list[str]
    L: float = math.nan  # m, largo total
    x_herraje: list[float] = field(default_factory=list)  # m, donde estaba cada herraje en el .ork
    m_herraje: float = 0.0  # kg, masa de los herrajes (se llevan al CG)
    x_CG_ork: float = math.nan  # m, CG tal como viene en el .ork (herraje donde estaba)

    @property
    def SM(self) -> float:
        return (self.x_CP - self.x_CG) / self.D_ref


def tolerancia_amarre_ork(cfg: ConfigOpt, ruta: str | Path, mach: float | None = None) -> AmarreOrk:
    """Abre un .ork con OpenRocket 24.12 y devuelve su tolerancia de amarre:

        tol = α_max q S_ref C_Nα (x_CP − x_CG) / (m g)

    con la masa y el CG de OpenRocket (estructura y Mass components), el CP y el C_Nα de Barrowman a
    α = 0 y el Mach, la presión dinámica y α_max de la configuración. No exige la topología del .ork
    tipo 2: sirve para cualquier cohete que OpenRocket abra.

    El herraje de remolque (Mass components cuyo nombre contiene "herraje") va siempre en el CG, como
    en el modelo: su masa cuenta, pero no su momento, esté donde esté en el .ork:

        x_CG = (m x_CG,ork − Σ m_h x_h) / (m − Σ m_h)"""
    from .lastre import tolerancia_amarre
    mach = cfg.vuelo.mach if mach is None else float(mach)
    with abrir_openrocket(cfg) as pr:
        rocket = pr.orl.load_doc(str(ruta)).getRocket()
        fc = rocket.getSelectedConfiguration()
        cond = pr.core.aerodynamics.FlightConditions(fc)
        cond.setMach(mach)
        cond.setAOA(0.0)
        ws = pr.core.logging.WarningSet()
        cp = pr.core.aerodynamics.BarrowmanCalculator().getCP(fc, cond, ws)
        rb = pr.core.masscalc.MassCalculator.calculateStructure(fc)
        m, x_CG = float(rb.getMass()), float(rb.getCM().x)
        x_CP, CNa = float(cp.x), float(cp.weight)
        S_ref, D_ref = float(cond.getRefArea()), float(cond.getRefLength())
        avisos = [str(w) for w in ws]
        L = float(fc.getLength())
        x_h, m_h, it = [], [], rocket.iterator(True)
        while it.hasNext():
            c = it.next()
            if "herraje" in str(c.getName()).lower():
                x_h.append(float(c.getComponentLocations()[0].x) + float(c.getComponentCG().x))
                m_h.append(float(c.getComponentMass()))
    x_CG_ork = x_CG
    if m_h:  # el herraje va en el CG: se quita su momento
        x_CG = (m * x_CG_ork - sum(mi * xi for mi, xi in zip(m_h, x_h))) / (m - sum(m_h))
    q, a_max = cfg.vuelo.q, cfg.remolque.alpha_max
    return AmarreOrk(m=m, x_CG=x_CG, x_CP=x_CP, CNa=CNa, S_ref=S_ref, D_ref=D_ref, mach=mach, q=q, alpha_max=a_max,
                     tol=tolerancia_amarre(m, x_CG, x_CP, q, S_ref, CNa, a_max), advertencias=avisos, L=L,
                     x_herraje=x_h, m_herraje=sum(m_h), x_CG_ork=x_CG_ork)
