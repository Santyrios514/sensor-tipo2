"""Verificación y calibración con OpenRocket 24.12 (orlab + JPype). Único módulo con Java.

`PuenteTipo2` abre `modelos/analisis_tipo2.ork`, ubica sus cinco componentes por tipo y orden
(los dos tubos se llaman igual, "Body Tube"), fija la geometría del candidato (nariz, cuerpo,
transición, tubo de cola, paredes con ρ_eq y aletas freeform sobre el tubo con offset *bottom* 0)
y agrega el lastre, la electrónica y las masas puntuales como *Mass components*.

Flujo (`verificar_y_calibrar`, el de dbf-sensor):
1. se verifican los N_verif mejores y una muestra de N_cal estratificada por (D, k, f_c, L_tc);
2. se ajusta x_OR ≈ α + β x_sust; si el residuo máximo supera tol_cp_mm se recalcula toda la
   malla con el CP calibrado y se verifican los nuevos mejores (hasta max_iter_calibracion);
3. cada verificado se vuelve a llenar con el CP y el C_Nα de OpenRocket (SM_or, J_or); gana el
   mejor J_or entre los factibles con OpenRocket.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .barrido import evaluar, grupos_desde_ranking, specs_de_fila
from .config import ConfigOpt, densidad_equivalente, espesor_total
from .exportar import Detalle, detalle, masas_del_llenado
from .objetivo import J, f_valores
from .perfiles import FORMA_OR
from .sustituto import IDENTIDAD, Calibracion, ajustar

MM, G = 1e-3, 1e-3
L_MASA = 1e-4  # m, largo de los Mass components (masas casi puntuales en su CG)
COMPONENTES = ("nariz", "cuerpo", "cola", "tubo_cola", "aletas")
TIPOS = {"nariz": "NoseCone", "cuerpo": "BodyTube", "cola": "Transition", "tubo_cola": "BodyTube",
         "aletas": "FreeformFinSet"}

VERIFICACION = [
    "cand_id", "rol", "iteracion", "factible_or", "motivos_or", "x_CP_sust_mm", "x_CP_cal_mm", "x_CP_or_mm",
    "dx_or_sust_mm", "dx_or_cal_mm", "CNa_or", "CNa_sust", "CD_or", "m_or_g", "x_CG_or_mm", "m_modelo_g",
    "x_CG_modelo_mm", "dif_masa_or_pct", "SM_or_cal", "J_or", "restriccion_activa_or", "advertencias",
]


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
            self.comp[k].setMaterial(self._material(f"PARED_EQ_{k}", densidad_equivalente(capas)))
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

    def guardar(self, ruta: Path):
        ruta.parent.mkdir(parents=True, exist_ok=True)
        self.orl.save_doc(str(ruta), self.doc)

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


# --------------------------------------------------------------------------- verificación de un candidato


def verificar_candidato(cfg: ConfigOpt, pr: PuenteTipo2, fila: pd.Series, cal: Calibracion = IDENTIDAD,
                        guardar: Path | None = None) -> dict:
    """CP, C_Nα, C_D, masas y CG de OpenRocket; llenado, SM y J con el CP de OpenRocket."""
    c, a = specs_de_fila(fila)
    det = detalle(cfg, c, a, cal)
    out = {"cand_id": fila["cand_id"], "x_CP_sust_mm": fila.get("x_CP_sust_mm"),
           "x_CP_cal_mm": float(cal.aplicar(fila["x_CP_sust_mm"] * MM)) / MM, "CNa_sust": fila.get("CN_alpha_total")}
    if det.g is None or det.cp is None:
        out.update({"factible_or": False, "motivos_or": ";".join(det.motivos)})
        return out
    pr.aplicar(cfg, det)
    aero = pr.aero(cfg.vuelo.mach)
    out.update({"x_CP_or_mm": aero.x_CP / MM, "dx_or_sust_mm": aero.x_CP / MM - out["x_CP_sust_mm"],
                "dx_or_cal_mm": aero.x_CP / MM - out["x_CP_cal_mm"], "CNa_or": aero.CNa, "CD_or": aero.CD,
                "advertencias": ";".join(aero.warnings)})
    for comp, (cna, x) in aero.por_componente.items():
        out[f"CNa_or_{comp}"], out[f"x_CP_or_{comp}_mm"] = cna, x / MM
    det_or = detalle(cfg, c, a, x_CP=aero.x_CP, CNa=aero.CNa)
    Ll = det_or.llenado
    if Ll is None or not math.isfinite(Ll.ell):
        out.update({"factible_or": False, "SM_or_cal": math.nan,
                    "motivos_or": ";".join(det_or.motivos) or "sin_llenado"})
        return out
    f = Ll.fila
    pr.agregar_masas(masas_del_llenado(cfg, det_or))
    mas = pr.masas()
    D = c.D
    out.update({"m_or_g": mas.m_total / G, "x_CG_or_mm": mas.x_CG / MM, "m_modelo_g": f["m_total_g"],
                "x_CG_modelo_mm": f["x_CG_mm"], "dif_masa_or_pct": 100 * (mas.m_total / G / f["m_total_g"] - 1),
                "SM_or_cal": (aero.x_CP - f["x_CG_mm"] * MM) / D, "m_total_or_g": f["m_total_g"],
                "restriccion_activa_or": Ll.restriccion_activa, "factible_or": Ll.ok,
                "motivos_or": ";".join(Ll.motivos), "tol_amarre_or_mm": f["tol_amarre_mm"]})
    for comp, (m, _) in mas.por_componente.items():
        out[f"m_or_{comp}_g"] = m / G
    F = f_valores(cfg, [f["m_total_g"] * G], [fila[cfg.objetivo.col_diametro] * MM], [fila["k_efectivo"]],
                  [out["SM_or_cal"]])
    out["J_or"] = float(J(cfg, F)[0])
    if guardar is not None:
        pr.guardar(guardar)
    return out


# --------------------------------------------------------------------------- muestra y calibración


def muestra_estratificada(df: pd.DataFrame, n: int, excluir: set[str], semilla: int) -> list[str]:
    """n candidatos factibles repartidos por estratos (D, k, f_c, L_tc), en ronda, con semilla fija,
    para cubrir cuerpos abombados y tubos largos (spec v2 §7)."""
    rng = np.random.default_rng(semilla)
    fac = df[df["factible"].astype(bool) & ~df["cand_id"].isin(excluir)]
    estratos = [list(g["cand_id"]) for _, g in fac.groupby(["D_mm", "k", "f_cil", "L_tc_mm"], sort=True)]
    for e in estratos:
        rng.shuffle(e)
    out: list[str] = []
    while len(out) < n and any(estratos):
        for e in estratos:
            if e and len(out) < n:
                out.append(e.pop())
    return out


@dataclass
class ResultadoVerificacion:
    ranking: pd.DataFrame
    verificacion: pd.DataFrame
    calibracion: Calibracion
    historial: list[dict] = field(default_factory=list)
    ganador: str | None = None


def verificar_y_calibrar(cfg: ConfigOpt, pr: PuenteTipo2, ranking: pd.DataFrame, procesos: int | None = None,
                         log=print) -> ResultadoVerificacion:
    vo = cfg.ejecucion.get("verificacion_or") or {}
    N_verif, N_cal = int(vo.get("N_verif", 20)), int(vo.get("N_cal", 30))
    tol = float(vo.get("tol_cp_mm", 2.0)) * MM
    max_iter = int(vo.get("max_iter_calibracion", 3))
    semilla = int(cfg.ejecucion.get("semilla", 12345))
    cal = IDENTIDAD
    df = ranking
    hechos: dict[str, dict] = {}
    historial = []

    def verificar(ids: list[str], rol: str, it: int):
        filas = df.set_index("cand_id")
        for cid in ids:
            if cid in hechos:
                if rol not in hechos[cid]["rol"]:
                    hechos[cid]["rol"] += f"+{rol}"
                continue
            fila = filas.loc[cid].copy()
            fila["cand_id"] = cid
            try:
                r = verificar_candidato(cfg, pr, fila)
            except ErrorAleta as e:
                r = {"cand_id": cid, "factible_or": False, "advertencias": str(e)}
            r.update({"rol": rol, "iteracion": it})
            hechos[cid] = r
            log(f"    {rol:13s} {cid}  x_sust {r.get('x_CP_sust_mm', math.nan):7.2f}  "
                f"x_OR {r.get('x_CP_or_mm', math.nan):7.2f} mm  SM_OR {r.get('SM_or_cal', math.nan):.2f}")

    top = list(df[df["factible"].astype(bool)].head(N_verif)["cand_id"])
    log(f"  Iteración 0: {len(top)} mejores + muestra de calibración de {N_cal}")
    verificar(top, "verificacion", 0)
    verificar(muestra_estratificada(df, N_cal, set(top), semilla), "calibracion", 0)
    for it in range(1, max_iter + 1):
        v = pd.DataFrame(hechos.values())
        if "x_CP_or_mm" not in v:
            break
        v = v[np.isfinite(v["x_CP_or_mm"].astype(float))]
        if len(v) < 2:
            break
        cal = ajustar(v["x_CP_sust_mm"] * MM, v["x_CP_or_mm"] * MM, iteracion=it)
        historial.append(cal.a_dict())
        log(f"  Calibración {it}: α = {cal.alpha / MM:+.3f} mm, β = {cal.beta:.5f}, R² = {cal.R2:.5f}, "
            f"residuo máx. {cal.res_max / MM:.3f} mm")
        if cal.res_max <= tol:
            log(f"  Residuo ≤ {tol / MM:g} mm: el sustituto basta, no se recalcula la malla.")
            break
        log(f"  Recalculando {len(df):,} candidatos con el CP calibrado…")
        nuevo = evaluar(cfg, grupos_desde_ranking(df), cal, procesos)
        if "refinamiento" in df:
            nuevo = nuevo.merge(df[["cand_id", "refinamiento"]], on="cand_id", how="left")
        df = nuevo
        nuevo_top = list(df[df["factible"].astype(bool)].head(N_verif)["cand_id"])
        if set(nuevo_top) <= set(hechos):
            log("  Los mejores no cambiaron: fin de la calibración.")
            break
        verificar([c for c in nuevo_top if c not in hechos], "verificacion", it)

    ver = pd.DataFrame(hechos.values()) if hechos else pd.DataFrame(columns=VERIFICACION)
    if len(ver) and "x_CP_sust_mm" in ver:
        ver["x_CP_cal_mm"] = cal.aplicar(ver["x_CP_sust_mm"].astype(float) * MM) / MM
        if "x_CP_or_mm" in ver:
            ver["dx_or_cal_mm"] = ver["x_CP_or_mm"] - ver["x_CP_cal_mm"]
    df = _anotar(df, ver)
    ok = ver[ver["factible_or"].fillna(False).astype(bool)] if "factible_or" in ver else ver.iloc[0:0]
    ganador = str(ok.sort_values(["J_or", "cand_id"]).iloc[0]["cand_id"]) if not ok.empty else None
    return ResultadoVerificacion(df, ver, cal, historial, ganador)


def _anotar(df: pd.DataFrame, ver: pd.DataFrame) -> pd.DataFrame:
    cols = ["x_CP_or_mm", "SM_or_cal", "dif_masa_or_pct", "J_or", "m_total_or_g", "x_CG_or_mm", "CD_or"]
    df = df.drop(columns=[c for c in cols if c in df])
    if len(ver):
        df = df.merge(ver[["cand_id"] + [c for c in cols if c in ver]], on="cand_id", how="left")
    df["verificado_or"] = df["cand_id"].isin(ver["cand_id"]) if len(ver) else False
    return df
