from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

RAIZ = Path(__file__).resolve().parents[1]
CONFIG = RAIZ / "config" / "optimizacion.yaml"
ORK = RAIZ / "modelos" / "analisis_tipo2.ork"


def raw_opt() -> dict:
    with open(CONFIG, encoding="utf-8") as f:
        return copy.deepcopy(yaml.safe_load(f))


def raw_ork() -> dict:
    """Config con la pared y las aletas tal como están en el .ork (capa única de 2 mm, ρ = 1342)."""
    raw = raw_opt()
    raw["materiales"]["MAT_PARED_EQ"] = 1342
    raw["pared"]["por_defecto"] = [{"material": "MAT_PARED_EQ", "espesor_mm": 2.0}]
    raw["geometria_fija"]["aletas"].update(material="MAT_PARED_EQ", rotacion_deg=0)
    raw["condiciones_vuelo"]["mach"] = 0.3
    raw["restricciones"]["L_n_rel_D_min"] = 0.0  # la nariz del .ork es de 0.93 D
    return raw


def raw_libre() -> dict:
    """Config por defecto sin los mínimos geométricos de la spec v2 (nariz, cuerda, tubo, tope de
    r_tip): para probar la física sobre geometrías arbitrarias, como la del .ork."""
    raw = raw_opt()
    raw["restricciones"].update(L_n_rel_D_min=0.0, c_r_min_mm=0.0, k_min=0.0, d_tc_min_mm=0.0, r_tip_rel_R_max=None)
    return raw


def malla_chica(raw: dict, **cambios) -> dict:
    """Malla reducida (con factibles bajo la spec v2) para tests rápidos."""
    raw["malla"] = {
        "L_total_mm": [400], "D_mm": [80, 90], "L_n_rel_D": [1.0], "f_cil": [0.0, 0.5],
        "cola_forma": [{"forma": "elipsoide", "parametro": None}, {"forma": "conica", "parametro": None}],
        "k": [0.175, 0.3], "L_tc_mm": [40, 125, 140],
        "mu_cr": [0.7], "gamma_ct": [0.7], "sigma_flecha": [1.0], "r_tip_rel_R": [1.0, 1.1, 1.2],
    }
    raw["malla"].update(cambios)
    raw["ejecucion"]["refinamiento"] = {"activar": True, "top_K": 2}
    return raw


@pytest.fixture
def raw() -> dict:
    return raw_opt()


@pytest.fixture(scope="session")
def cfg():
    from sensor_tipo2.config import cargar
    return cargar(CONFIG)


@pytest.fixture(scope="session")
def cfg_libre():
    from sensor_tipo2.config import cargar
    return cargar(raw_libre())


@pytest.fixture(scope="session")
def cfg_ork():
    from sensor_tipo2.config import cargar
    return cargar(raw_ork())


# mejor de la malla por defecto de la spec v2: cuerpo abombado, aletas en 1.2 R, electrónica de 20 mm
FACTIBLE = ((400, 90, 1.0, "conica", None, 0.0, 15.75, 140), (0.7, 0.7, 1.0, 1.2))


@pytest.fixture(scope="session")
def base_ork():
    """Cuerpo y aleta del .ork tipo 2."""
    from sensor_tipo2.config import AletaSpec, CuerpoSpec
    return CuerpoSpec(400, 70, 65 / 70, "elipsoide", None, 220 / 285, 30, 50), AletaSpec(1.0, 0.7, 1.0, 1.0)
