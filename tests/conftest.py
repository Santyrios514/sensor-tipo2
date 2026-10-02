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
    return raw


def malla_chica(raw: dict, **cambios) -> dict:
    """Malla reducida para tests rápidos."""
    raw["malla"] = {
        "L_total_mm": [400], "D_mm": [70, 90], "L_n_rel_D": [0.4],
        "cola_forma": [{"forma": "elipsoide", "parametro": None}, {"forma": "conica", "parametro": None}],
        "L_t_rel_D": [1.25, 1.5], "d_tc_mm": [20, 25], "L_tc_mm": [120, 150],
        "mu_cr": [0.55], "gamma_ct": [0.85], "sigma_flecha": [1.0], "r_tip_rel_R": [1.0, 1.35, 1.41421356, 1.6],
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
def cfg_ork():
    from sensor_tipo2.config import cargar
    return cargar(raw_ork())


# candidato factible con aletas dentro del D acostado (r_tip/R = √2) y la electrónica de 20 mm
FACTIBLE = ((400, 90, 0.4, "conica", None, 1.5, 20, 150), (0.55, 0.85, 1.0, 1.41421356))


@pytest.fixture(scope="session")
def base_ork():
    """Cuerpo y aleta del .ork tipo 2."""
    from sensor_tipo2.config import AletaSpec, CuerpoSpec
    return CuerpoSpec(400, 70, 65 / 70, "elipsoide", None, 65 / 70, 30, 50), AletaSpec(1.0, 0.7, 1.0, 1.0)
