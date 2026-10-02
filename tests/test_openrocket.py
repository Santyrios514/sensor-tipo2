"""Integración con OpenRocket 24.12 (orlab + JPype). Se salta si no hay java, orlab o jar."""

import math
import shutil

import pandas as pd
import pytest

from .conftest import ORK

orlab = pytest.importorskip("orlab")
pytest.importorskip("jpype")
if shutil.which("java") is None:
    pytest.skip("sin java", allow_module_level=True)

MM = 1e-3


@pytest.fixture(scope="module")
def puente():
    from sensor_tipo2.verificacion import PuenteTipo2
    try:
        inst = orlab.OpenRocketInstance(log_level="ERROR")
        inst.__enter__()
    except Exception as e:  # sin jar en caché
        pytest.skip(f"OpenRocket no disponible: {e}")
    yield PuenteTipo2(inst, ORK)
    inst.__exit__(None, None, None)


def test_base_ork_regresion(puente, cfg_ork, base_ork):
    """El .ork tipo 2 reescrito por el puente da lo mismo que el archivo original (M = 0.3):
    CP = −133.646 mm, C_Nα = 1.25208, masa 215.611 g; y el sustituto coincide."""
    from sensor_tipo2.exportar import detalle
    det = detalle(cfg_ork, *base_ork)
    puente.aplicar(cfg_ork, det)
    ae = puente.aero(0.3)
    assert ae.x_CP / MM == pytest.approx(-133.646, abs=0.01)
    assert ae.CNa == pytest.approx(1.252078, rel=1e-5)
    assert puente.masas().m_total * 1e3 == pytest.approx(215.611, abs=0.01)
    assert det.cp.x_CP / MM == pytest.approx(ae.x_CP / MM, abs=0.05)
    m_mod = sum(m for m, _ in det.cu.masas_por_estacion().values()) + det.g.masa
    assert m_mod == pytest.approx(puente.masas().m_total, rel=0.005)


@pytest.mark.parametrize("cuerpo, aleta", [  # frontera masa–diámetro (D acostado = D, convención del .ork)
    ((400, 115, 0.2, "conica", None, 0.8, 23.75, 95), (0.85, 1.0, 0.5, 1.41421356)),
    ((400, 100, 0.2, "conica", None, 0.9, 20, 100), (0.85, 0.85, 1.0, 1.41421356)),
    ((400, 90, 0.2, "conica", None, 1.0, 20, 100), (0.85, 0.85, 1.0, 1.41421356)),
])
def test_candidatos_contra_openrocket(puente, cfg, cuerpo, aleta):
    from sensor_tipo2.barrido import evaluar_cuerpo
    from sensor_tipo2.config import AletaSpec, CuerpoSpec
    from sensor_tipo2.verificacion import verificar_candidato
    fila = pd.Series(evaluar_cuerpo(cfg, CuerpoSpec(*cuerpo), [AletaSpec(*aleta)])[0])
    assert fila["factible"]
    r = verificar_candidato(cfg, puente, fila)
    assert abs(r["dx_or_sust_mm"]) < 0.5
    assert abs(r["dif_masa_or_pct"]) < 0.5
    assert abs(r["x_CG_or_mm"] - r["x_CG_modelo_mm"]) < 0.5
    assert r["factible_or"] and math.isfinite(r["J_or"])
