import numpy as np
import pytest

from eqlat import (
    equivalent_latitude_piecewise,
    equivalent_latitude_roi,
    grid_cell_areas,
)
from eqlat.utils import total_sphere_area


@pytest.fixture
def zonal_field():
    """Zonally symmetric PV increasing with latitude: eqlat must equal latitude."""
    lat = np.linspace(-89.5, 89.5, 180)
    lon = np.arange(0.0, 360.0, 1.0)
    pv = np.repeat(np.sin(np.deg2rad(lat))[:, None], lon.size, axis=1)
    return pv, lat, lon


def test_grid_cell_areas_sum_to_sphere(zonal_field):
    _, lat, lon = zonal_field
    areas = grid_cell_areas(lat, lon)
    assert areas.sum() == pytest.approx(total_sphere_area(), rel=1e-3)


@pytest.mark.parametrize("method", [equivalent_latitude_piecewise, equivalent_latitude_roi])
def test_zonal_field_recovers_latitude(zonal_field, method):
    pv, lat, lon = zonal_field
    thresholds = np.sin(np.deg2rad(np.arange(-80.0, 81.0, 10.0)))
    res = method(pv, lat, lon, pv_thresholds=thresholds)
    expected = np.rad2deg(np.arcsin(thresholds))
    np.testing.assert_allclose(res["eqlat"], expected, atol=1.0)
