import numpy as np
import catwoman
import batman

import pytest

abs_tol = 1e-5


def test_quadratic_ld():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.1
    params.phi = -45

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]

    model = catwoman.TransitModel(params, t)
    model2 = batman.TransitModel(params, t)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_power2_ld():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.1
    params.phi = -45

    params.limb_dark = 'power2'
    params.u = [0.3, 0.2]

    model = catwoman.TransitModel(params, t)
    model2 = batman.TransitModel(params, t)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_logarithmic_ld():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.1
    params.phi = -45

    params.limb_dark = 'logarithmic'
    params.u = [0.3, 0.2]

    model = catwoman.TransitModel(params, t)
    model2 = batman.TransitModel(params, t)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_exponential_ld():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.1
    params.phi = -45

    params.limb_dark = 'exponential'
    params.u = [0.3, 0.2]

    model = catwoman.TransitModel(params, t)
    model2 = batman.TransitModel(params, t)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_nonlinear_ld():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.1
    params.phi = -45

    params.limb_dark = 'nonlinear'
    params.u = [0.2, 0.1, 0.01, 0.01]

    model = catwoman.TransitModel(params, t)
    model2 = batman.TransitModel(params, t)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_custom_ld():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.1
    params.phi = -45

    params.limb_dark = 'custom'
    params.u = [0.2, 0.1, 0.01, 0.01, 0.0, 0.0]

    model = catwoman.TransitModel(params, t)
    model2 = batman.TransitModel(params, t)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)