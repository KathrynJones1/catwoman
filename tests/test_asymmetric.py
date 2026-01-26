import numpy as np
import catwoman
try:
    import twoface
except:
    print('For this test to run you need to install twoface:')
    print('pip install twoface')
import pytest

abs_tol = 1e-5


def test_0():
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]

    t = np.linspace(-0.7, 0.7, 10000)
    model = catwoman.TransitModel(params, t)

    t = np.linspace(-0.7 + 10*params.per, 0.7 + 10*params.per, 10000)
    model2 = catwoman.TransitModel(params, t)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_1():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = -90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_2():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = -45

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_3():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 0.

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_4():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 45

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_5():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_6():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 72.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_7():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 55.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_8():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.2
    params.w = 45.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_9():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.1
    params.w = 180.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_10():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 72.
    params.ecc = 0.1
    params.w = 180.
    params.rp = 0.1
    params.rp2 = 0.07
    params.phi = 90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)


def test_11():
    t = np.linspace(-7, 7, 10000)
    params = catwoman.TransitParams()
    params.per = 3.
    params.t0 = 0.
    params.a = 1.5
    params.inc = 90.
    params.ecc = 0.0
    params.w = 90.
    params.rp = 0.3
    params.rp2 = 0.0
    params.phi = 90

    params.limb_dark = 'quadratic'
    params.u = [0.3, 0.2]
    params.num_spots = 0

    model = catwoman.TransitModel(params, t)
    model2 = twoface.TransitModel(params, t, num_rings=3000)

    result = model.light_curve(params)
    expected = model2.light_curve(params)
    assert result == pytest.approx(expected, abs=abs_tol)
