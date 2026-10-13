import numpy as np

from ase2sprkkr.physics.broadening import create_gaussian_broadener


def _sorted_reference(x, y, sigma=None, fwhm=None, **kwargs):
    order = np.argsort(x)
    x_sorted = np.asarray(x)[order]
    y_sorted = np.asarray(y)[..., order]
    sigma_sorted = None if sigma is None or np.ndim(sigma) == 0 else np.asarray(sigma)[order]
    fwhm_sorted = None if fwhm is None or np.ndim(fwhm) == 0 else np.asarray(fwhm)[order]
    broadener = create_gaussian_broadener(
        x_sorted,
        sigma=sigma if sigma_sorted is None else sigma_sorted,
        fwhm=fwhm if fwhm_sorted is None else fwhm_sorted,
        **kwargs,
    )
    return broadener(y_sorted)


def test_gaussian_unsorted_grid_reorders_data_and_sigma():
    x = np.array([2.0, 0.0, 1.0, 3.0])
    y = np.array([20.0, 1.0, 7.0, 3.0])
    sigma = np.array([0.6, 0.15, 0.3, 0.8])

    actual = create_gaussian_broadener(x, sigma=sigma, method="dense")(y)
    expected_sorted = _sorted_reference(x, y, sigma=sigma, method="dense")

    np.testing.assert_allclose(actual[np.argsort(x)], expected_sorted)


def test_gaussian_unsorted_grid_reorders_fwhm():
    x = np.array([1.5, -0.5, 0.5, 2.5])
    y = np.array([2.0, 8.0, 1.0, 4.0])
    fwhm = np.array([0.9, 0.25, 0.5, 1.1])

    actual = create_gaussian_broadener(x, fwhm=fwhm, method="dense")(y)
    expected_sorted = _sorted_reference(x, y, fwhm=fwhm, method="dense")

    np.testing.assert_allclose(actual[np.argsort(x)], expected_sorted)


def test_gaussian_broadens_last_axis_and_supports_distinct_output_grid():
    x = np.array([2.0, 0.0, 1.0])
    x_out = np.array([1.5, -0.5, 0.5, 2.5])
    y = np.array([[20.0, 1.0, 7.0], [2.0, 4.0, 8.0]])

    actual = create_gaussian_broadener(
        x,
        sigma=0.4,
        x_out=x_out,
        method="dense",
    )(y)

    order = np.argsort(x)
    reference = create_gaussian_broadener(
        x[order],
        sigma=0.4,
        x_out=x_out,
        method="dense",
    )(y[..., order])

    assert actual.shape == (2, len(x_out))
    np.testing.assert_allclose(actual, reference)


def test_gaussian_rejects_wrong_energy_axis_length():
    broaden = create_gaussian_broadener([0.0, 1.0, 2.0], sigma=0.2)
    with np.testing.assert_raises(ValueError):
        broaden(np.ones((2, 4)))
