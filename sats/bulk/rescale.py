"""
rescale.py

Estimate the bulk linear expansion at different state points.
"""

from __future__ import division

from sats.core.elements import Element
from sats.ui.log import warning, error

properties = {
    'W': {
        'lattice_constant': {
            'bcc': 3.18,
            'fcc': 3.95,
            'hcp': 2.79,
            'a15': 5.06
        },
        'thermal_expansion': 4.6e-6,
        'thermal_expansion_ref': "10.6028/nbsscipaper.203",
        'melting_point': 3695,
        'melting_density': "17.6 g.cm-3",
        'bulk_modulus': 310.691,
        'bulk_modulus_ref': '10.1016/S0022-3697(03)00155-0',
        'bulk_modulus_derivative': 3.94,  # approx for alpha
        'bulk_modulus_derivative_ref': '10.1016/S0022-3697(03)00155-0',
    },
    'Fe': {
        'lattice_constant': {
            'bcc': 2.87,
            'fcc': 3.45,
            'hcp': 2.46,
            'a15': 4.48
        },
        'thermal_expansion': 14.5e-6,
        'thermal_expansion_ref': "10.1098/rspa.1955.0102",
        'melting_point': 1808,
        'melting_density': 6.98,
        'bulk_modulus': 170.0,  # approx
        'bulk_modulus_ref': '10.1029/JB091iB05p04677',
        'bulk_modulus_derivative': 5.0,  # approx for alpha
        'bulk_modulus_derivative_ref': '10.1029/JB091iB05p04677',

    }
}


def bulk_rescale_factor(material=None, temperature=0, pressure=0,
                        volume_percent=0):
    """
    Estimated linear bulk scaling factor at the given state point. Uses
    linear thermal expansion with some checks for melting and a Birch
    Murnaghan equation of state for compression.

    Parameters
    ----------
    material : int or str or None
        Material for which to estimate change in lattice. If None then thermal
        and mechanical changes will be zero.
    temperature : float
        Temperature of state point to calculate expansion for.
    pressure : float
        Pressure of state point to calculate change for.
    volume_percent : float
        Percent change in the volume of the material.


    Returns
    -------
    linear_scale_factor : float
        Linear scale factor for expansion or contraction under the given
        conditions.
    """

    try:
        species = Element(material)
    except ValueError:
        species = material

    thermal_factor = 1.0
    pressure_factor = 1.0
    percent_factor = (1.0 + volume_percent/100)**(1/3)

    if temperature > 0:
        if material is None or material not in properties:
            warning("Unknown material {0}; no thermal contributions "
                    "included".format(material))
        else:
            thermal_factor = thermal_expansion(
                temperature, properties[species]['thermal_expansion'])

    if pressure != 0:
        if material is None or material not in properties:
            warning("Unknown material {0}; no pressure contributions "
                    "included".format(material))
        else:
            pressure_factor = birch_murnaghan_volume(
                pressure, properties[species]['bulk_modulus'],
                properties[species]['bulk_modulus_derivative'])**(1/3)

    # Can multiply as unspecified values will be 1.0 anyway
    return thermal_factor*pressure_factor*percent_factor


def thermal_expansion(temperature, coefficient):
    """
    Parameters
    ----------
    coefficient : float
        Coefficient of linear thermal expansion, in K^{-1}
    temperature : float
        Temperature of the material in Kelvin.

    Returns
    -------
    expansion : float
        Linear expansion as l/l0 fraction
    """

    return coefficient*temperature + 1


def birch_murnaghan(gamma, bulk_modulus, bulk_modulus_derivative):
    """
    Birch-Murnaghan equation of state in terms of gamma = V/V0.

    Parameters
    ----------
    gamma : float
        Ratio of volume to zero pressure volume, V/V0
    bulk_modulus : float
        Bulk modulus of the material, B0; units should be the same as the
        derivative.
    bulk_modulus_derivative : float
        Derivative of the bulk modulus with respect to pressure, B0'

    Returns
    -------
    pressure : float
        pressure in the same units as the input parameters.

    """
    # equation is in terms of V0/V
    gamma = 1/gamma

    pressure = (
        (3*bulk_modulus/2) *
        (gamma**(7/3) - gamma**(5/3)) *
        (1 + (3/4)*(bulk_modulus_derivative - 4) * (gamma**(2/3) - 1)))

    return pressure


def birch_murnaghan_volume(pressure, bulk_modulus, bulk_modulus_derivative):
    """
    Birch-Murnaghan equation of state in terms of the volume.

    Parameters
    ----------
    pressure : float
        pressure in the same units as the bulk parameters.
    bulk_modulus : float
        Bulk modulus of the material, B0; units should be the same as the
        derivative.
    bulk_modulus_derivative : float
        Derivative of the bulk modulus with respect to pressure, B0'.

    Returns
    -------
    gamma : float
        Ratio of volume to zero pressure volume, V/V0. No units.
    """

    # using a secant solver
    max_iter = 100
    tolerance = 1.0e-6

    gamma = 1.0
    trial_gamma = 1.1

    new_error = birch_murnaghan(trial_gamma, bulk_modulus,
                                bulk_modulus_derivative) - pressure

    for _i in range(max_iter):
        if abs(gamma - trial_gamma) < tolerance:
            # Converged
            return gamma

        past_gamma, trial_gamma, past_error = trial_gamma, gamma, new_error
        new_error = birch_murnaghan(trial_gamma, bulk_modulus,
                                    bulk_modulus_derivative) - pressure
        gamma = (trial_gamma -
                 new_error*(trial_gamma-past_gamma)/(new_error-past_error))
    else:
        # Did not converge!
        error("Max iterations reached finding Birch Murnaghan volume.")
        return gamma
