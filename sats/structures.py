"""
structures.py

Helpers to create different starting environments.

"""

import re

import quippy.structures

from sats.elements import Element


def bcc_interstitial_dumbell(lattice_parameter, species, direction='111',
                             interstitial_species=None, supercell=(1, 1, 1),
                             index=-1):
    """
    Create a BCC crystal structure with a dumbbell interstitial atom along
    the given lattice direction.

    Parameters
    ----------
    lattice_parameter : float
        BCC lattice parameter in A.
    species : int or str
        Identity of species used to construct the BCC lattice.
    direction: str
        Lattice direction along which to create the dumbbell.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    supercell : tuple of int, optional
        Supercell size to place the defect in; if not specified, uses a unit
        cell.
    index : int
        Index of atom with which to create the dumbbell. Added atom is always
        at the end.

    Returns
    -------
    quippy.Atoms
        BCC lattice with a dumbbell interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    bulk = quippy.structures.bcc(a=lattice_parameter, z=species.atomic_number)
    # catch all False values for supercell
    if supercell:
        bulk = quippy.structures.supercell(bulk, supercell[0], supercell[1],
                                           supercell[2])

    if index < 0:
        index = bulk.n + index

    # parses the digits from any string, ignoring brackets etc.
    direction = [int(x) for x in re.findall('[+-]*[0-9]', direction)]
    displacement = [1.0/x if x else 0.0 for x in direction]
    # Move 1/3 of the <111> atom separation in each direction
    magnitude = lattice_parameter/(12*sum(x**2 for x in displacement))**0.5
    displacement = [x*magnitude for x in displacement]
    bulk.add_atoms(pos=bulk[index].position + displacement,
                   z=interstitial_species.z)
    bulk[index].position -= displacement

    return bulk

