"""
structures.py

Helpers to create different starting environments.

"""

import quippy.structures

from sats.elements import Element


def bcc_interstitial_111(lattice_parameter, species, interstitial_species=None,
                         supercell=(1, 1, 1), index=-1):
    """
    Create a BCC crystal structure with a <111> dumbbell interstitial atom.

    Parameters
    ----------
    lattice_parameter : float
        BCC lattice parameter in A.
    species : int or str
        Identity of species used to construct the BCC lattice.
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
        BCC lattice with a 111 dumbbell interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    bulk = quippy.structures.bcc(a=lattice_parameter, z=species.atomic_number)
    # catch all false values for supercell
    if supercell:
        bulk = quippy.structures.supercell(bulk, supercell[0], supercell[1],
                                           supercell[2])

    if index < 0:
        index = bulk.n + index

    # bcc 111 dumbbell
    magnitude = (1.0/6)*lattice_parameter
    displacement = [magnitude, magnitude, magnitude]
    bulk.add_atoms(pos=bulk[index].position + displacement,
                   z=interstitial_species.z)
    bulk[index].position -= displacement

    return bulk
