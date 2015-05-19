"""
surfaces.py

Surfaces and atom-by-atom layers to create all steps.

"""

import ase.lattice.surface

from sats.core.elements import Element


def intermediate_surfaces(species, lattice_parameter, surface='bcc100',
                          depth=12, vacuum=6):
    """
    Generate five 2x2 surfaces with top and bottom atoms progressively removed.
    Produces structures with vacancies, steps and adatoms.

    Parameters
    ----------
    species : str or int
        Species used to create surface.
    lattice_parameter : float
        Crystal lattice parameter, in A.
    surface : str
        Surface to create. Valid surfaces are bcc100, bcc110, bcc111, fcc100,
        fcc110, fcc111, fcc211, hcp0001.
    depth : int, optional
        Number of layers in the surface with all adatoms removed. Flat surfaces
        will have depth+2 layers total.
    vacuum : float, optional
        Height of vacuum to be included, in A.

    Returns
    -------
    surfaces : list of ase.atoms.Atoms
        Five surfaces with the arrangement of occupied (X) and vacant (O)
        sites on top, with an inversion for bottom sites:
            XX  XX  OX  OX  OX  b
            XX  OX  OX  XO  OO  ^>a
    """

    species = Element(species)
    surface = surface.lower()

    standard_surfaces = {
        'bcc100': ase.lattice.surface.bcc100,
        'bcc110': ase.lattice.surface.bcc110,
        'bcc111': ase.lattice.surface.bcc111,
        'fcc100': ase.lattice.surface.fcc100,
        'fcc110': ase.lattice.surface.fcc110,
        'fcc111': ase.lattice.surface.fcc111,
        'fcc211': ase.lattice.surface.fcc211,
        'hcp0001': ase.lattice.surface.hcp0001
    }

    if surface in standard_surfaces:
        flat_surface = standard_surfaces[surface](
            species, size=(2, 2, depth), a=lattice_parameter, vacuum=vacuum)
    else:
        raise NotImplementedError('Cannot create {} surfaces'.format(surface))

    # sort indexes along c, then a and b; this assumes that the supercell
    # has all c identical for images (seems to be fine)
    # Indexes are:
    #     Top:  23    Bottom:  67
    #           01             45

    sites = [atom.index for atom in
             sorted(flat_surface, key=lambda a: (a.c, a.a, a.b))]
    sites = sites[:4] + sites[-4:]

    # Surface 1:
    # XX
    # XX
    new_surfaces = [flat_surface]

    # Surface 2:
    # XX
    # OX
    tmp_surface = flat_surface.copy()
    del tmp_surface[[sites[0], sites[7]]]
    new_surfaces.append(tmp_surface)

    # Surface 3:
    # OX
    # OX
    tmp_surface = flat_surface.copy()
    del tmp_surface[[sites[0], sites[2], sites[5], sites[7]]]
    new_surfaces.append(tmp_surface)

    # Surface 4:
    # OX
    # XO
    tmp_surface = flat_surface.copy()
    del tmp_surface[[sites[1], sites[2], sites[4], sites[7]]]
    new_surfaces.append(tmp_surface)

    # Surface 5:
    # OX
    # O0
    tmp_surface = flat_surface.copy()
    del tmp_surface[[sites[0], sites[1], sites[2], sites[5], sites[6], sites[7]]]
    new_surfaces.append(tmp_surface)

    return new_surfaces
