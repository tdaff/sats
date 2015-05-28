"""
vacancy.py

Make crystals with vacancies or vacancy rows.

"""

from numpy import cross
from numpy.linalg import norm

from sats.core.elements import Element
from sats.configurations.bulk import bulk


def vacancy(lattice, species, n_vacancies=1, min_atoms=0, direction=None):
    """
    Helper function to create a vacancy in a lattice. Multiple vacancies are
    created in a line along close packed directions.

    Parameters
    ----------
    lattice : str
        Lattice type for the bulk structure. Generated as a sats bulk.
    species : int or str
        Identity of species used to construct the lattice.
    n_vacancies : int, optional
        Number of vacancies to create in a row.
    min_atoms : int, optional
        Construct a supercell to include at least this many atoms.
    direction : tuple of int, optional
        Create vacancy row along this direction. If none is specified then
        the close packed direction is used.

    Returns
    -------
    vacancy : ase.Atoms
        Lattice with the requested defect.
    """

    species = Element(species)

    vbulk = bulk(species, lattice, min_atoms)

    if direction is None:
        direction = {
            'bcc': (1, 1, 1),
            'hcp': (1/3.0**0.5, 1, 0),
            'fcc': (1, 1, 0)
        }[lattice]

    # order from lowest x, y and z to get an origin
    sites = [atom.index for atom in
             sorted(vbulk, key=lambda a: (a.c, a.a, a.b))]

    origin = vbulk[sites[0]]

    line_atoms = []
    # index all the atoms along the line from origin
    for atom in vbulk:
        separation = atom.position - origin.position
        if norm(cross(separation, direction)) < 1.0e-5:
            line_atoms.append((norm(separation), atom.index))

    for atom in sorted(line_atoms, reverse=True)[-n_vacancies:]:
        del vbulk[atom[1]]

    return vbulk
