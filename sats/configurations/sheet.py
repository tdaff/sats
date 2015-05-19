"""
sheet.py


Sheet like structures for 2d simulation with no periodicity perpendicular
to the plane.

"""

from math import sin, cos, radians, ceil

from ase.atoms import Atoms

from sats.core.elements import Element
from sats.ui.options import add_option

# Not using as a decorator, just add all the sheet types here
add_option('structure', 'types', [], str, 'Type of structure to work with.',
           multiple=True, valid_values=['db', 'crw', 'tet', 'oct'])

def sheet(species, alat, blat=None, gamma=90.0, supercell=(1, 1),
          orthorhombic=False, depth=10.0):
    """
    Create a sheet of atoms in a 2d periodic lattice.

    Parameters
    ----------
    species : int or str
        Identity of species used to construct the sheet.
    alat : float
        Length of the a lattice parameter.
    blat : float, optional
        Length of the b lattice parameter, if this is not given it defaults
        to the same value as a.
    gamma : float, optional
        Angle between the a and b lattice parameters, defaults to rectangular
        lattice.
    supercell : tuple of int, optional
        Number of times to repeat the structure in the a and b directions.
    orthorhombic : bool
        If True, will try to generate orthogonal axes for the system. This
        only works for hexagonal and centred rectangular systems. If successful
        the supercell will be halved in the b direction to keep the same number
        of atoms.
    depth : float
        Size of the cell in the direction out of the plane. The system is not
        periodic in this direction, so this will only matter calculations
        that cannot handle 2D periodicity.

    Returns
    -------
    ase.atoms.Atoms
        sheet structure
    """

    species = [Element(species)]

    # Auto square
    if blat is None:
        blat = alat

    # Supercell always takes a sequence of three values, third is always 1
    if isinstance(supercell, int):
        supercell = (supercell, supercell, 1)
    else:
        supercell = (supercell[0], supercell[1], 1)

    # These work for any lattice
    cell = [[alat, 0.0, 0.0],
            [blat*cos(radians(gamma)), blat*sin(radians(gamma)), 0.0],
            [0.0, 0.0, depth]]
    scaled_positions = [(0.0, 0.0, 0.0)]

    # Overwrite the lattice if we prefer an orthogonal lattice
    if orthorhombic and gamma != 90.0:
        # if parallel component of b is half of a, can make orthogonal
        if abs(2*blat*cos(radians(gamma))) - abs(alat) < 1e-5:
            cell = [[alat, 0.0, 0.0],
                    [0, 2*blat*sin(radians(gamma)), 0.0],
                    [0.0, 0.0, depth]]
            scaled_positions = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]
            species *= 2
            if supercell[1] > 1:
                # Halve the supercell in the b direction to keep
                # approximately the same number of atoms (round up)
                supercell = (supercell[0], int((supercell[1]+1)/2.0), 1)
        else:
            print("I don't know how to make this orthogonal "
                  "a={0}, b={1}, gamma={2}".format(alat, blat, gamma))

    new_sheet = Atoms(symbols=species, scaled_positions=scaled_positions,
                      cell=cell, pbc=(True, True, False))

    super_sheet = new_sheet.repeat(supercell)

    return super_sheet


def disc(species, radius, alat, blat=None, gamma=90.0, margin=10.0):
    """
    Create a disc of atoms as an isolated cluster.

    Parameters
    ----------
    species : int or str
        Identity of species used to construct the sheet.
    radius : float
        Radius of disc of atoms to create, atoms on the lattice outside the
        radius will be removed.
    alat : float
        Length of the a lattice parameter.
    blat : float, optional
        Length of the b lattice parameter, if this is not given it defaults
        to the same value as a.
    gamma : float, optional
        Angle between the a and b lattice parameters, defaults to rectangular
        lattice.
    margin : float
        Margin to include around the atoms in all directions

    Returns
    -------
    ase.atoms.Atoms
        disc structure
    """

    species = [Element(species)]

    # Auto square
    if blat is None:
        blat = alat

    # These work for any lattice
    cell = [[alat, 0.0, 0.0],
            [blat*cos(radians(gamma)), blat*sin(radians(gamma)), 0.0],
            [0.0, 0.0, margin]]
    scaled_positions = [(0.0, 0.0, 0.0)]

    # Overwrite the lattice if we prefer an orthogonal lattice
    if gamma != 90.0:
        # if parallel component of b is half of a, can make orthogonal
        if abs(2*blat*cos(radians(gamma))) - abs(alat) < 1e-5:
            cell = [[alat, 0.0, 0.0],
                    [0, 2*blat*sin(radians(gamma)), 0.0],
                    [0.0, 0.0, margin]]
            scaled_positions = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]
            species *= 2

    new_disc = Atoms(symbols=species, scaled_positions=scaled_positions,
                     cell=cell, pbc=(False, False, False))

    # Multiply by 2.1 to get any at the tip
    repeat = (int(ceil(2.2*radius/new_disc.cell[0, 0])),
              int(ceil(2.2*radius/new_disc.cell[1, 1])), 1)

    new_disc = new_disc.repeat(repeat)

    print("Before atoms: {}".format(len(new_disc)))
    # Remove atoms outside disc
    for idx in range(len(new_disc)-1, -1, -1):
        atom = new_disc[idx]
        if ((atom.x-(radius*1.1))**2 + (atom.y-(radius*1.1))**2) > radius**2:
            del new_disc[idx]

    print("After atoms: {}".format(len(new_disc)))
    new_disc.set_cell([2*(radius+margin), 2*(radius+margin), margin])

    return new_disc
