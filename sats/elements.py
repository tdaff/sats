"""
elements.py

Various data relating to different elements.

"""


class Element(object):
    """
    Representation of a single element species.
    """
    def __init__(self, species):
        """
        Interpret species as an element type.

        Parameters
        ----------
        species : int or str
            The identity of the atomic species, either as the atomic number
            or the element symbol.

        Raises
        ------
        ValueError
            Unidentified species
        """

        try:
            self.z = int(species)
        except ValueError:
            # Assume it is an element symbol
            self.z = ATOMIC_NUMBER.index(species)

    @property
    def atomic_number(self):
        """Atomic number."""
        return self.z

    @atomic_number.setter
    def atomic_number(self, value):
        """Atomic number."""
        self.z = int(value)

    @property
    def symbol(self):
        """Element symbol."""
        return ATOMIC_NUMBER[self.z]

    @symbol.setter
    def symbol(self, value):
        """Element symbol."""
        self.z = ATOMIC_NUMBER.index(value)


# If this is a nice list we can just .index or [slice] to get atomic numbers
ATOMIC_NUMBER = [
    "ZERO", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uuo"]