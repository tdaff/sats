"""
structures.py

Helpers to create different starting environments.

"""


from quippy.structures import bcc, fcc, supercell


def bcc_intersitital_111(basis, supercell=(1, 1, 1)):
    """
    Create a crystal structure with an interstitial atom.

    :param basis: structure to use to generate interstital
    :param supercell: generate a supercell before
    :return: quippy atoms object with interstitial atom
    """