"""
io.py

File reading and writing utilities.

"""

from collections import OrderedDict

import numpy as np
from ase.constraints import FixAtoms, FixCartesian
from quippy import Atoms
from quippy.castep import CastepCell, CastepParam

from sats.util import kpoint_spacing_to_mesh

BASE_CELL = {'kpoint_mp_spacing': 0.03,
             'SPECIES_POT': ['W /work/e89/e89/tdd20/Tungsten/W_OTF.usp']}

BASE_PARAM = {'task': 'SinglePoint',
              'calculate_stress': 'true',
              'opt_strategy': 'speed',
              'xc_functional': 'PBE',
              'cut_off_energy': '600.0',
              'grid_scale': '2.0',
              'finite_basis_corr': 'none',
              'smearing_width': '0.1',
              'perc_extra_bands': '100.0',
              'max_scf_cycles': '100',
              'elec_dump_file': 'NULL',
              'num_dump_cycles': '0',
              'backup_interval': '0',
              'write_checkpoint': 'NONE',
              'popn_calculate': 'false',
              'page_wvfns': '0',
              'write_bib': 'false',
              'bs_write_eigenvalues': 'false',
              'write_cell_structure': 'true',
              'write_cif_structure': 'true'}


def castep_write(atoms, filename='default.cell', optimise=False,
                 fix_lattice=False, kpoint_spacing=0.03):
    """Write atoms a castep cell file and param file."""

    atoms = Atoms(atoms)
    local_cell = BASE_CELL.copy()

    # Set some extra parameters
    local_cell['kpoint_mp_spacing'] = kpoint_spacing

    if fix_lattice:
        local_cell['fix_all_cell'] = 'true'

    if filename.endswith('.cell'):
        cell_filename = filename
        param_filename = filename[:-5] + '.param'
    else:
        cell_filename = filename + '.cell'
        param_filename = filename + '.param'

    local_param = BASE_PARAM.copy()
    if optimise:
        local_param['task'] = 'GeometryOptimization'

    cell = CastepCell(cellfile=local_cell, atoms=atoms)
    cell.update_from_atoms(atoms)
    cell.write(cell_filename)

    param = CastepParam(paramfile=local_param, atoms=atoms)
    param.write(param_filename)


def espresso_write(structure, prefix=None, kpoint_spacing=0.03,
                   custom_pwi=None):
    """Write a simple espresso .pwi file. Not very custom at the moment"""

    # File to write
    pwi = []

    if prefix is None:
        prefix = structure.info['name']

    # mass, pseudopot, magnetisation
    default_pseudo = {
        'Fe': 'Fe.pbe-spn-rrkjus_psl.0.2.1.UPF',
        'H': 'H.pbe-rrkjus_psl.0.1.UPF'}

    default_magmom = {
        'Fe': 0.32,
        'H': 0.0}

    pwi_params = OrderedDict(
        [('CONTROL', OrderedDict([
            ('prefix', prefix),
            ('calculation', 'scf'),
            ('restart_mode', 'from_scratch'),
            ('pseudo_dir', './pseudo/'),
            ('outdir', './{0}_scratch/'.format(prefix)),
            ('verbosity', 'high'),
            ('tprnfor', True),
            ('tstress', True),
            ('disk_io', 'low'),
            ('wf_collect', True),
            ('max_seconds', 82800)])),
        ('SYSTEM', OrderedDict([
            ('ibrav', 0),
            ('nat', len(structure)),
            ('ntyp', len(set(structure.numbers))),
            ('ecutwfc', 90),
            ('ecutrho', 1080),
            ('occupations', 'smearing'),
            ('smearing', 'marzari-vanderbilt'),
            ('degauss', 0.01),
            ('nspin', 2),
            ('nosym', True)])),
        ('ELECTRONS', OrderedDict([
            ('electron_maxstep', 300),
            ('mixing_beta', 0.05),
            ('conv_thr', 1e-9)])),
        ('IONS', OrderedDict([
            ('ion_dynamics', 'bfgs')]))])

    # Nx3 array of constraints matches what QE uses
    constraint_mask = np.ones((len(structure), 3), dtype='int')
    for constraint in structure.constraints:
        if isinstance(constraint, FixAtoms):
            constraint_mask[constraint.index] = 0
        elif isinstance(constraint, FixCartesian):
            constraint_mask[constraint.a] = constraint.mask

    # Make different types for different magnetic moments
    # Rememeber: magnetisation uses 1 based indexes
    atomic_species = OrderedDict()
    atomic_species_str = []
    atomic_positions_str = []
    if 'magmoms' in structure.arrays:
        # Spin has been set manually
        for atom, magmom in zip(structure, structure.arrays['magmoms']):
            if (atom.symbol, magmom) not in atomic_species:
                sidx = len(atomic_species) + 1
                atomic_species[(atom.symbol, magmom)] = sidx
                mag_str = 'starting_magnetization({0})'.format(sidx)
                pwi_params['SYSTEM'][mag_str] = float(magmom)
                atomic_species_str.append(
                    "{species}{sidx} {mass} {pseudo}\n".format(
                        species=atom.symbol, sidx=sidx, mass=atom.mass,
                        pseudo=default_pseudo[atom.symbol]))
            # lookup sidx to append to name
            sidx = atomic_species[(atom.symbol, magmom)]
            atomic_positions_str.append(
                "{atom.symbol}{sidx} "
                "{atom.x:.10f} {atom.y:.10f} {atom.z:.10f} "
                "{mask[0]} {mask[1]} {mask[2]}\n".format(
                    atom=atom, sidx=sidx, mask=constraint_mask[atom.index]))

        # different magnetisms means different types
        pwi_params['SYSTEM']['ntyp'] = len(atomic_species)
    else:
        # Apply default spin for atoms if they need it
        for atom in structure:
            if atom.symbol not in atomic_species:
                sidx = len(atomic_species) + 1
                atomic_species[atom.symbol] = sidx
                if default_magmom[atom.symbol] != 0.0:
                    mag_str = 'starting_magnetization({0})'.format(sidx)
                    pwi_params['SYSTEM'][mag_str] = default_magmom[atom.symbol]
                atomic_species_str.append(
                    "{species} {mass} {pseudo}\n".format(
                        species=atom.symbol, mass=atom.mass,
                        pseudo=default_pseudo[atom.symbol]))
            atomic_positions_str.append(
                "{atom.symbol} "
                "{atom.x:.10f} {atom.y:.10f} {atom.z:.10f} "
                "{mask[0]} {mask[1]} {mask[2]}\n".format(
                    atom=atom, mask=constraint_mask[atom.index]))

    if custom_pwi:
        for section in custom_pwi:
            if section.upper() in pwi_params:
                pwi_params[section.upper()].update(custom_pwi[section])
            else:
                pwi_params[section.upper()] = custom_pwi[section]

    for section in pwi_params:
        pwi.append("&{0}\n".format(section))
        for key, value in pwi_params[section].items():
            if value is True:
                pwi.append('   {0:16} = .true.\n'.format(key))
            elif value is False:
                pwi.append('   {0:16} = .false.\n'.format(key))
            else:
                # repr format to get quotes around strings
                pwi.append('   {0:16} = {1!r:}\n'.format(key, value))
        pwi.append("/\n")
    pwi.append('\n')

    # Pseudopotentials
    pwi.append("ATOMIC_SPECIES\n")
    pwi.extend(atomic_species_str)
    pwi.append("\n")

    # KPOINTS
    pwi.append("K_POINTS automatic\n")
    pwi.append("{0[0]} {0[1]} {0[2]}  1 1 1\n".format(kpoint_spacing_to_mesh(structure, kpoint_spacing)))
    pwi.append("\n")

    # CELL block
    if pwi_params['SYSTEM']['ibrav'] == 0:
        pwi.append("CELL_PARAMETERS angstrom\n")
        pwi.append("{cell[0][0]:.16f} {cell[0][1]:.16f} {cell[0][2]:.16f}\n"
                   "{cell[1][0]:.16f} {cell[1][1]:.16f} {cell[1][2]:.16f}\n"
                   "{cell[2][0]:.16f} {cell[2][1]:.16f} {cell[2][2]:.16f}\n"
                   "".format(cell=structure.cell))
        pwi.append("\n")

    # Positions
    pwi.append("ATOMIC_POSITIONS angstrom\n")
    pwi.extend(atomic_positions_str)
    pwi.append("\n")

    if 'external_force' in structure.arrays:
        pwi.append("ATOMIC_FORCES\n")
        pwi.append(
            '\n'.join(
                ' '.join('{0}'.format(fxyz) for fxyz in force)
                          for force in structure.arrays['external_force']))
        pwi.append('\n')

    with open("{0}.pwi".format(prefix), 'w') as out:
        out.write("".join(pwi))
