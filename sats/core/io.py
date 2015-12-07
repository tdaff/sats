"""
io.py

File reading and writing utilities.

"""

from collections import OrderedDict

from ase.constraints import FixAtoms
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


def espresso_write(structure, prefix=None, kpoint_spacing=0.05, custom_pwi=None):
    """Write a simple espresso .pwi file. Not very custom at the moment"""

    # File to write
    pwi = []

    if prefix is None:
        prefix = structure.info['name']

    default_species = {
        'Fe': (55.845, 'Fe.pbe-spn-rrkjus_psl.0.2.1.UPF'),
        'H': (1.008, 'H.pbe-rrkjus_psl.0.1.UPF')}

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

    # Make sure that magnetization is applied to the iron
    for idx, species in enumerate(set(atom.symbol for atom in structure), 1):
        if species == 'Fe':
            mag_str = 'starting_magnetization({0})'.format(idx)
            pwi_params['SYSTEM'][mag_str] = 0.32

    if custom_pwi:
        for section in custom_pwi:
            if section.upper() in pwi_params:
                pwi_params[section.upper()].update(custom_pwi[section])
            else:
                pwi_params[section.upper()] = custom_pwi[section]
        #raise NotImplementedError("PWI customisation not implemented yet!")

    for section in pwi_params:
        pwi.append("&{0}\n".format(section))
        for key, value in pwi_params[section].items():
            if value is True:
                pwi.append('   {0:16} = .true.\n'.format(key))
            else:
                # repr format to get quotes around strings
                pwi.append('   {0:16} = {1!r:}\n'.format(key, value))
        pwi.append("/\n")
    pwi.append('\n')

    # Pseudopotentials
    pwi.append("ATOMIC_SPECIES\n")
    for species in set(atom.symbol for atom in structure):
        pwi.append("{0} {1[0]} {1[1]}\n".format(species, default_species[species]))
    pwi.append("\n")

    # KPOINTS
    pwi.append("K_POINTS automatic\n")
    pwi.append("{0[0]} {0[1]} {0[2]}  1 1 1\n".format(kpoint_spacing_to_mesh(structure, kpoint_spacing)))
    pwi.append("\n")

    # CELL block
    pwi.append("CELL_PARAMETERS angstrom\n")
    pwi.append("{cell[0][0]:.16f} {cell[0][1]:.16f} {cell[0][2]:.16f}\n"
               "{cell[1][0]:.16f} {cell[1][1]:.16f} {cell[1][2]:.16f}\n"
               "{cell[2][0]:.16f} {cell[2][1]:.16f} {cell[2][2]:.16f}\n"
               "".format(cell=structure.cell))
    pwi.append("\n")

    # Positions
    mask = [False]*len(structure)
    for constraint in structure.constraints:
        if isinstance(constraint, FixAtoms):
            mask = [any([x, y]) for x, y in zip(mask, constraint.index)]

    pwi.append("ATOMIC_POSITIONS angstrom\n")
    for atom, masked in zip(structure, mask):
        pwi.append("{atom.symbol} {atom.x:.10f} {atom.y:.10f} {atom.z:.10f}"
                   "".format(atom=atom))
        if masked:
            # Fix position
            pwi.append(" 0 0 0")
        pwi.append("\n")
    pwi.append("\n")

    with open("{0}.pwi".format(prefix), 'w') as out:
        out.write("".join(pwi))
