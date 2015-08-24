"""
io.py

File reading and writing utilities.

"""

from quippy import Atoms
from quippy.castep import CastepCell, CastepParam

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
                 fix_lattice=False):
    """Write atoms a castep cell file and param file."""

    atoms = Atoms(atoms)
    local_cell = BASE_CELL.copy()
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
