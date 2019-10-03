""" utilites
"""
import os
import stat
import subprocess
import warnings
import autofile
import automol
import elstruct


def run_qchem_par(prog, method):
    """ dictionary of parameters for different electronic structure codes
    """
    if prog == 'gaussian09':
        sp_script_str = ("#!/usr/bin/env bash\n"
                         "g09 run.inp run.out >> stdout.log &> stderr.log")
        opt_script_str = sp_script_str
        kwargs = {
            'memory': 20,
            'machine_options': ['%NProcShared=10'],
            'gen_lines': ['# int=ultrafine'],
        }
        opt_kwargs = {
            'memory': 20,
            'machine_options': ['%NProcShared=10'],
            'gen_lines': ['# int=ultrafine'],
            'feedback': True,
#            'job_options': ['verytight'],
#            'job_options': ['verytight'],
            'errors': [
                elstruct.Error.OPT_NOCONV
            ],
            'options_mat': [
                [{},
                 {},
                 {},
                 {'job_options': ['calcfc']},
                 {'job_options': ['calcfc']},
                 {'job_options': ['calcall']}]
            ],
        }

    if prog == 'psi4':
        sp_script_str = ("#!/usr/bin/env bash\n"
                         "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = sp_script_str
        kwargs = {}
        opt_kwargs = {}

    if prog == 'molpro2015':
        sp_script_str = ("#!/usr/bin/env bash\n"
                         "molpro -n 8 run.inp -o run.out >> stdout.log &> stderr.log")
        if method == 'caspt2':
            opt_script_str = ("#!/usr/bin/env bash\n"
                              "molpro -n 8 run.inp -o run.out >> stdout.log &> stderr.log")
        else:
            opt_script_str = ("#!/usr/bin/env bash\n"
                              "molpro --mppx -n 12 run.inp -o run.out >> stdout.log &> stderr.log")
        if method == 'caspt2' or method == 'caspt2c':
            kwargs = {
                'memory': 10,
                'corr_options': ['shift=0.2'],
                'mol_options': ['nosym'],
            }
            opt_kwargs = {
                'memory': 5,
                'corr_options': ['shift=0.2'],
                'mol_options': ['nosym'],
                'options_mat': [
                    [{},
                     {},
                     {'job_options': ['numhess']},
                     {'job_options': ['numhess=10']},
                     {'job_options': ['numhess=1']}]
                ],
            }
        elif method == 'caspt2i':
            kwargs = {
                'memory': 10,
                'corr_options': ['shift=0.2', 'ipea=0.25'],
                'mol_options': ['nosym'],
            }
            opt_kwargs = {
                'memory': 5,
                'corr_options': ['shift=0.2', 'ipea=0.25'],
                'mol_options': ['nosym'],
                'options_mat': [
                    [{},
                     {},
                     {'job_options': ['numhess']},
                     {'job_options': ['numhess=10']},
                     {'job_options': ['numhess=1']}]
                ],
            }
        else:
            kwargs = {
                'memory': 10,
                'mol_options': ['nosym'],
            }
            opt_kwargs = {
                'memory': 5,
                'mol_options': ['nosym'],
                'options_mat': [
                    [{},
                     {},
                     {'job_options': ['numhess']},
                     {'job_options': ['numhess=10']},
                     {'job_options': ['numhess=1']}]
                ],
            }

    if prog == 'qchem':
        sp_script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = sp_script_str
        kwargs = {
            'memory': 20,
        }
        opt_kwargs = {}

    if prog == 'cfour':
        sp_script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = sp_script_str
        kwargs = {
            'memory': 20,
        }
        opt_kwargs = {}

    if prog == 'orca':
        sp_script_str = ("#!/usr/bin/env bash\n"
                      "molpro -i run.inp -o run.out >> stdout.log &> stderr.log")
        opt_script_str = sp_script_str
        kwargs = {
            'memory': 20,
        }
        opt_kwargs = {}

    return sp_script_str, opt_script_str, kwargs, opt_kwargs


def orbital_restriction(spc_info, thy_level):
    """ orbital restriction logical
    """
    mul = spc_info[2]
    if thy_level[3] == 'RR':
        orb_restr = True
    elif thy_level[3] == 'UU':
        orb_restr = False
    elif thy_level[3] == 'RU':
        if mul == 1:
            orb_restr = True
        else:
            orb_restr = False
    return orb_restr


def geometry_dictionary(geom_path):
    """ read in dictionary of saved geometries
    """
    geom_dct = {}
    for dir_path, _, file_names in os.walk(geom_path):
        for file_name in file_names:
            file_path = os.path.join(dir_path, file_name)
            if file_path.endswith('.xyz'):
                xyz_str = autofile.file.read_file(file_path)
                geo = automol.geom.from_xyz_string(xyz_str)
                ich = automol.geom.inchi(geo)
                if ich in geom_dct:
                    print('Warning: Dupilicate xyz geometry for ', ich)
                geom_dct[ich] = geo
    return geom_dct


def min_energy_conformer_locators(cnf_save_fs):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs.leaf.existing()
    if cnf_locs_lst:
        cnf_enes = [cnf_save_fs.leaf.file.energy.read(locs)
                    for locs in cnf_locs_lst]
        min_cnf_locs = cnf_locs_lst[cnf_enes.index(min(cnf_enes))]
    else:
        min_cnf_locs = None
    return min_cnf_locs


def min_dist_conformer_zma(dist_name, cnf_save_fs):
    """ locators for minimum energy conformer """
    cnf_locs_lst = cnf_save_fs.leaf.existing()
    cnf_zmas = [cnf_save_fs.leaf.file.zmatrix.read(locs)
                for locs in cnf_locs_lst]
    min_dist = 100.
    min_zma = []
    for zma in cnf_zmas:
        dist = automol.zmatrix.values(zma)[dist_name]
        if dist < min_dist:
            min_dist = dist
            min_zma = zma
    min_zma = [min_zma]
    return min_zma


def locs_sort(save_fs):
    """ sort trajectory file according to energies
    """
    locs_lst = save_fs.leaf.existing()
    if locs_lst:
        enes = [save_fs.leaf.file.energy.read(locs)
                for locs in locs_lst]
        sorted_locs = []
        for _, loc in sorted(zip(enes, locs_lst), key=lambda x: x[0]):
            sorted_locs.append(loc)
    return sorted_locs


def traj_sort(save_fs):
    """ sort trajectory file according to energies
    """
    locs_lst = save_fs.leaf.existing()
    if locs_lst:
        enes = [save_fs.leaf.file.energy.read(locs)
                for locs in locs_lst]
        geos = [save_fs.leaf.file.geometry.read(locs)
                for locs in locs_lst]
        traj = []
        for ene, geo, locs in sorted(zip(enes, geos, locs_lst), key=lambda x: x[0]):
            comment = 'energy: {0:>15.10f} \t {1}'.format(ene, locs[0])
            traj.append((comment, geo))
        traj_path = save_fs.trunk.file.trajectory.path()
        print("Updating trajectory file at {}".format(traj_path))
        save_fs.trunk.file.trajectory.write(traj)


def nsamp_init(nsamp_par, ntaudof):
    """ determine nsamp for given species"""
    if nsamp_par[0]:
        nsamp = min(nsamp_par[1] + nsamp_par[2] * nsamp_par[3]**ntaudof,
                    nsamp_par[4])
    else:
        nsamp = nsamp_par[5]
    return nsamp


def reaction_energy(save_prefix, rxn_ich, rxn_chg, rxn_mul, thy_level):
    """ reaction energy """
    rct_ichs, prd_ichs = rxn_ich
    rct_chgs, prd_chgs = rxn_chg
    rct_muls, prd_muls = rxn_mul
    rct_enes = reagent_energies(
        save_prefix, rct_ichs, rct_chgs, rct_muls, thy_level)
    prd_enes = reagent_energies(
        save_prefix, prd_ichs, prd_chgs, prd_muls, thy_level)
    return sum(prd_enes) - sum(rct_enes)


def reagent_energies(save_prefix, rgt_ichs, rgt_chgs, rgt_muls, thy_level):
    """ reagent energies """
    enes = []
    for rgt_ich, rgt_chg, rgt_mul in zip(rgt_ichs, rgt_chgs, rgt_muls):
        spc_save_fs = autofile.fs.species(save_prefix)
        rgt_info = [rgt_ich, rgt_chg, rgt_mul]
        spc_save_path = spc_save_fs.leaf.path(rgt_info)

        orb_restr = orbital_restriction(rgt_info, thy_level)
        thy_lvl = thy_level[0:3]
        thy_lvl.append(orb_restr)
        thy_save_fs = autofile.fs.theory(spc_save_path)
        thy_save_path = thy_save_fs.leaf.path(thy_lvl[1:4])
        cnf_save_fs = autofile.fs.conformer(thy_save_path)
        min_cnf_locs = min_energy_conformer_locators(cnf_save_fs)
        ene = cnf_save_fs.leaf.file.energy.read(min_cnf_locs)
        enes.append(ene)
    return enes


def ts_mul_from_reaction_muls(rcts, prds, spc_dct):
    """
    evaluate the ts multiplicity from the multiplicities of the reactants and products
    """
    nrcts = len(rcts)
    nprds = len(prds)
    rct_spin_sum = 0
    prd_spin_sum = 0
    rad_rad = True
    rct_muls = []
    prd_muls = []
    if nrcts == 1 and nprds == 1:
        ts_mul = max(spc_dct[rcts[0]]['mul'], spc_dct[prds[0]]['mul'])
        rad_rad = False
    else:
        for rct in rcts:
            rct_spin_sum += (spc_dct[rct]['mul'] - 1.)/2.
            rct_muls.append(spc_dct[rct]['mul'])
        for prd in prds:
            prd_spin_sum += (spc_dct[prd]['mul'] - 1.)/2.
            prd_muls.append(spc_dct[prd]['mul'])
        if (min(rct_muls) == 1 or nrcts == 1) and (min(prd_muls) == 1 or nprds == 1):
            rad_rad = False
        ts_mul = min(rct_spin_sum, prd_spin_sum)
        ts_mul = int(round(2*ts_mul + 1))
    return ts_mul, rad_rad


def run_script(script_str, run_dir):
    """ run a program from a script
    """

    script_name = 'build.sh'
    with _EnterDirectory(run_dir):
        # write the submit script to the run directory
        try:
            os.remove('build.sh')
            with open(script_name, 'w') as script_obj:
                script_obj.write(script_str)
        except:
            with open(script_name, 'w') as script_obj:
                script_obj.write(script_str)

        # make the script executable
        os.chmod(script_name, mode=os.stat(script_name).st_mode | stat.S_IEXEC)

        # call the program
        try:
            subprocess.check_call('./{:s}'.format(script_name))
        except subprocess.CalledProcessError as err:
            # if the program failed, continue with a warning
            warnings.warn("run failed in {}".format(run_dir))


class _EnterDirectory():

    def __init__(self, directory):
        assert os.path.isdir(directory)
        self.directory = directory
        self.working_directory = os.getcwd()

    def __enter__(self):
        os.chdir(self.directory)

    def __exit__(self, _exc_type, _exc_value, _traceback):
        os.chdir(self.working_directory)
