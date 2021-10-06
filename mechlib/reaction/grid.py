""" Build the grid for a transition state search
"""

import automol
import autofile


# Main callable function
def grid_maximum_zmatrices(typ, ts_zma, scan_grids, scan_names, scn_save_fs,
                           mod_thy_info, constraint_dct,
                           series='max', include_endpts=True):
    """ Parses grid(s) of points run along a reaction
        coordinates for the maxima to be able to return guess Z-Matrices
        used in subsequent saddle point optimizations.

        Currently, all reaction scans are ran along a 1D grid, with the
        exception of elimination reactions, which are ran with a 2D grid.

        :param typ: reaction class type
        :type typ: str
        :param grid: set of points that comprise the grid
        :type grid: tuple(numpy.ndarray)
        :param ts_zma: initial ts_zma from ???
        :type ts_zma:
        :param scan_name: name of coordinate in zma along which scan conducted
        :type scan_name: str
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
    """

    print('Assessing scan for a potential maxima...')
    if typ != automol.par.ReactionClass.Typ.ELIMINATION:
        grid, name = scan_grids[0], scan_names[0]
        max_zmas = _find_max_1d(typ, grid, ts_zma, name,
                                mod_thy_info, scn_save_fs, constraint_dct,
                                series=series, include_endpts=include_endpts)
    else:
        grid1, grid2 = scan_grids
        name1, name2 = scan_names
        max_zmas = _find_max_2d(grid1, grid2, name1, name2,
                                mod_thy_info, scn_save_fs, constraint_dct)

    return max_zmas


# Max Finder Functions
def _find_max_1d(typ, grid, ts_zma, scan_name,
                 mod_thy_info, scn_save_fs, constraint_dct,
                 series='max', include_endpts=True):
    """ Parses a one-dimensional grid of points run along a reaction
        coordinates for the maxima to be able to return guess Z-Matrices
        used in subsequent saddle point optimizations.

        For hydrogen migrations, we return an additional zma that involves
        taking the original guess zma and setting the value of the
        scan coordinate to its value at the maximum of the grid. Deals
        with unphysical coords caused by serial optimizations.

        Also return a sequence of zmas {Z(n): 0 <= n <= nmax} where Z(n)
        is a Z-matrix at a point of a scan grid and nmax is index for where
        energy of the grid is maximum.

        :param typ: reaction class type
        :type typ: str
        :param grid: set of points that comprise the grid
        :type grid: tuple(numpy.ndarray)
        :param ts_zma: initial ts_zma from ???
        :type ts_zma:
        :param scan_name: name of coordinate in zma along which scan conducted
        :type scan_name: str
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
        :rtype: tuple(automol.zmat object)
    """

    # Get the locs and energies along the grid
    locs_lst, enes_lst = _grid_vals(
        grid, scan_name, scn_save_fs,
        mod_thy_info, constraint_dct)

    # Get the index where the max energy is found
    max_idx = automol.pot.find_max1d(enes_lst, include_endpts=include_endpts)

    if max_idx is not None:
        print(f'Found maximum at {locs_lst[max_idx]}')

        # Build a series of Z-Matrices based on what series is requested
        max_zmas = ()
        if series == 'max':
            # Get zma at maximum
            max_zma = scn_save_fs[-1].file.zmatrix.read(locs_lst[max_idx])
            max_zmas += (max_zma,)

            # Add second guess zma for migrations:
            # ZMA = original guess zma with val of scan coord at max
            if typ == automol.par.ReactionClass.Typ.HYDROGEN_MIGRATION:
                max_grid_val = grid[max_idx]
                mig_zma = automol.zmat.set_values_by_name(
                    ts_zma, {scan_name: max_grid_val})
                max_zmas += (mig_zma,)
        elif series == 'max-n1':
            for idx in range(max_idx+1):
                max_zmas += (scn_save_fs[-1].file.zmatrix.read(locs_lst[idx]),)
            # Flip list to proceed from max(R) in decreasing order
            max_zmas = max_zmas[::-1]
    else:
        print('No maxima found along the potential')
        max_zmas = None

    return max_zmas


def _find_max_2d(grid1, grid2, scan_name1, scan_name2,
                 mod_thy_info, scn_save_fs, constraint_dct):
    """ Parses a two-dimensional grid of points run along a reaction
        coordinates for the maxima to be able to return guess Z-Matrices
        used in subsequent saddle point optimizations.

        Currently only returns the maxima along the scan.
        Place ZMA in list for generality, compatiability of 1D grid finder ret
        Might add a second option of getting a ZMA at some point

        :param grid1: set of points that comprise first dimension of grid
        :type: tuple(numpy.ndarray)
        :param grid2: set of points that comprise second dimension of grid
        :type: tuple(numpy.ndarray)
        :param scan_name1: name of coordinate in zma along dirst dimension
        :type scan_name1: str
        :param scan_name2: name of coordinate in zma along second dimension
        :type scan_name2: str
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
    """

    # Find the maximum along the 2D grid
    enes_lst = []
    locs_lst_lst = []
    for grid_val_j in grid1:
        locs_list = []
        for grid_val_i in grid2:
            if constraint_dct is None:
                locs_list.append([[scan_name1, scan_name2],
                                  [grid_val_j, grid_val_i]])
            else:
                locs_list.append([constraint_dct, [scan_name1, scan_name2],
                                  [grid_val_i, grid_val_j]])
        enes = []
        locs_lst = []
        for locs in locs_list:
            # use filesys ene read call
            if scn_save_fs[-1].exists(locs):
                scn_path = scn_save_fs[-1].path(locs)
                sp_save_fs = autofile.fs.single_point(scn_path)
                enes.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
                locs_lst.append(locs)
        locs_lst_lst.append(locs_lst)
        if enes:
            enes_lst.append(enes)

    max_enes = []
    max_locs = []
    for idx_j, enes in enumerate(enes_lst):
        max_ene = -10000.
        max_loc = ''
        for idx_i, ene in enumerate(enes):
            if ene > max_ene:
                max_ene = ene
                max_loc = locs_lst_lst[idx_j][idx_i]
        max_enes.append(max_ene)
        max_locs.append(max_loc)

    min_ene = 10000.
    locs = []
    for idx_j, ene in enumerate(max_enes):
        if ene < min_ene:
            min_ene = ene
            locs = max_locs[idx_j]

    # Use the max locs to determine the max_zma, ret as tuple
    max_locs = locs
    max_zma = scn_save_fs[-1].file.zmatrix.read(max_locs)

    return (max_zma,)


def _grid_vals(grid, scan_name, scn_save_fs,
               mod_thy_info, constraint_dct):
    """ Build a list of the filesystem locators for each optimized
        structure as well as a list of their corresponding energies.

        :param grid: set of points that comprise the grid
        :type grid: tuple(numpy.ndarray)
        :param scan_name: name of coordinate in zma along which scan conducted
        :type scan_name: str
        :param scn_save_fs: SCAN/CSCAN object with save filesys prefix
        :type scn_save_fs: autofile.fs.scan or autofile.fs.cscan object
        :param mod_thy_info: ???
        :type mod_thy_info: ???
        :param constraint_dct: values of coordinates to constrain during scan
        :type constraint_dct: dict[str: float]
    """

    # Initialize the lists
    locs_lst = []
    enes_lst = []

    # Build the lists of all the locs for the grid
    grid_locs = []
    for grid_val_i in grid:
        if constraint_dct is None:
            grid_locs.append([[scan_name], [grid_val_i]])
        else:
            grid_locs.append([constraint_dct, [scan_name], [grid_val_i]])

    # Get the energies along the grid
    for locs in grid_locs:
        if scn_save_fs[-1].exists(locs):
            scn_path = scn_save_fs[-1].path(locs)
            sp_save_fs = autofile.fs.single_point(scn_path)
            enes_lst.append(sp_save_fs[-1].file.energy.read(mod_thy_info[1:4]))
            locs_lst.append(locs)

    return locs_lst, enes_lst
