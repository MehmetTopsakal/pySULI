import os
import sys
import shutil
import scipy
import time
import glob
import ipympl
import matplotlib.pyplot as plt

import numpy as np
import xarray as xr

import random
import string






import functions_lib_version

from IPython.display import clear_output





# Tools for navigating the gpx object revcieved from GSAS-II analysis

def print_dict_as_tree(struct, indent=0, depth=0):
    """
    Recursively navigate a dictionary and print it as a tree structure.

    Args:
        struct (dict_like): dictionary-like structure to be printed
        indent (int): current level of indentation, use as a "memory" of recursive depth
        depth (int): current depth in the dictionary, 0 represents the first level
    """
    for key, value in struct.items():
        print(' ' * indent + str(depth) + ". " + str(key))
        if isinstance(value, dict):
            print_dict_as_tree(value, indent + 4, depth + 1)
        else:
            print(' ' * (indent + 4) + str(depth) + ". " + str(value))

def find_key(struct, key):
    """
    Recursively navigate a dictionary-like structure and search for a key. This
    doesn't work great if the key is redundant, but is ideal for keys which occur
    only once in the gpx structure.

    Args:
        struct (dict_like): any dictionary, works with GPX
        key (str): the key to find

    Returns:
        result (unknown): the value of the key
    """
    if isinstance(struct, dict):
        if key in struct:
            return struct[key]
        else:
            for k, v in struct.items():
                result = find_key(v, key)
                if result is not None:
                    return result
    elif isinstance(struct, list):
        for item in struct:
            result = find_key(item, key)
            if result is not None:
                return result

def find_allof_key(d, key, results=None, depth=0, min_depth=0, max_depth=np.inf):
    """
    Recursively navigate a dictionary-like structure and search for each instance of a
    key. This returns an array of all the structures found stemming from the key. If
    it reaches a structure at the lowest level (ie, no subdictionaries are in the
    structure), containing the key, returns the structure containing the key. The
    min_depth and max_depth thresholds are useful if there are many results, and you
    need to filter them.

    Args:
        d (dict_like): any dictionary, works with GPX
        key (str): the key to find
        results (array_like): tuple list which grows with each instance of key
        depth (int): depth of a given appended structure
        min_depth (int): sets a minimum depth for collected structures
        min_depth (int): analogous to min depth

    Returns:
        results (array_like): list of tuples, where each first index represents tells the depths, and the second index is the structure found; to access an output structue structure at a given index, do `find_allof_key(...)[index][1]`
    """
    # instantiates results list
    if results is None:
        results = []

    # checks if we have a dict; if so looks for target key in dict, and continues recursively searching dict
    if isinstance(d, dict):
        if key in d and depth >= min_depth and depth <= max_depth:
            results.append(("Depth " + str(depth), d))
        for key, value in d.items():
            if isinstance(value, dict):
                find_allof_key(value, key, results, depth=depth+1, min_depth=min_depth, max_depth=max_depth)
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        find_allof_key(item, key, results, min_depth=min_depth, max_depth=max_depth)

    # checks if we have a list; if so looks for dictionaries in the list and searches them
    elif isinstance(d, list):
        for item in d:
            if isinstance(item, dict):
                find_key(item, key, results, min_depth=min_depth, max_depth=max_depth)

    # trys to perform search on dictionary-like object or list
    else:
        try:
            if key in d and depth >= min_depth and depth <= max_depth:
                results.append(("Depth " + str(depth), d))
            for key, value in d.items():
                if isinstance(value, dict):
                    find_key(value, key, results, depth=depth+1, min_depth=min_depth, max_depth=max_depth)
                elif isinstance(value, list):
                    for item in value:
                        if isinstance(item, dict):
                            find_key(item, key, results, depth=depth+1, min_depth=min_depth, max_depth=max_depth)
        except:
            print("Not a valid dictionary.")
            return 0
    return results

def get_valid_phases(gpx):
    """
    Args:
        gpx (dict_like): any dictionary, but meant for GPX

    Returns:
        phases (array_like): list of strs represent the keys corresponding to each phase fitted in the GSAS-II refinement
    """
    try:
        phases = gpx['Phases']
        return list(phases.keys())[1:]
    except:
        print("Invalid gpx object.")
        return 0

def get_cell_consts(gpx, phase, include_angles=False):
    """
    Navigates a gpx data output to search for a phase, then creates and
    returns a dictionary containing lattice constant information.

    Args:
        gpx (dict_like): any dictionary, but meant for GPX
        phase (str): the phase tag
        include_angles (bool, optional): flag which specified whether to include angle constants. Defaults to False.

    Returns:
        const_dict (dict_like): dictionary with lattice constants as keys
    """
    try:
        phase_info = gpx['Phases'][phase]
        constants = find_key(phase_info, 'Cell')
    except:
        print("Invalid gpx object.")
        return 0
    const_dict = {}
    const_dict['a'] = constants[1]
    const_dict['b'] = constants[2]
    const_dict['c'] = constants[3]
    if include_angles:
        const_dict['alpha'] = constants[4]
        const_dict['beta'] = constants[5]
        const_dict['gamma'] = constants[6]
    return const_dict

def print_all_cell_consts(gpx):
    """
        Prints out each phase, and its respective lattice constants.

        gpx (dict_like): any dictionary, but meant for GPX
    """
    valid_phases = get_valid_phases(gpx)
    for phase in valid_phases:
        print(phase)
        cell_constants = get_cell_consts(gpx, phase)
        for const in cell_constants:
            print(' ' * 4 + const + ':')
            print(' ' * 8 + str(cell_constants[const]))


def gpx_plotter(
    gpx,
    line_axes,
    stem_axes,
    phases,
    radial_range,
    marker="o",
    stem=True,
    unit="d",
    plot_unit="q_A^-1",
    y_shift=0.1
):
    """
    Methods which plots intensity points at the given d location (X) with
    normalized intensity (Y). It gets these values from the gpx object
    created by the GSAS-II extender functions.

    Args:
        gpx (dict_like): GSAS project object, has dictionary structure
        line_axes (_type_): axes where to plot points representing peaks
        stem_axes (_type_): axes where to plot vertical lines representing peaks
        phases (_type_): string array of phases (can get from get_valid_phases(GPX))
        radial_range (_type_): range for the x axis.
        marker (str, optional): marker style. Defaults to "o".
        stem (bool, optional): whether to add stem plot. Defaults to True.
        unit (str, optional): the unit of the peak positions in gpx file. Defaults to "d".
        plot_unit (str, optional): the unit for x axis. Defaults to "q_A^-1".
        y_shift (float, optional): y shift between phase labels. Defaults to 0.1.
    """
    for ep, phase in enumerate(phases):
        # getting intensity peak positions (X) and relative magnitudes (Y)
        X = np.array([])
        Y = np.array([])
        for value in gpx['PWDR data.xy']['Reflection Lists'][phase]['RefList']:
            # takes the d column from the Reflection Lists dataframe-like array
            X = np.append(X, value[4])
        # Calculates the intensity maxes and normalized with 100 being the max
        peaks = gpx['PWDR data.xy']['Reflection Lists'][phase]
        Y = peaks['RefList'].T[8]*np.array([refl[11] for refl in peaks['RefList']])
        Imax = np.max(Y)
        if Imax:
            Y *= 100.0/Imax

        if unit == plot_unit:
            continue
        elif unit=="2th_deg" and plot_unit=="q_A^-1":
            X = twotheta_to_q(X, gpx['PWDR data.xy']['Instrument Parameters'][0]['Lam'][ep])
        elif unit=="d" and plot_unit=="q_A^-1":
            X = d_to_q(X)

        label=phase
        label_x=1.02
        color="C%d" % ep
        label_y=ep * y_shift

        for axl in line_axes:
            for i in X:
                axl.axvline(x=i, lw=0.6, color=color)
                axl.set_xlim([radial_range[0], radial_range[1]])

        for axs in stem_axes:
            axs_stem = axs.twinx()
            if stem:
                markerline, stemlines, baseline = axs_stem.stem(X, Y, markerfmt=".")
                plt.setp(stemlines, linewidth=0.5, color=color)
                plt.setp(markerline, color=color)
            axs_stem.set_xlim([radial_range[0], radial_range[1]])
            axs_stem.set_yticks([])
            axs_stem.set_ylim(bottom=0.1)
            axs_stem.text(
                label_x, label_y, label, color=color, transform=axs_stem.transAxes
            )











class HiddenPrints:
    """
    This class hides print outputs from functions. It is useful for processes like refinement which produce a lot of text prints.
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def set_GSAS_II_user(user):
    """
    This method allows you to append your GSAS-II directory to your python path. It is not used in the Refiner class below because you would have to modify the pySULI source code yourself in order to use this method. However, it is useful if you have multiple users, as you can simply match cases as shown below.

    Args:
        user (str): string representing the user to match; you can see the 3 example strings of the developers below, but you would have to add cases as needed

    Returns:
        G2sc, pybaselines: pointers for using GSASIIscriptable and pybaslines modules
    """
    match user:
        case "mtopsakal":
            sys.path += ['/home/mt/software/miniforge3/envs/GSASII/GSAS-II/GSASII']
        case "pmeshkov":
            sys.path += ['/opt/anaconda3/envs/GSASII/GSAS-II/GSASII']
        case "kmorell":
            sys.path += ['/Users/kevinmorell/Downloads/anaconda3/envs/GSASII/GSAS-II/GSASII']
        case "user":
            sys.path += ['your path here']

    # we then import GSASIIscriptable
    import GSASIIscriptable as G2sc
    import pybaselines # this comes with gsas2_package

    return G2sc, pybaselines

def set_GSAS_II_path(GSASII_path):
    """
    Set the path to the GSASII package.

    Args:
        GSASII_path (str): string representing the path to the GSASII package; should be of the form [your path]/GSAS-II/GSASII

    Returns:
        G2sc, pybaselines: pointers for the GSASIIscriptable and pybaslines modules
    """
    sys.path += [GSASII_path]

    # we then import GSASIIscriptable
    import GSASIIscriptable as G2sc
    import pybaselines # this comes with gsas2_package

    return G2sc, pybaselines

class Refiner:
    def __init__(self,
                 nc_path,
                 phases,
                 instrument_parameters_file,
                 gsas2_scratch = None,
                 GSASII_path = None,
                 da_input_bkg = None,
                 q_range = (1,10),
                 verbose = False
                 ):
        """
        This is the refiner object. It allows you to match certain phases (from their corresponding cif files) to a nc file representing the diffraction pattern. It also requires a path to an instrument parameters files which determined the initial instrument parameters. The gsas2_scratch parameter represents the directory where the refiner will create and modify files (data.xy, gsas.bak0.gpx, gsas.gpx, gsas.instprm, gsas.lst) throughout its refinement process. The rest of the parameters will be described below. 
        
        Upon instantiation, basic processing of the nc file, and then background fitting is conducted. Then the gpx and ds objects/files are created. 

        Args:
            nc_path (str): the path to your starting measurment file, the '.nc' file
            phases (dict_array): an array of phases, each element a dictionary of the form {'cif_abs_path':'LaB6.cif','phase_name':'LaB6','scale':1}
            instrument_parameters_file (str): a string representing the path to the '.instprm' file
            gsas2_scratch (str, optional): See description above for more details, represents path to gsas2_scratch directory. Defaults to None, and if default, will create a GSAS-II directory on its own.
            GSASII_path (str, optional): a string representing the path of the GSAS-II module containing directory, should be of the form '[your path]/GSASII/GSAS-II/GSASII'. Defaults to None, only use default if you have your GSAS-II path set up, in a way that you can directly import GSAS-II modules using 'import [module]'.
            da_input_bkg (xarray.DataArray, optional): xarray DataArray representing the expected background spectrum. Defaults to None, only use default if the background spectrum is unknown; providing an expected background spectrum will yield better refinement.
            q_range (tuple, optional): Range of q values to consider (in units 1/A0, although it depends on your data units). Defaults to (1,10).
        Returns:
            Refiner: Refiner object with its associated functions, as seen following the initialization function.
        """

        if GSASII_path is not None:
            G2sc, pybaselines = set_GSAS_II_path(GSASII_path)
        else:
            try:
                import GSASIIscriptable as G2sc
                import pybaselines
            except:
                raise ModuleNotFoundError("GSAS-II path not found.\n Install GSAS-II, and insert the GSAS-II_path into the refiner class.")


        if gsas2_scratch is None:
            if not os.path.exists('.gsas2_scratch'):
                os.makedirs('.gsas2_scratch',exist_ok=True)
            randstr = ''.join(random.choices(string.ascii_uppercase+string.digits, k=7))
            gsas2_scratch = '%s/%.2f_%s.tmp'%('.gsas2_scratch',time.time(),randstr)
            os.makedirs(gsas2_scratch)
            self.gsas2_scratch = gsas2_scratch
        else:
            self.gsas2_scratch = gsas2_scratch


        # set instance variables
        self.nc_path = nc_path
        self.phases = phases
        self.instrument_parameters_file = instrument_parameters_file
        self.q_range = q_range
        self.da_input_bkg = da_input_bkg
        self.da_i2d = None
        self.da_i2d_m = None
        self.bkg_auto = None
        self.name = nc_path
        self.verbose = verbose


        with xr.open_dataset(self.nc_path) as ds:

            #auto set q_range
            if self.q_range is None:
                self.q_range = [float(ds.X_in_q.min()), float(ds.X_in_q.max())]
                
            for k in ['Y_obs','Y_calc','Y_bkg_auto','Y_bkg_gsas','Y_bkg','Y_bkg_arpls','Y_bkg_auto']:
                if k in ds.keys():
                    del ds[k]

            self.da_i2d = ds.i2d.sel(radial=slice(self.q_range[0],self.q_range[1])).astype('float32')
            self.da_i2d_m = self.da_i2d.mean(dim='azimuthal')

            # set background function
            if self.da_input_bkg is None:
                arpls_ = pybaselines.Baseline(x_data=self.da_i2d_m.radial.values).arpls((self.da_i2d_m).values, lam=1e5)[0]
                shift_ = min((self.da_i2d_m).values - arpls_)
                bkg_arpls = (arpls_+shift_)
                self.bkg_auto = bkg_arpls

            else:
                self.da_input_bkg = self.da_input_bkg.sel(radial=slice(self.q_range[0],self.q_range[1]))
                blank_scale = (self.da_i2d_m[0] / self.da_input_bkg[0]).values
                while (min((self.da_i2d_m.values-blank_scale*self.da_input_bkg.values)) < 0):
                    blank_scale = blank_scale*0.95
                da_input_bkg_scaled = blank_scale*self.da_input_bkg
                
                # sometimes the x_data values (the radial values) of the data and the background disagree; if each disagreement is less than 10^-6, we will ignore the dissagreement, and will treat replace the radial values of the background with those of the data(the choise to replace those of the background is arbitrary)
                try:
                    arpls_ = pybaselines.Baseline(x_data=self.da_i2d_m.radial.values).arpls((self.da_i2d_m-da_input_bkg_scaled).values, lam=1e5)[0]
                except:
                    diffs = self.da_i2d_m.radial.data != da_input_bkg_scaled.radial.data
                    diff_data = self.da_i2d_m.radial.data[diffs]
                    diff_bkg = da_input_bkg_scaled.radial.data[diffs]
                    if abs(max(diff_data - diff_bkg)) < 0.000001:
                        da_input_bkg_scaled = da_input_bkg_scaled.assign_coords(radial=self.da_i2d_m.radial)
                        arpls_ = pybaselines.Baseline(x_data=self.da_i2d_m.radial.values).arpls((self.da_i2d_m-da_input_bkg_scaled).values, lam=1e5)[0]
                    else:
                        raise ValueError("The background's radial values are two different from the observed data's radial values (there is at least 1 point with a difference above 10e-6)")
                shift_ = min((self.da_i2d_m-da_input_bkg_scaled).values - arpls_)
                bkg_arpls = (arpls_+shift_)
                self.bkg_auto = bkg_arpls + da_input_bkg_scaled.values

            Y_to_gsas = self.da_i2d_m.values-self.bkg_auto

            self.y_scale = 1000/max(Y_to_gsas)
            self.y_baseline = 10

            # set gsas intensity profile and save as xy file in gsas2_scratch directory
            Y_to_gsas = self.y_baseline+Y_to_gsas*self.y_scale
            X_to_gsas = np.rad2deg(functions_lib_version.q_to_twotheta(self.da_i2d_m.radial.values, wavelength=(ds.attrs['wavelength']*1.0e10)))
            np.savetxt('%s/data.xy'%self.gsas2_scratch, np.column_stack( (X_to_gsas,Y_to_gsas) ), fmt='%.4f %.4f')

            # save ds as an instance variable
            self.ds = ds

            # save gpx as an instance variable
            self.gpx = G2sc.G2Project(newgpx='%s/gsas.gpx'%self.gsas2_scratch)
            self.gpx.data['Controls']['data']['max cyc'] = 100
            self.gpx.add_powder_histogram('%s/data.xy'%self.gsas2_scratch,instrument_parameters_file)

            hist = self.gpx.histograms()[0]
            for p in phases:
                self.gpx.add_phase(p['cif_abs_path'],phasename=p['phase_name'], histograms=[hist],fmthint='CIF')
            self.gpx.save()

    def save_refinement(self):
        """
        This method is to be run after the refinements are conducted (ex, after myrefiner.refine_background(5)), and it saves all the refinement results to the ds files. The call of this method is necessary for the plotting methods to work.
        
        Returns:
            array: a tuple (length 2 array), with the gpx G2Project object and the ds_list array [self.gpx, self.ds_list].
        """
        # apply refinement recipy to the instance variable gpx
        # possibly make this simpler later

        # get x, yobs, ycalc, and background from gpx histograms
        histogram = self.gpx.histograms()[0]

        gsas_X  = histogram.getdata('x').astype('float32')
        gsas_Yo = histogram.getdata('yobs').astype('float32')
        gsas_Yc = histogram.getdata('ycalc').astype('float32')
        gsas_B  = histogram.getdata('Background').astype('float32')

        # now convert Y back
        Y_Obs_from_gsas  = ((gsas_Yo -self.y_baseline  )/self.y_scale )
        Y_Bkg_from_gsas  = ((gsas_B  -self.y_baseline  )/self.y_scale )
        Y_calc_from_gsas = ((gsas_Yc -self.y_baseline  )/self.y_scale ) # this also includes background from GSAS

        # add data to the instance variable ds
        self.ds = self.ds.assign_coords(
            {"X_in_q": self.da_i2d_m.radial.values.astype('float32')},
            )
        self.ds = self.ds.assign_coords(
            {"X_in_tth": gsas_X.astype('float32')},
            )
        self.ds = self.ds.assign_coords(
            {"X_in_d": functions_lib_version.q_to_d(self.da_i2d_m.radial.values.astype('float32'))},
            )

        self.ds['Y_obs'] = xr.DataArray(
                                    data=self.da_i2d_m.values,
                                    dims=['X_in_q'],
                                )
        self.ds['Y_calc'] = xr.DataArray(
                                    data=Y_calc_from_gsas-Y_Bkg_from_gsas,
                                    dims=['X_in_q'],
                                )# by Ycalc, we mean only calculated peaks, no background

        self.ds['Y_bkg_gsas'] = xr.DataArray(
                                    data=Y_Bkg_from_gsas,
                                    dims=['X_in_q'],
                                )
        self.ds['Y_bkg_auto'] = xr.DataArray(
                                    data=self.bkg_auto,
                                    dims=['X_in_q'],
                                )
        file_path = os.path.join(self.gsas2_scratch, 'gsas.lst')

        # Check if the directory and file exist
        if not os.path.exists(self.gsas2_scratch):
            print(f"Directory does not exist: {self.gsas2_scratch}")
        elif self.verbose:
            print(f"Directory exists: {self.gsas2_scratch}")

        if not os.path.isfile(file_path):
            print(f"File not found: {file_path}")
        else:
            # Open and read the file if it exists
            with open(file_path) as lst_file:
                gsas_lst = lst_file.read()
                self.ds.attrs['gsas_lst'] = gsas_lst
        
        #with open('%s/gsas.lst'%(self.gsas2_scratch)) as lst_file:
        #    gsas_lst = lst_file.read()
        #    self.ds.attrs['gsas_lst'] = gsas_lst

        return [self.gpx,self.ds]

    def update_ds_file(self):
        """
        Updates the ds file in nc_path
        """
        self.ds.to_netcdf(self.nc_path+'.new.nc',engine="scipy")
        time.sleep(0.1)
        shutil.move(self.nc_path+'.new.nc',self.nc_path)

    def update_gpx(self):
        """
        Updates the gpx file in nc_path
        """
        shutil.copyfile(self.gpx['Controls']['data']['LastSavedAs'],'%s.gpx'%self.nc_path[:-3])

    def get_gpx(self):
        """
        Returns:
            G2Project: gpx project object storing the relevant variables for refinement
        """
        return self.gpx

    def get_ds(self):
        """
        Returns:
            nCDF: nCDF object storing the experimental data and refinement fits
        """
        return self.ds

    def get_y_vals(self):
        """
        Returns:
            X, Yobs, Ycalc, Ybkg: returns four xarray DataArrays:
            0. X (q values of data point for a given index)
            1. Yobs (observed intensity of data point for a corresponding index in the X array)
            2. Ycalc (calculated intensity attributed to the bragg diffraction peaks of a data point given corresponding index in the X array)
            3. Ybkg (calculated intensity attributed to the background of a data point given corresponding index in the X array; note that to get the full intensity profile, fitted to Yobs, you must sum Ycalc and Ybkg)
        """
        X, Yobs, Ycalc, Ybkg= self.ds['X_in_q'].values, self.ds['Y_obs'].values, self.ds['Y_calc'].values, (self.ds['Y_bkg_auto']+self.ds['Y_bkg_gsas']).values
        return X, Yobs, Ycalc, Ybkg

    def plot_refinement_results(self, plt_range=None, clearplt=False):
        if clearplt:
            plt.clf()
        
        # if x-axis range isn't given, just use original q_range
        if plt_range is None:
            plt_range = self.q_range
        # otherwise, ensure the plt_range is valid, and if it isn't, cut to original q_range
        elif plt_range[0] < self.q_range[0]:
            plt_range[0] = self.q_range[0]
        elif plt_range[1] > self.plt_range[1]:
            plt_range[1] = self.q_range[1]

        # set x ticks
        tick_ct = 12.0
        tick_size = (plt_range[1]-plt_range[0])/tick_ct
        x_ticks = np.round(np.arange(plt_range[0],plt_range[1]+tick_size,tick_size),3)

        # create a mosaic and plot based on its pattern
        fig = plt.figure(figsize=(8,10),dpi=96)
        mosaic = """
        RRR
        YYY
        YYY
        PPP
        DDD
        """
        ax_dict = fig.subplot_mosaic(mosaic, sharex=True)

        # plot the cake pattern in subplot 1 from the top
        ax = ax_dict["R"]

        np.log(self.da_i2d).plot.imshow(ax=ax,robust=True,add_colorbar=False,cmap='Greys',vmin=0)
        ax.set_ylabel(self.ds.i2d.ylabel)
        ax.set_xlabel(None)
        ax.set_xlim([plt_range[0],plt_range[1]])
        Rwp, GoF = self.gpx.data['Covariance']['data']['Rvals']['Rwp'], self.gpx.data['Covariance']['data']['Rvals']['GOF']
        ax.set_title('%s\n(R$_{wp}$=%.3f GoF=%.3f, Temp(C)=%s)'%(self.name.split('/')[-1],Rwp,GoF,self.ds.attrs['HAB_temperature']))
        ax.set_facecolor('#FFFED4')

        # plot observed and calculated I vs 2Theta in subplot 2 from the top
        ax = ax_dict["Y"]
        X, Yobs, Ycalc, Ybkg= self.ds['X_in_q'].values, self.ds['Y_obs'].values, self.ds['Y_calc'].values, (self.ds['Y_bkg_auto']+self.ds['Y_bkg_gsas']).values
        ax.plot(X, np.log( Yobs ),label='Y$_{obs}$',lw=2,color='k')
        ax.plot(X, np.log( Ycalc+Ybkg ),label='Y$_{calc}$+Y$_{bkg}$',lw=1,color='y')
        ax.fill_between(X, np.log( Ybkg ),label='Y$_{bkg}$',alpha=0.2)
        ax.set_ylim(bottom=np.log(min(min(Ybkg),min(Ybkg)))-0.1)
        ax.set_ylabel('Log$_{10}$(counts) (a.u.)')
        ax.set_xlim([plt_range[0],plt_range[1]])
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks)

        # find lattice constants for each phase and save them to use as label strings in subplot 2
        phases = functions_lib_version.get_valid_phases(self.gpx)
        phase_ct = len(phases)
        label_strs = [None] * phase_ct
        label_colors = [None] * phase_ct
        for ep, phase in enumerate(phases):
            label_strs[ep] = phase
            consts = iter(['a', 'b', 'c'])
            for value in np.unique(list(functions_lib_version.get_cell_consts(self.gpx, phase).values())):
                label_strs[ep] = label_strs[ep] + "\n(" + next(consts) + " = " + str(round(value,6)) + ")"
            label_colors[ep] = "C%d" % ep

        # create Line2D objects to use as labels for subplot 2 legend
        from matplotlib.lines import Line2D
        custom_handles = []
        for i in np.arange(0,phase_ct,1):
            custom_handles.append(Line2D([0], [0], marker='o',color=label_colors[i], label='Scatter',
                          markerfacecolor=label_colors[i], markersize=5))

        # get set legend handles from the plot and set legend with all handles and labels
        handles, labels = ax.get_legend_handles_labels()
        all_handles = handles + custom_handles
        print(label_strs)
        all_labels = labels + label_strs
        ax.legend(all_handles, all_labels)

        # use gpx_plotter() to plot intensities as points, stems, in subplot 3, subplot 2, respectively
        ax = ax_dict["Y"]
        ax_bottom = ax_dict["P"]
        ax_bottom.sharex(ax_dict["R"])
        functions_lib_version.gpx_plotter(
            self.gpx,
            line_axes=[ax,ax_bottom],
            stem_axes=[ax_bottom],
            radial_range=plt_range,
            phases=functions_lib_version.get_valid_phases(self.gpx),
            marker="o",
            stem=True,
            unit="d",
            plot_unit="q_A^-1",
            y_shift=0.1
        )

        # plot difference in observed vs calculated diffraction pattern in subplot 4
        ax = ax_dict["D"]
        ax.plot(X, Yobs-Ycalc-Ybkg,label='diff.',lw=1,color='r')
        ax.set_xlabel('Scattering vector $q$ ($\AA^{-1}$)')
        ax.set_ylabel('counts (a.u.)')
        ax.set_xlim([plt_range[0],plt_range[1]])
        ax.legend(loc='best')
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks)

        # ensuring 2nd and 4th plots only have x-axis ticks
        for ax in ["Y", "D"]:
            ax_dict[ax].xaxis.set_tick_params(which='both', labelbottom=True)

    def plot_corr_matrix(self, clearplt=False):
        """
        Plot the correlation matrix of the most recently refined variables as a heatmap. It is best to run this method after each refinement, to see highly correlated variables, and try to reduce their high correlations, for a more honest fit.

        Args:
            clearplt (bool, optional): toggle for clearing the plots; useful if using something like %matplotlib widget in jupyter nb, and plots are being overwritted by other plots, which is an occasional glitch. Note, when true this will allow just one plot to be shown in the document at once, but with no overwritting. Defaults to False.
        """
        import seaborn as sns
        
        covmat = np.array(self.gpx['Covariance']['data']['covMatrix'])

        sig = np.array(self.gpx['Covariance']['data']['sig'])
        cor_denom =  np.outer(sig, sig)

        labels = self.gpx['Covariance']['data']['varyList']
        corr_matrix = np.multiply(covmat, 1.0/cor_denom)
        
        for i,row in enumerate(corr_matrix):
            for j,elem in enumerate(row):
                if elem > 1 or elem < -1:
                    cor_denom[i][j]=0
                    
        # Clear any existing plots
        if clearplt:
            plt.clf()

        # Create the heatmap
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(corr_matrix, annot=False, fmt=".2f", xticklabels=labels, yticklabels=labels, cmap='viridis', linewidths=1, linecolor='black')


        plt.title('Correlation Matrix Heatmap')
        plt.show()
            
    def tmp_delete(self, mode='file'):
        """
        Deletes files generated by refinement stored in .tmp folder, or the .tmp folder itself, which is stored in the .gsas2_scratch folder

        Args:
            mode (str, optional): delete either .tmp folder ('tree') or files within it ('file'). Defaults to 'file'.
        """
        # modes: 'file' deletes all contents of gsas2_scratch except instprm; 'tree' deletes entire gsas2_scratch folder. default mode='file'
        if mode=='tree' and os.path.isdir(self.gsas2_scratch):
           shutil.rmtree(self.gsas2_scratch, ignore_errors=True)
        elif mode=='file':
            for file in os.listdir(self.gsas2_scratch):
                if file != 'gsas.instprm':
                    if os.path.isfile('%s/%s'%(self.gsas2_scratch,file)):
                        os.remove('%s/%s'%(self.gsas2_scratch,file))
    
    def instprm_updater(self):
        """
        Writes new instprm file with refined instprm's at instrument_parameters_file directory, name format 'old_instprm_name.new'
        """
        instprm_dict = self.gpx['PWDR data.xy']['Instrument Parameters'][0]

        with open(self.instrument_parameters_file+'.new', 'w') as f:
            f.write('#GSAS-II instrument parameter file; do not add/delete items!\n')
            f.write('Type:PXC\n')
            f.write('Bank:1.0\n')
            f.write('Lam:%s\n'%(instprm_dict['Lam'][1]))
            f.write('Polariz.:%s\n'%(instprm_dict['Polariz.'][1]))
            f.write('Azimuth:%s\n'%(instprm_dict['Azimuth'][1]))
            f.write('Zero:%s\n'%(instprm_dict['Zero'][1]))
            f.write('U:%s\n'%(instprm_dict['U'][1]))
            f.write('V:%s\n'%(instprm_dict['V'][1]))
            f.write('W:%s\n'%(instprm_dict['W'][1]))
            f.write('X:%s\n'%(instprm_dict['X'][1]))
            f.write('Y:%s\n'%(instprm_dict['Y'][1]))
            f.write('Z:%s\n'%(instprm_dict['Z'][1]))
            f.write('SH/L:%s\n'%(instprm_dict['SH/L'][1]))
    
    def set_new_data(self, nc_path):
        """
        This method should only be used if you want to continue with refinement using the parameters of the last refinement, with a new data file; otherwise, make a new refiner object.

        Args:
            nc_path (_type_): _description_

        Raises:
            ValueError: _description_
        """
        
        try:
            import GSASIIscriptable as G2sc
            import pybaselines
        except:
            raise ModuleNotFoundError("GSAS-II path not found.\n Install GSAS-II, and insert the GSAS-II_path into the refiner class.")
        
        # save previous wavelength, in case this data file has no wavelength
        wavelength = self.ds.attrs["wavelength"]
        
        #update name for plotting
        self.name = nc_path
        
        with xr.open_dataset(nc_path) as ds:
            for k in ['Y_obs','Y_calc','Y_bkg_auto','Y_bkg_gsas','Y_bkg','Y_bkg_arpls','Y_bkg_auto']:
                if k in ds.keys():
                    del ds[k]
                    
            self.da_i2d = ds.i2d.sel(radial=slice(self.q_range[0],self.q_range[1])).astype('float32')
            self.da_i2d_m = self.da_i2d.mean(dim='azimuthal')
            
            arpls_ = pybaselines.Baseline(x_data=self.da_i2d_m.radial.values).arpls((self.da_i2d_m).values, lam=1e5)[0]
            shift_ = min((self.da_i2d_m).values - arpls_)
            bkg_arpls = (arpls_+shift_)
            self.bkg_auto = bkg_arpls
            Y_to_gsas = self.da_i2d_m.values-self.bkg_auto

            self.y_scale = 1000/max(Y_to_gsas)
            self.y_baseline = 10

            # set gsas intensity profile and save as xy file in gsas2_scratch directory
            Y_to_gsas = self.y_baseline+Y_to_gsas*self.y_scale

            # try to find wavelength; if it is not present, reset it
            try:
                ds.attrs['wavelength']
            except:
                if self.verbose:
                    print("Note, wavelength not found; using previous wavelength")
                ds.attrs["wavelength"] = wavelength
            
            X_to_gsas = np.rad2deg(functions_lib_version.q_to_twotheta(self.da_i2d_m.radial.values, wavelength=(ds.attrs['wavelength']*1.0e10)))
            np.savetxt('%s/data.xy'%self.gsas2_scratch, np.column_stack( (X_to_gsas,Y_to_gsas) ), fmt='%.4f %.4f')

            # save ds as an instance variable
            self.ds = ds

            # save gpx as an instance variable
            #self.gpx = G2sc.G2Project(newgpx='%s/gsas.gpx'%self.gsas2_scratch)
            #self.gpx.data['Controls']['data']['max cyc'] = 100
            #self.gpx.add_powder_histogram('%s/data.xy'%self.gsas2_scratch,self.instrument_parameters_file)

            #hist   = self.gpx.histograms()[0]
            #for p in self.phases:
            #    self.gpx.add_phase(p['cif_abs_path'],phasename=p['phase_name'], histograms=[hist],fmthint='CIF')
            self.gpx.save()
    
    def get_q_range(self):
        """
        Returns:
            list: [q minimum, q maximum]
        """
        return self.q_range

    def print_wR(self, header=''):
        """
        Print the name of the PWDR data histogram associated to this refiner along with its overall weighted profile R factor.
        
        Args:
            header (str, optional): adds 'header' before wR. Defaults to ''.
        """
        clear_output()
        hist = self.gpx.histograms()[0]
        print('\n'+header)
        print("\t{:20s}: {:.2f}".format(hist.name,hist.get_wR()))
        print("")

    def set_LeBail(self, LeBail=False):
        """
        Toggle LeBail analysis for refinements

        Args:
            LeBail (bool, optional): if True LeBail analysis is performed. Defaults to False.
        """
        self.gpx.set_refinement({"set":{'LeBail': LeBail}})
        self.gpx.save()
        if(self.verbose):
            print('\nLeBail is set to %s\n '%(LeBail))

    def refine_background(self,num_coeffs):
        """
        Refines background. Increased coeffs generates tighter fit. However, overfitting by adding too many variables fits the
        curve to noise and unwanted information.
        
        Args:
            num_coeffs (int): number of coefficients used in refinement
        """
        try:
            rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        except:
            rwp_old = 'na'

        ParDict = {'set': {'Background': {'refine': True,
                                        'type': 'chebyschev-1',
                                        'no. coeffs': num_coeffs
                                        }}}
        self.gpx.set_refinement(ParDict)
        if self.verbose:
            self.gpx.refine(); rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        else:
            with HiddenPrints():
                self.gpx.refine(); rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']

        self.gpx.set_refinement({'set': {'Background': {'refine': False}}})
        self.gpx.save()

        if self.verbose:
            try:
                print('\n\n\nBackground is refined: Rwp=%.3f (was %.3f)\n\n\n '%(rwp_new,rwp_old))
            except:
                print('\n\n\nBackground is refined: Rwp=%.3f \n\n\n '%(rwp_new))


    def refine_cell_params(self,phase_ind=None):
        """
        Refines cell parameters.

        Args:
            phase_ind (index, optional): index of phase to be refined. If None, refines all. Defaults to None.
        """
        rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        phases = self.gpx.phases()
        for e,p in enumerate(phases):
            if phase_ind is None:
                self.gpx['Phases'][p.name]['General']['Cell'][0]= True
            else:
                if e == phase_ind:
                    self.gpx['Phases'][p.name]['General']['Cell'][0]= True
        
        self.gpx.save()
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()

        phases = self.gpx.phases()
        for p in phases:
            self.gpx['Phases'][p.name]['General']['Cell'][0]= False
            
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        self.gpx.save()

        if self.verbose:
            if phase_ind is None:
                print('\n\n\nCell parameters of all phases are refined simultaneously: Rwp=%.3f (was %.3f)\n\n\n '%(rwp_new,rwp_old))
            else:
                print('\n\n\nCell parameters of phase #d is refined: Rwp=%.3f (was %.3f)\n\n\n '%(phase_ind,rwp_new,rwp_old))


    def refine_strain_broadening(self):
        """
        Refines strain broadening
        """
        rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        ParDict = {'set': {'Mustrain': {'type': 'isotropic',
                                        'refine': True
                                        }}}
        self.gpx.set_refinement(ParDict)
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()

        ParDict = {'set': {'Mustrain': {'type': 'isotropic',
                                'refine': False
                                }}}
        self.gpx.set_refinement(ParDict)
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        self.gpx.save()

        if self.verbose:
            print('\n\n\nStrain broadeing is refined: Rwp=%.3f (was %.3f)\n\n\n '%(rwp_new,rwp_old))


    def refine_size_broadening(self):
        """
        Refines size broadening
        """
        rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        ParDict = {'set': {'Size': {'type': 'isotropic',
                                        'refine': True
                                        }}}
        self.gpx.set_refinement(ParDict)
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()

        ParDict = {'set': {'Size': {'type': 'isotropic',
                                'refine': False
                                }}}
        self.gpx.set_refinement(ParDict)
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        self.gpx.save()

        if self.verbose:
            print('\nSize broadening is refined: Rwp=%.3f (was %.3f)\n '%(rwp_new,rwp_old))


    def refine_inst_parameters(self,inst_pars_to_refine=['X', 'Y', 'Z', 'Zero', 'SH/L', 'U', 'V', 'W']):
        """
        Refines instrument paramters

        Args:
            inst_pars_to_refine (list, optional): which instrument parameters to refine. Defaults to ['X', 'Y', 'Z', 'Zero', 'SH/L', 'U', 'V', 'W'].
        """
        rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        self.gpx.set_refinement({"set": {'Instrument Parameters': inst_pars_to_refine}})
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()

        ParDict = {"clear": {'Instrument Parameters': ['X', 'Y', 'Z', 'Zero', 'SH/L', 'U', 'V', 'W']}}
        self.gpx.set_refinement(ParDict)
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        self.gpx.save()

        if self.verbose:
            print('\nInstrument parameters %s are refined: Rwp=%.3f (was %.3f)\n '%(inst_pars_to_refine,rwp_new,rwp_old))
            
    def reset_verbose(self, verbose=True):
        """
        Reset the self.verbose variable for future uses.

        Args:
            verbose (boolean): value to set self.verbose to. Defaults to True.
        """
        self.verbose = verbose        
    
class seqRefiner:

    def __init__(self,
                 nc_paths,
                 phases,
                 instrument_parameters_file,
                 gsas2_scratch = None,
                 GSASII_path = None,
                 da_input_bkg = None,
                 q_range = None,
                 verbose = True
                 ):
        """
        This is the sequential refiner object. It allows you to perform refinement, just like with the refiner object, except sequentially. Initialization is identical, except instead of a single nc_path, we now have an array nc_paths. The rest of the parameters will be described below. 
        
        Upon instantiation, the nc files are read and their data is extracted. Simple background fitting is conducted for each one. The ds variable is no longer a single xarray, as it was in Refiner class, but instead is an array of xarrays, with the ith index corresponding to the ith file at nc_paths[i]. 

        Args:
            nc_paths (array_like): the paths to your measurment files, the '.nc' files
            phases (array_like): an array of phases, each element a dictionary of the form {'cif_abs_path':'LaB6.cif','phase_name':'LaB6','scale':1}
            instrument_parameters_file (str): a string representing the path to the '.instprm' file
            gsas2_scratch (str, optional): See description above for more details, represents path to gsas2_scratch directory. Defaults to None, and if default, will create a GSAS-II directory on its own.
            GSASII_path (str, optional): a string representing the path of the GSAS-II module containing directory, should be of the form '[your path]/GSASII/GSAS-II/GSASII'. Defaults to None, only use default if you have your GSAS-II path set up, in a way that you can directly import GSAS-II modules using 'import [module]'.
            da_input_bkg (xarray.DataArray, optional): xarray DataArray representing the expected background spectrum. Defaults to None, only use default if the background spectrum is unknown; providing an expected background spectrum will yield better refinement.
            q_range (tuple, optional): Range of q values to consider (in units 1/A0, although it depends on your data units). Defaults to (1,10).
            verbose (boolean): Toggle print statements. On by default; set to False for asthetic purposes.
        Returns:
            seqRefiner: seqRefiner object with its associated functions, as seen following the initialization function.
        """

        if GSASII_path is not None:
            G2sc, pybaselines = set_GSAS_II_path(GSASII_path)
        else:
            try:
                import GSASIIscriptable as G2sc
                import pybaselines
            except:
                raise ModuleNotFoundError("GSAS-II path not found.\n Install GSAS-II, and insert the GSAS-II_path into the refiner class.")


        if gsas2_scratch is None:
            if not os.path.exists('.gsas2_scratch'):
                os.makedirs('.gsas2_scratch',exist_ok=True)
            randstr = ''.join(random.choices(string.ascii_uppercase+string.digits, k=7))
            gsas2_scratch = '%s/%.2f_%s.tmp'%('.gsas2_scratch',time.time(),randstr)
            os.makedirs(gsas2_scratch)
            self.gsas2_scratch = gsas2_scratch
        else:
            self.gsas2_scratch = gsas2_scratch


        # set instance variables
        self.nc_paths = nc_paths
        self.names = nc_paths
        self.phases = phases
        self.instrument_parameters_file = instrument_parameters_file
        self.q_range = q_range
        self.da_input_bkg = da_input_bkg
        
        self.count = len(nc_paths)
        self.ds_list = [None] * self.count
        self.bkg_auto_list = [None] * self.count
        
        self.wavelength = None
        self.verbose = verbose
        
        # save gpx as an instance variable
        self.gpx = G2sc.G2Project(newgpx='%s/gsas.gpx'%self.gsas2_scratch)
        for p in phases:
            self.gpx.add_phase(p['cif_abs_path'],phasename=p['phase_name'],fmthint='CIF')
        #self.gpx.data['Controls']['data']['max cyc'] = 100

        for i, nc_path in enumerate(self.nc_paths):
            if verbose:
                print("Saving file", i, "to GPX")
            
            with xr.open_dataset(nc_path) as ds:
                if self.wavelength is None:
                    try:
                        self.wavelength = ds.attrs["wavelength"]
                    except:
                        print("No wavelength present in the first .nc file.")
                
                #auto set q_range
                if self.q_range is None:
                    self.q_range = [float(ds.radial.min()), float(ds.radial.max())]
                    
                for k in ['Y_obs','Y_calc','Y_bkg_auto','Y_bkg_gsas','Y_bkg','Y_bkg_arpls','Y_bkg_auto']:
                    if k in ds.keys():
                        del ds[k]

                da_i2d = ds.i2d.sel(radial=slice(self.q_range[0],self.q_range[1])).astype('float32')
                da_i2d_m = da_i2d.mean(dim='azimuthal')

                # set background function
                if self.da_input_bkg is None:
                    arpls_ = pybaselines.Baseline(x_data=da_i2d_m.radial.values).arpls(da_i2d_m.values, lam=1e5)[0]
                    shift_ = min((da_i2d_m).values - arpls_)
                    bkg_arpls = (arpls_+shift_)
                    self.bkg_auto_list[i] = bkg_arpls

                else:
                    self.bkg_auto_list[i] = self.da_input_bkg.sel(radial=slice(self.q_range[0],self.q_range[1]))
                    blank_scale = (da_i2d_m[0] / self.bkg_auto_list[i]).values
                    while (min((da_i2d_m.values-blank_scale*self.bkg_auto_list[i].values)) < 0):
                        blank_scale = blank_scale*0.95
                    da_input_bkg_scaled = blank_scale*self.bkg_auto_list[i]
                    
                    # sometimes the x_data values (the radial values) of the data and the background disagree; if each disagreement is less than 10^-6, we will ignore the dissagreement, and will treat replace the radial values of the background with those of the data(the choise to replace those of the background is arbitrary)
                    try:
                        arpls_ = pybaselines.Baseline(x_data=da_i2d_m.radial.values).arpls((da_i2d_m-da_input_bkg_scaled).values, lam=1e5)[0]
                    except:
                        diffs = da_i2d_m.radial.data != da_input_bkg_scaled.radial.data
                        diff_data = da_i2d_m.radial.data[diffs]
                        diff_bkg = da_input_bkg_scaled.radial.data[diffs]
                        if abs(max(diff_data - diff_bkg)) < 0.000001:
                            da_input_bkg_scaled = da_input_bkg_scaled.assign_coords(radial=da_i2d_m.radial)
                            arpls_ = pybaselines.Baseline(x_data=da_i2d_m.radial.values).arpls((da_i2d_m-da_input_bkg_scaled).values, lam=1e5)[0]
                        else:
                            raise ValueError("The background's radial values are two different from the observed data's radial values (there is at least 1 point with a difference above 10e-6)")
                    shift_ = min((da_i2d_m-da_input_bkg_scaled).values - arpls_)
                    bkg_arpls = (arpls_+shift_)
                    self.bkg_auto_list[i] = bkg_arpls + da_input_bkg_scaled.values
                
                Y_to_gsas = da_i2d_m.values-self.bkg_auto_list[i]

                self.y_scale = 1000/max(Y_to_gsas)
                self.y_baseline = 10

                # set gsas intensity profile and save as xy file in gsas2_scratch directory
                Y_to_gsas = self.y_baseline+Y_to_gsas*self.y_scale
                X_to_gsas = np.rad2deg(functions_lib_version.q_to_twotheta(da_i2d_m.radial.values, wavelength=(self.wavelength*1.0e10)))
                np.savetxt('%s/data.xy'%self.gsas2_scratch, np.column_stack( (X_to_gsas,Y_to_gsas) ), fmt='%.4f %.4f')

                # save ds in the array of all the datasets
                self.ds_list[i] = ds
                self.gpx.add_powder_histogram('%s/data.xy'%self.gsas2_scratch,instrument_parameters_file, phases='all')

                hist = self.gpx.histograms()[-1]
                self.gpx.set_Frozen(variable = None, histogram = hist)
                self.gpx.save()
        self.gpx.set_Controls('sequential',self.gpx.histograms())
        self.gpx.set_Controls('cycles',100)
        self.gpx.set_Controls('seqCopy',True)
        self.gpx.save()

    def save_refinement(self):
        """
        This method is to be run after the refinements are conducted (ex, after myrefiner.refine_background(5)), and it saves all the refinement results to the ds files. The call of this method is necessary for the plotting methods to work. The ds_list object which is returned may be useful as it stores the histograms for every individual nc_file, and is much easier to work with than the gpx object; however, the gpx object stores all the same information as the elements in the ds_list object.
        
        Returns:
            array_like: a tuple (length 2 array), with the gpx G2Project object and the ds_list array [self.gpx, self.ds_list].
        """
        

        for i in range(self.count):
            # get x, yobs, ycalc, and background from gpx histograms
            histogram = self.gpx.histograms()[i]

            gsas_X  = histogram.getdata('x').astype('float32')
            gsas_Yo = histogram.getdata('yobs').astype('float32')
            gsas_Yc = histogram.getdata('ycalc').astype('float32')
            gsas_B  = histogram.getdata('Background').astype('float32')

            # now convert Y back
            Y_Obs_from_gsas  = ((gsas_Yo -self.y_baseline  )/self.y_scale )
            Y_Bkg_from_gsas  = ((gsas_B  -self.y_baseline  )/self.y_scale )
            Y_calc_from_gsas = ((gsas_Yc -self.y_baseline  )/self.y_scale ) # this also includes background from GSAS

            da_i2d = self.ds_list[i].i2d.sel(radial=slice(self.q_range[0],self.q_range[1])).astype('float32')
            da_i2d_m = da_i2d.mean(dim='azimuthal')

            # add data to the variable ds
            self.ds_list[i] = self.ds_list[i].assign_coords(
                {"X_in_q": da_i2d_m.radial.values.astype('float32')},
                )
            self.ds_list[i] = self.ds_list[i].assign_coords(
                {"X_in_tth": gsas_X.astype('float32')},
                )
            self.ds_list[i] = self.ds_list[i].assign_coords(
                {"X_in_d": functions_lib_version.q_to_d(da_i2d_m.radial.values.astype('float32'))},
                )

            self.ds_list[i]['Y_obs'] = xr.DataArray(
                                        data=da_i2d_m.values,
                                        dims=['X_in_q'],
                                    )
            self.ds_list[i]['Y_calc'] = xr.DataArray(
                                        data=Y_calc_from_gsas-Y_Bkg_from_gsas,
                                        dims=['X_in_q'],
                                    )# by Ycalc, we mean only calculated peaks, no background

            self.ds_list[i]['Y_bkg_gsas'] = xr.DataArray(
                                        data=Y_Bkg_from_gsas,
                                        dims=['X_in_q'],
                                    )
            self.ds_list[i]['Y_bkg_auto'] = xr.DataArray(
                                        data=self.bkg_auto_list[i],
                                        dims=['X_in_q'],
                                    )
            
            file_path = os.path.join(self.gsas2_scratch, 'gsas.lst')

            # Check if the directory and file exist
            if not os.path.exists(self.gsas2_scratch):
                print(f"Directory does not exist: {self.gsas2_scratch}")
            elif self.verbose:
                print(f"Directory exists: {self.gsas2_scratch}")

            if not os.path.isfile(file_path):
                print(f"File not found: {file_path}")
            else:
                # Open and read the file if it exists
                with open(file_path) as lst_file:
                    gsas_lst = lst_file.read()
                    self.ds_list[i].attrs['gsas_lst'] = gsas_lst
        
        #with open('%s/gsas.lst'%(self.gsas2_scratch)) as lst_file:
        #    gsas_lst = lst_file.read()
        #    self.ds.attrs['gsas_lst'] = gsas_lst

        return [self.gpx,self.ds_list]

    def get_histograms(self):
        """
        Returns an array of the histograms (G2PwdrData objects) associated with each data file. A wrapper on the G2Project.histogram() method, so that you don't need to access gpx.

        Returns:
            self.gpx.histograms() (G2PwdrData array): array storing the histograms of the gpx object.
        """
        return self.gpx.histograms()

    def update_ds(self):
        """
        Updates the ds file in nc_path
        """
        self.ds.to_netcdf(self.nc_path+'.new.nc',engine="scipy")
        time.sleep(0.1)
        shutil.move(self.nc_path+'.new.nc',self.nc_path)
        
    def get_ds(self):
        """
        Returns:
            self.ds_list (nCDF list): list of nCDF objects storing the experimental data and refinement fits
        """
        return self.ds_list

    def update_gpx(self):
        """
        Updates the gpx file in nc_path
        """
        shutil.copyfile(self.gpx['Controls']['data']['LastSavedAs'],'%s.gpx'%self.nc_path[:-3])

    def get_gpx(self):
        """
        Returns:
            self.gpx (G2Project): gpx project object storing the relevant variables for refinement
        """
        return self.gpx
    
    def get_q_range(self):
        """
        Returns:
            self.q_range (int array): tuple (array of length 2) storing (min, max) for the refinement range
        """
        
        return self.q_range

    def get_i2d(self, i):
        """
        A method for getting the cake pattern (reorganized detector image where the radial direction is on the xaxis), and the azimuthal mean intensity vs 2theta profile, from a given nc_file at index i.

        Args:
            i (int): index for the corresponding element in files to extract data from

        Returns:
            xarray: 2d cake pattern xarray
            xarray: 1d intensity vs 2theta profile
        """
        da_i2d = self.ds_list[i].i2d.sel(radial=slice(self.q_range[0],self.q_range[1])).astype('float32')
        da_i2d_m = da_i2d.mean(dim='azimuthal')
        
        return da_i2d, da_i2d_m

    def get_y_vals(self, i):
        """
        Returns:
            X, Yobs, Ycalc, Ybkg: returns four xarray DataArrays:
            0. X (q values of data point for a given index)
            1. Yobs (observed intensity of data point for a corresponding index in the X array)
            2. Ycalc (calculated intensity attributed to the bragg diffraction peaks of a data point given corresponding index in the X array)
            3. Ybkg (calculated intensity attributed to the background of a data point given corresponding index in the X array; note that to get the full intensity profile, fitted to Yobs, you must sum Ycalc and Ybkg)
        """
        X, Yobs, Ycalc, Ybkg= self.ds_list[i]['X_in_q'].values, self.ds_list[i]['Y_obs'].values, self.ds_list[i]['Y_calc'].values, (self.ds_list[i]['Y_bkg_auto']+self.ds_list[i]['Y_bkg_gsas']).values
        return X, Yobs, Ycalc, Ybkg
    
    def get_ref_scores(self, i):
        """
        Method for getting the refinement scores for an individual file.

        Args:
            i (int): index for the corresponding element in files to extract data from.

        Returns:
            Rwp, GOF (double, Double): tuple (2 element array) of the Rwp and GOF fit score metrics.
        """
        if not i == 0:
            Rwp, GOF = self.gpx.data['Sequential results']['data']['PWDR data.xy_'+str(i)]['Rvals']['Rwp'], self.gpx.data['Sequential results']['data']['PWDR data.xy_'+str(i)]['Rvals']['GOF']
        else:
            Rwp, GOF = self.gpx.data['Sequential results']['data']['PWDR data.xy']['Rvals']['Rwp'], self.gpx.data['Sequential results']['data']['PWDR data.xy']['Rvals']['GOF']
        return Rwp, GOF

    def plot_one(self, index=0, plt_range=None, clearplt=False):
        """
        Method for plotting just one of the fits, after sequential refinement. Plots the cake pattern (reorganized detector image where the radial direction is on the xaxis) on the top, followed by the observed and calculated intensity vs 2theta plot, followed by theoretical diffraction peaks from the cif file, followed by the difference between the calculated and observed intensities. 
        
        Args:
            plt_range (int array): tuple (array of length 2) storing (min, max) for the plotting range
            index (int): index of the file array to be plotted
            clearplt (boolean): toggle for clearing the plots; useful if using something like %matplotlib widget in jupyter nb, and plots are being overwritted by other plots, which is an occasional glitch. Note, when true this will allow just one plot to be shown in the document at once, but with no overwritting.
        """
        if clearplt:
            plt.clf()

        # if x-axis range isn't given, just use original q_range
        if plt_range is None:
            plt_range = self.q_range
        # otherwise, ensure the plt_range is valid, and if it isn't, cut to original q_range
        elif plt_range[0] < self.q_range[0]:
            plt_range[0] = self.q_range[0]
        elif plt_range[1] > self.plt_range[1]:
            plt_range[1] = self.q_range[1]

        # set x ticks
        tick_ct = 12.0
        tick_size = (plt_range[1]-plt_range[0])/tick_ct
        x_ticks = np.round(np.arange(plt_range[0],plt_range[1]+tick_size,tick_size),3)

        # create a mosaic and plot based on its pattern
        fig = plt.figure(figsize=(8,10),dpi=96)
        mosaic = """
        RRR
        YYY
        YYY
        PPP
        DDD
        """
        ax_dict = fig.subplot_mosaic(mosaic, sharex=True)

        # plot the cake pattern in subplot 1 from the top
        ax = ax_dict["R"]
        
        da_i2d = self.ds_list[index].i2d.sel(radial=slice(self.q_range[0],self.q_range[1])).astype('float32')
        da_i2d_m = da_i2d.mean(dim='azimuthal')

        np.log(da_i2d).plot.imshow(ax=ax,robust=True,add_colorbar=False,cmap='Greys',vmin=0)
        ax.set_ylabel(self.ds_list[index].i2d.ylabel)
        ax.set_xlabel(None)
        ax.set_xlim([plt_range[0],plt_range[1]])
        if not index == 0:
            Rwp, GoF = self.gpx.data['Sequential results']['data']['PWDR data.xy_'+str(index)]['Rvals']['Rwp'], self.gpx.data['Sequential results']['data']['PWDR data.xy_'+str(index)]['Rvals']['GOF']
        else:
            Rwp, GoF = self.gpx.data['Sequential results']['data']['PWDR data.xy']['Rvals']['Rwp'], self.gpx.data['Sequential results']['data']['PWDR data.xy']['Rvals']['GOF']
        ax.set_title('%s\n(R$_{wp}$=%.3f GoF=%.3f, Temp(C)=%s)'%(self.names[index].split('/')[-1],Rwp,GoF, self.ds_list[index].attrs['HAB_temperature']))
        ax.set_facecolor('#FFFED4')

        # plot observed and calculated I vs 2Theta in subplot 2 from the top
        ax = ax_dict["Y"]
        X, Yobs, Ycalc, Ybkg=  self.ds_list[index]['X_in_q'].values,  self.ds_list[index]['Y_obs'].values,  self.ds_list[index]['Y_calc'].values, ( self.ds_list[index]['Y_bkg_auto']+ self.ds_list[index]['Y_bkg_gsas']).values
        ax.plot(X, np.log( Yobs ),label='Y$_{obs}$',lw=2,color='k')
        ax.plot(X, np.log( Ycalc+Ybkg ),label='Y$_{calc}$+Y$_{bkg}$',lw=1,color='y')
        ax.fill_between(X, np.log( Ybkg ),label='Y$_{bkg}$',alpha=0.2)
        ax.set_ylim(bottom=np.log(min(min(Ybkg),min(Ybkg)))-0.1)
        ax.set_ylabel('Log$_{10}$(counts) (a.u.)')
        ax.set_xlim([plt_range[0],plt_range[1]])
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks)

        # find lattice constants for each phase and save them to use as label strings in subplot 2
        phases = functions_lib_version.get_valid_phases(self.gpx)
        phase_ct = len(phases)
        label_strs = [None] * phase_ct
        label_colors = [None] * phase_ct
        for ep, phase in enumerate(phases):
            label_strs[ep] = phase
            consts = iter(['a', 'b', 'c'])
            for value in np.unique(list(functions_lib_version.get_cell_consts(self.gpx, phase).values())):
                label_strs[ep] = label_strs[ep] + "\n(" + next(consts) + " = " + str(round(value,6)) + ")"
            label_colors[ep] = "C%d" % ep

        # create Line2D objects to use as labels for subplot 2 legend
        from matplotlib.lines import Line2D
        custom_handles = []
        for i in np.arange(0,phase_ct,1):
            custom_handles.append(Line2D([0], [0], marker='o',color=label_colors[i], label='Scatter', markerfacecolor=label_colors[i], markersize=5))

        # get set legend handles from the plot and set legend with all handles and labels
        handles, labels = ax.get_legend_handles_labels()
        all_handles = handles + custom_handles
        all_labels = labels + label_strs
        ax.legend(all_handles, all_labels)

        # use gpx_plotter() to plot intensities as points, stems, in subplot 3, subplot 2, respectively
        ax = ax_dict["Y"]
        ax_bottom = ax_dict["P"]
        ax_bottom.sharex(ax_dict["R"])
        functions_lib_version.gpx_plotter(
            self.gpx,
            line_axes=[ax,ax_bottom],
            stem_axes=[ax_bottom],
            radial_range=plt_range,
            phases=functions_lib_version.get_valid_phases(self.gpx),
            marker="o",
            stem=True,
            unit="d",
            plot_unit="q_A^-1",
            y_shift=0.1
        )

        # plot difference in observed vs calculated diffraction pattern in subplot 4
        ax = ax_dict["D"]
        ax.plot(X, Yobs-Ycalc-Ybkg,label='diff.',lw=1,color='r')
        ax.set_xlabel('Scattering vector $q$ ($\AA^{-1}$)')
        ax.set_ylabel('counts (a.u.)')
        ax.set_xlim([plt_range[0],plt_range[1]])
        ax.legend(loc='best')
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks)

        # ensuring 2nd and 4th plots only have x-axis ticks
        for ax in ["Y", "D"]:
            ax_dict[ax].xaxis.set_tick_params(which='both', labelbottom=True)
            
        plt.show()
            
    def plot_all(self, file_range = None, temp_read = True, plt_range=None, clearplt=False):
        """
        Creates a gif of plots identical to the plot_one() method. Plots the gif using pyplot but also saves the gif as sequential_ref[random int between 0 and 9999].gif.
        
        Args:
            file_range (int array): range of files to select from for the gif; leave as None to select all files
            temp_read (Boolean): toggle for reading and displaying the temperature on a plot
            plt_range (int array): tuple (2d array) representing (min, max) for the x axis of the plot
            clearplt (Boolean): toggle for clearing the plots; useful if using something like %matplotlib widget in jupyter nb, and plots are being overwritted by other plots, which is an occasional glitch. Note, when true this will allow just one plot to be shown in the document at once, but with no overwritting.

        Returns:
            _type_: _description_
        """
        if clearplt:
            plt.clf()
            
        # plot 0th frame first, before creating gif
        if file_range is None:
            filenames = self.nc_paths
        elif file_range[0] >= 0 and file_range[1] <= self.count:
            filenames = self.nc_paths[file_range[0]:file_range[1]]
        elif self.verbose:
            print("Invalid file range, will default to first 10 files to save resources. Set to `None` for all files.")
            filenames = self.nc_paths[0:11]
        
        # if x-axis range isn't given, just use original q_range
        if plt_range is None:
            plt_range = self.q_range
        # otherwise, ensure the plt_range is valid, and if it isn't, cut to original q_range
        elif plt_range[0] < self.q_range[0]:
            plt_range[0] = self.q_range[0]
        elif plt_range[1] > self.plt_range[1]:
            plt_range[1] = self.q_range[1]

        # set x ticks
        tick_ct = 12.0
        tick_size = (plt_range[1]-plt_range[0])/tick_ct
        x_ticks = np.round(np.arange(plt_range[0],plt_range[1]+tick_size,tick_size),3)

        # create a mosaic and plot based on its pattern
        fig = plt.figure(figsize=(8,10),dpi=96)
        mosaic = """
        RRR
        YYY
        YYY
        PPP
        DDD
        """
        ax_dict = fig.subplot_mosaic(mosaic, sharex=True)

        # plot the cake pattern in subplot 1 from the top
        ax = ax_dict["R"]
        
        da_i2d = self.ds_list[0].i2d.sel(radial=slice(self.q_range[0],self.q_range[1])).astype('float32')

        i2dim = np.log(da_i2d).plot.imshow(ax=ax,robust=True,add_colorbar=False,cmap='Greys',vmin=0)
        ax.set_ylabel(self.ds_list[0].i2d.ylabel)
        ax.set_xlabel(None)
        ax.set_xlim([plt_range[0],plt_range[1]])
        
        Rwp, GoF = self.gpx.data['Sequential results']['data']['PWDR data.xy']['Rvals']['Rwp'], self.gpx.data['Sequential results']['data']['PWDR data.xy']['Rvals']['GOF']
        ax.set_title('%s\n(R$_{wp}$=%.3f GoF=%.3f, Temp(C)=%s, FrameNo=%s)'%(filenames[0].split('/')[-1],Rwp,GoF, self.ds_list[0].attrs['HAB_temperature'],0))
        ax.set_facecolor('#FFFED4')

        # plot observed and calculated I vs 2Theta in subplot 2 from the top
        ax = ax_dict["Y"]
        X, Yobs, Ycalc, Ybkg=  self.ds_list[0]['X_in_q'].values,  self.ds_list[0]['Y_obs'].values,  self.ds_list[0]['Y_calc'].values, ( self.ds_list[0]['Y_bkg_auto']+ self.ds_list[0]['Y_bkg_gsas']).values
        yobs, = ax.plot(X, np.log( Yobs ),label='Y$_{obs}$',lw=2,color='k')
        yfit, = ax.plot(X, np.log( Ycalc+Ybkg ),label='Y$_{calc}$+Y$_{bkg}$',lw=1,color='y')
        ax.fill_between(X, np.log( Ybkg ),label='Y$_{bkg}$',alpha=0.2)
        ax.set_ylim(bottom=np.log(min(min(Ybkg),min(Ybkg)))-0.1)
        ax.set_ylabel('Log$_{10}$(counts) (a.u.)')
        ax.set_xlim([plt_range[0],plt_range[1]])
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks)

        # find lattice constants for each phase and save them to use as label strings in subplot 2
        phases = functions_lib_version.get_valid_phases(self.gpx)
        phase_ct = len(phases)
        label_strs = [None] * phase_ct
        label_colors = [None] * phase_ct
        for ep, phase in enumerate(phases):
            label_strs[ep] = phase
            consts = iter(['a', 'b', 'c'])
            for value in np.unique(list(functions_lib_version.get_cell_consts(self.gpx, phase).values())):
                label_strs[ep] = label_strs[ep] + "\n(" + next(consts) + " = " + str(round(value,6)) + ")"
            label_colors[ep] = "C%d" % ep

        # create Line2D objects to use as labels for subplot 2 legend
        from matplotlib.lines import Line2D
        custom_handles = []
        for i in np.arange(0,phase_ct,1):
            custom_handles.append(Line2D([0], [0], marker='o',color=label_colors[i], label='Scatter',
                          markerfacecolor=label_colors[i], markersize=5))

        # get set legend handles from the plot and set legend with all handles and labels
        handles, labels = ax.get_legend_handles_labels()
        all_handles = handles + custom_handles
        all_labels = labels + label_strs
        ax.legend(all_handles, all_labels)

        # use gpx_plotter() to plot intensities as points, stems, in subplot 3, subplot 2, respectively
        ax = ax_dict["Y"]
        ax_bottom = ax_dict["P"]
        ax_bottom.sharex(ax_dict["R"])
        functions_lib_version.gpx_plotter(
            self.gpx,
            line_axes=[ax,ax_bottom],
            stem_axes=[ax_bottom],
            radial_range=plt_range,
            phases=functions_lib_version.get_valid_phases(self.gpx),
            marker="o",
            stem=True,
            unit="d",
            plot_unit="q_A^-1",
            y_shift=0.1
        )

        # plot difference in observed vs calculated diffraction pattern in subplot 4
        ax = ax_dict["D"]
        diff, = ax.plot(X, Yobs-Ycalc-Ybkg,label='diff.',lw=1,color='r')
        ax.set_xlabel('Scattering vector $q$ ($\AA^{-1}$)')
        ax.set_ylabel('counts (a.u.)')
        ax.set_xlim([plt_range[0],plt_range[1]])
        ax.legend(loc='best')
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks)

        # ensuring 2nd and 4th plots only have x-axis ticks
        for ax in ["Y", "D"]:
            ax_dict[ax].xaxis.set_tick_params(which='both', labelbottom=True)
        
        # now create the gif
        
        # Update function for the animation
        def update(frame):
            # assigning variables to be used below            
            if not frame == 0:
                Rwp, GoF = self.gpx.data['Sequential results']['data']['PWDR data.xy_'+str(frame)]['Rvals']['Rwp'], self.gpx.data['Sequential results']['data']['PWDR data.xy_'+str(frame)]['Rvals']['GOF']
            else:
                Rwp, GoF = self.gpx.data['Sequential results']['data']['PWDR data.xy']['Rvals']['Rwp'], self.gpx.data['Sequential results']['data']['PWDR data.xy']['Rvals']['GOF']
            nc_path = filenames[frame]
            X, Yobs, Ycalc, Ybkg= self.get_y_vals(frame)
            difval = Yobs-Ycalc-Ybkg
            
            ax = ax_dict["R"]
            if temp_read:
                try:
                    ax.set_title('%s\n(R$_{wp}$=%.3f GoF=%.3f, Temp(C)=%s, FrameNo=%s)'%(nc_path.split('/')[-1],Rwp,GoF, self.get_ds()[frame].attrs['HAB_temperature'],frame))
                except:
                    raise KeyError('No temperature data in file number:', frame)
            else:
                ax.set_title('%s\n(R$_{wp}$=%.3f GoF=%.3f, FrameNo=%s)'%(nc_path.split('/')[-1],Rwp,GoF,frame))
            ax = ax_dict["Y"]
            ax.set_ylim(bottom=min(np.log(Yobs)))
            yobs.set_ydata(np.log(Yobs) )
            yfit.set_ydata(np.log(Ycalc+Ybkg))
            for collection in ax.collections:
                collection.remove()
            ax.fill_between(X, np.log( Ybkg ),label='Y$_{bkg}$',alpha=0.2, color="c")
            
            ax = ax_dict["D"]
            ax.set_ylim(bottom=-max(np.abs(difval)), top=max(np.abs(difval)))
            diff.set_ydata(difval)
            
            data = np.log(self.ds_list[frame].i2d.sel(radial=slice(self.q_range[0],self.q_range[1])).astype('float32'))
            i2dim.set_array(data)
    
            progress_bar.update(1)  # Update the progress bar

            return [i2dim], yobs, yfit, diff

        frames = len(filenames)
        interval = 100  
        
        from random import randint
        filename = 'sequential_ref' + str(randint(0, 9999)) + '.gif'
        
        from matplotlib.animation import FuncAnimation
        from tqdm import tqdm
        with tqdm(total=frames) as progress_bar:
            ani = FuncAnimation(fig, update, frames=frames, interval=interval)

            # Save the animation as a GIF
            ani.save(filename, writer='pillow')
        if self.verbose:
            print("Gif saved as", filename, "in", os.getcwd())
        plt.show()
    
    def plot_corr_matrix(self, clearplt=False):
        """
        Plot the correlation matrix of the most recently refined variables as a heatmap. It is best to run this method after each refinement, to see highly correlated variables, and try to reduce their high correlations, for a more honest fit.

        Args:
            clearplt (bool, optional): toggle for clearing the plots; useful if using something like %matplotlib widget in jupyter nb, and plots are being overwritted by other plots, which is an occasional glitch. Note, when true this will allow just one plot to be shown in the document at once, but with no overwritting. Defaults to False.
        """
        import seaborn as sns
        
        covmat = np.array(self.gpx['Covariance']['data']['covMatrix'])

        sig = np.array(self.gpx['Covariance']['data']['sig'])
        cor_denom =  np.outer(sig, sig)

        labels = self.gpx['Covariance']['data']['varyList']
        corr_matrix = np.multiply(covmat, 1.0/cor_denom)
        
        for i,row in enumerate(corr_matrix):
            for j,elem in enumerate(row):
                if elem > 1 or elem < -1:
                    cor_denom[i][j]=0
                    
        # Clear any existing plots
        if clearplt:
            plt.clf()

        # Create the heatmap
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(corr_matrix, annot=False, fmt=".2f", xticklabels=labels, yticklabels=labels, cmap='viridis', linewidths=1, linecolor='black')


        plt.title('Correlation Matrix Heatmap')
        plt.show()
    
    def tmp_delete(self, mode='tree'):
        """
        Method for deleting from the gsas2_scratch directory. 
        modes: 
        1. 'file' deletes all contents of gsas2_scratch except instprm; 
        2. 'tree' deletes entire gsas2_scratch folder.

        Args:
            mode (str, optional): deletion mode, as explained above. Defaults to 'Tree'.
        """
        if mode=='tree' and os.path.isdir(self.gsas2_scratch):
           shutil.rmtree(self.gsas2_scratch, ignore_errors=True)
        elif mode=='file':
            for file in os.listdir(self.gsas2_scratch):
                if file != 'gsas.instprm':
                    if os.path.isfile('%s/%s'%(self.gsas2_scratch,file)):
                        os.remove('%s/%s'%(self.gsas2_scratch,file))
    
    def instprm_updater(self):
        """
        Creates and writes to a new instrument parameters file which contains the refined instrument parameters corresponding to every single input file. The file name is the same as the initial instrument parameter file, with '.new' added onto the end.
        """

        with open(self.instrument_parameters_file+'.new', 'w') as f:
            for i in range(self.count):
                if not i == 0:
                    instprm_dict = self.gpx['PWDR data.xy']['Instrument Parameters'][0]
                else:
                    instprm_dict = self.gpx['PWDR data.xy_'+str(i)]['Instrument Parameters'][0]
                    
                f.write('#GSAS-II instrument parameter file; do not add/delete items!\n')
                f.write('Type:PXC\n')
                f.write('Bank:1.0\n')
                f.write('Lam:%s\n'%(instprm_dict['Lam'][1]))
                f.write('Polariz.:%s\n'%(instprm_dict['Polariz.'][1]))
                f.write('Azimuth:%s\n'%(instprm_dict['Azimuth'][1]))
                f.write('Zero:%s\n'%(instprm_dict['Zero'][1]))
                f.write('U:%s\n'%(instprm_dict['U'][1]))
                f.write('V:%s\n'%(instprm_dict['V'][1]))
                f.write('W:%s\n'%(instprm_dict['W'][1]))
                f.write('X:%s\n'%(instprm_dict['X'][1]))
                f.write('Y:%s\n'%(instprm_dict['Y'][1]))
                f.write('Z:%s\n'%(instprm_dict['Z'][1]))
                f.write('SH/L:%s\n'%(instprm_dict['SH/L'][1]))
                f.write('\n')

    def print_wR(self, index=None, header=''):
        """
        Print the name of the PWDR data histogram associated to this refiner along with its overall weighted profile R factor.
        
        Args:
            index (int): specific index to select, if None, just prints all.
            header (str, optional): adds 'header' before wR. Defaults to ''.
        """
        clear_output()
        if index is None:
            for hist in self.gpx.histograms():
                print('\n'+header)
                print("\t{:20s}: {:.2f}".format(hist.name,hist.get_wR()))
                print("")
        else:
            hist = self.gpx.histograms()[0]
            print('\n'+header)
            print("\t{:20s}: {:.2f}".format(hist.name,hist.get_wR()))
            print("")

    def set_LeBail(self, LeBail=True):
        """
        Enable LeBail analysis

        Args:
            LeBail (bool, optional): if True LeBail analysis is performed. Defaults to False.
        """
        for p in self.gpx.phases():
            p.set_refinements({"LeBail": LeBail})
        if self.verbose:
            print('\nLeBail is set to %s\n '%(LeBail))
            if self.verbose:
                self.gpx.refine()
            else:
                with HiddenPrints():
                    self.gpx.refine()
        else:
            with HiddenPrints():
                if self.verbose:
                    self.gpx.refine()
                else:
                    with HiddenPrints():
                        self.gpx.refine()
        

    def refine_background(self, num_coeffs, set_to_false=True):
        """
        Refines background. Increased coeffs generates closer fit.
        
        Args:
            num_coeffs (int): number of coefficients used in background refinement.
            set_to_false (bool, optional): When true, refines onces then sets background refinement to false for future refinements. When false, leaves background refinement on true for future refinements. Defaults to True.
        """
        try:
            rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        except:
            rwp_old = 'na'
            
        for hist in self.gpx.histograms():
            hist.set_refinements({
                'Background': {'type': 'chebyschev-1', 'refine': True, 'no. coeffs': num_coeffs}
            })
        
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()
        
        if(set_to_false):
            for hist in self.gpx.histograms():
                hist.clear_refinements({
                    'Background': {'type': 'chebyschev-1', 'refine': False, 'no. coeffs': num_coeffs}
                })
            self.gpx.save()
        
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        if self.verbose:
            try:
                print('\n\n\nBackground is refined: Rwp=%.3f (was %.3f)\n\n\n '%(rwp_new,rwp_old))
            except:
                print('\n\n\nBackground is refined: Rwp=%.3f \n\n\n '%(rwp_new))

    def refine_cell_params(self, set_to_false=True):
        """
        Refines cell parameters. This is done differently than for the individual refinement class, due to the GSASII sequential refinement approach.
        
        Args:
            set_to_false (bool, optional): When true, refines onces then sets background refinement to false for future refinements. When false, leaves background refinement on true for future refinements. Defaults to True.
        """
        try:
            rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        except:
            rwp_old = 'na'
            
        for p in self.gpx.phases():
            p.set_HAP_refinements({
                'HStrain': True
            })
        
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()
        
        if set_to_false:
            for p in self.gpx.phases():
                p.set_HAP_refinements({
                    'HStrain': False
                })
            self.gpx.save()
        
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        
        if self.verbose:
            try:
                print('\n\nCell parameters are refined: Rwp=%.3f (was %.3f)\n\n\n '%(rwp_new,rwp_old))
            except:
                print('\n\nCell parameters are refined: Rwp=%.3f \n\n\n '%(rwp_new))

    def refine_strain_broadening(self, type = 'isotropic', set_to_false = True):
        """
        Refines strain broadening.

        Args:
            type (str, optional): The type of strain broadening to select. Defaults to 'isotropic'.
            set_to_false (bool, optional): When true, refines onces then sets background refinement to false for future refinements. When false, leaves background refinement on true for future refinements. Defaults to True.
        """
        
        try:
            rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        except:
            rwp_old = 'na'
            
        for p in self.gpx.phases():
            p.set_HAP_refinements({
                'Mustrain': {'type': type, 'refine': True}
            })
        
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()
        
        if set_to_false:
            for p in self.gpx.phases():
                p.set_HAP_refinements({
                    'Mustrain': {'type': type, 'refine': False}
                })
            self.gpx.save()
        
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        if self.verbose:
            try:
                print('\n\n\nStrain broadening is refined: Rwp=%.3f (was %.3f)\n\n\n '%(rwp_new,rwp_old))
            except:
                print('\n\n\nStrain broadening is refined: Rwp=%.3f \n\n\n '%(rwp_new))


    def refine_size_broadening(self, type = 'isotropic', set_to_false = True):
        """
        Refines size broadening.

        Args:
            type (str, optional): The type of size broadening to select. Defaults to 'isotropic'.
            set_to_false (bool, optional): When true, refines onces then sets background refinement to false for future refinements. When false, leaves background refinement on true for future refinements. Defaults to True.
        """
        try:
            rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        except:
            rwp_old = 'na'
            
        for p in self.gpx.phases():
            p.set_HAP_refinements({
                "Size": {'type': type, 'refine': True},
            })
        
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()
        
        if set_to_false:
            for p in self.gpx.phases():
                p.set_HAP_refinements({
                    "Size": {'type': type, 'refine': False},
                })
            self.gpx.save()
        
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        if self.verbose:
            try:
                print('\n\n\nSize broadening is refined: Rwp=%.3f (was %.3f)\n\n\n '%(rwp_new,rwp_old))
            except:
                print('\n\n\nSize broadening is refined: Rwp=%.3f \n\n\n '%(rwp_new))

    def refine_inst_parameters(self,inst_pars_to_refine=['X', 'Y', 'Z', 'Zero', 'SH/L', 'U', 'V', 'W'], set_to_false = True):
        """
        Refines a given list of instrument parameters. If no list is provided, uses defaults. Parameters must be valid in accordance with the GSASII package.

        Args:
            inst_pars_to_refine (list, optional): The list of instrument parameters to refine. Defaults to ['X', 'Y', 'Z', 'Zero', 'SH/L', 'U', 'V', 'W'].
            set_to_false (bool, optional): When true, refines onces then sets background refinement to false for future refinements. When false, leaves background refinement on true for future refinements. Defaults to True.
        """
        try:
            rwp_old = self.gpx['Covariance']['data']['Rvals']['Rwp']
        except:
            rwp_old = 'na'
        
        for hist in self.gpx.histograms():
            hist.set_refinements({ 
                'Instrument Parameters': inst_pars_to_refine
            })
        
        if self.verbose:
            self.gpx.refine()
        else:
            with HiddenPrints():
                self.gpx.refine()
        
        if set_to_false:
            for hist in self.gpx.histograms():
                hist.clear_refinements({ 
                'Instrument Parameters': inst_pars_to_refine
                })
            self.gpx.save()
        
        rwp_new = self.gpx['Covariance']['data']['Rvals']['Rwp']
        if self.verbose:
            try:
                print('\n\n\nInstance parameters are refined: Rwp=%.3f (was %.3f) '%(rwp_new,rwp_old))
                print('Parameters', inst_pars_to_refine, 'were refined.\n\n\n')
            except:
                print('\n\n\nSize broadening is refined: Rwp=%.3f'%(rwp_new))
                print('Parameters', inst_pars_to_refine, 'were refined.\n\n\n')
    def reset_verbose(self, verbose=True):
        """
        Reset the self.verbose variable for future uses.

        Args:
            verbose (boolean): value to set self.verbose to. Defaults to True.
        """
        self.verbose = verbose
