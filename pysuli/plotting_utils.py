import matplotlib.pyplot as plt




# plotting tools for plotting from cif files

def hkl_plotter(
    line_axes=None,
    stem_axes=None,
    mp_id=None,
    final=False,
    structure=None,
    str_file=None,
    label=None,
    color="C0",
    label_x=0.9,
    label_y=0.9,
    unit="2th_deg",
    radial_range=(1, 16),
    wl=0.77,
    scale=1,
    scale_a=1,
    scale_b=1,
    scale_c=1,
    export_cif_as=None,
    stem=True,
):
    """
    Method for plotting an intensity vs two-theta plot on top of a theoretical 
    hkl theoretical bragg peak stem plot. 

    Args:
        line_axes (Axes, optional): matplotlib.axes object for top plot. Defaults to None.
        stem_axes (Axes, optional): matplotlib.axes object for bottom plot. Defaults to None.
        mp_id (str, optional): material project API id. Defaults to None.
        final (bool, optional): whether to get final structure from mp API. Defaults to False.
        structure (Structure, optional): object containing structure information. Defaults to None.
        str_file (str, optional): CIF or other structure containing file. Defaults to None.
        label (str, optional): plot label. Defaults to None.
        color (str, optional): color for the plot. Defaults to "C0".
        label_x (float, optional): x label for stem plot. Defaults to 0.9.
        label_y (float, optional): y label for stem plot. Defaults to 0.9.           
        unit (str, optional): the unit for x axis. Defaults to "2th_deg".
        radial_range (tuple, optional): range for the x axis. Defaults to (1, 16).
        
        bottom (float, optional): _description_. Defaults to -0.2.
        wl (float, optional): _description_. Defaults to 0.77.
        scale (int, optional): _description_. Defaults to 1.
        scale_a (int, optional): _description_. Defaults to 1.
        scale_b (int, optional): _description_. Defaults to 1.
        scale_c (int, optional): _description_. Defaults to 1.
        export_cif_as (_type_, optional): _description_. Defaults to None.
        stem (bool, optional): _description_. Defaults to True.
        stem_logscale (bool, optional): _description_. Defaults to True.
        da_visible (_type_, optional): _description_. Defaults to None.
    """

    if mp_id is not None:
        from mp_api.client import MPRester

        mpr = MPRester("dHgNQRNYSpuizBPZYYab75iJNMJYCklB")  ###
        structure = mpr.get_structure_by_material_id(mp_id, final=final)[0]
    elif structure is None:
        structure = Structure.from_file(str_file)

    structure.lattice = Lattice.from_parameters(
        a=structure.lattice.abc[0] * scale * scale_a,
        b=structure.lattice.abc[1] * scale * scale_b,
        c=structure.lattice.abc[2] * scale * scale_c,
        alpha=structure.lattice.angles[0],
        beta=structure.lattice.angles[1],
        gamma=structure.lattice.angles[2],
    )

    xrdc = XRDCalculator(wavelength=wl)  ###computes xrd pattern given wavelength , debye scherrer rings, and symmetry precision

    if unit == "q_A^-1":
        ps = xrdc.get_pattern(
            structure,
            scaled=True,
            two_theta_range=np.rad2deg(q_to_twotheta(radial_range, wl)),
        )
        X, Y = twotheta_to_q(np.deg2rad(ps.x), wl), ps.y
    elif unit == "2th_deg":
        ps = xrdc.get_pattern(structure, scaled=True, two_theta_range=radial_range)
        X, Y = ps.x, ps.y
    else:
        ps = xrdc.get_pattern(structure, scaled=True, two_theta_range=radial_range)
        X, Y = ps.x, ps.y

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

    if export_cif_as is not None:
        structure.to(fmt="cif", filename=export_cif_as)
        
def phase_plotter(
        wl,
        line_axes=[],
        stem_axes=[],
        radial_range=(1, 16),
        stem=True,
        y_shift=0.1,
        phases=[],
        unit="2th_deg",
        mp_ids=[],
):
    """
    Method for plotting an array of phases using the hkl_plotter() method.

    Args:
        wl (float): wavelength
        line_axes (Axes, optional): matplotlib.axes object for top plot. Defaults to [].
        stem_axes (Axes, optional): matplotlib.axes object for bottom plot. Defaults to [].
        radial_range (tuple, optional): range for the x axis. Defaults to (1, 16).
        stem (bool, optional): whether to add stem plot. Defaults to True.
        y_shift (float, optional): y shift between phase labels. Defaults to 0.1.
        phases (array_like, optional): list of dictionaries each representing a phase, in the format `[{"cif": '_cifs/CeO2.cif', "label": "CeO$_2$", "scale": 1}]`. Defaults to [].
        unit (str, optional): the unit for x axis. Defaults to "2th_deg".
        mp_ids (array_like, optional): list of mp_id strings, to be passed to the hkl_plotter method; if for a given phase there is a ciff provided, set corresponding mp_id=None in the mp_id array
    """
    for ep, phase in enumerate(phases):
        hkl_plotter(
            line_axes=line_axes,
            stem_axes=stem_axes,
            str_file=phase["cif_abs_path"],
            label=phase["phase_name"],
            scale=phase["scale"],
            marker="o",
            color="C%d" % ep,
            label_x=1.02,
            label_y=ep * y_shift,
            unit=unit,
            radial_range=radial_range,
            bottom=-0.2,
            wl=wl,
            stem=stem,
            stem_logscale=False,
            mp_id=mp_ids[ep]
        )





