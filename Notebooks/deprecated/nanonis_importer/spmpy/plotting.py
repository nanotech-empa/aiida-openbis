import os

import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import numpy as np


def matlab_color(num: int):
    """
    Returns an RGB color code according to standard MatLab colors

    Parameters
    ----------
    num : int
        Value iterating between 0 and 6.

    Returns
    -------
    list
        List with RGB percentages.

    """
    if num > 6:
        num = np.mod(num, 7)

    colors = [
        [0, 0.4470, 0.7410],
        [0.8500, 0.3250, 0.0980],
        [0.9290, 0.6940, 0.1250],
        [0.4940, 0.1840, 0.5560],
        [0.4660, 0.6740, 0.1880],
        [0.3010, 0.7450, 0.9330],
        [0.6350, 0.0780, 0.1840],
    ]

    return colors[np.round(num)]


def specs_plot(specs: list, **params):
    """
    Plots several spectra in line plot.

    Parameters
    ----------
    specs : list
        List of Spm objects that are generated from DAT files.
    **params : TODO

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        Spectra figure.

    """
    # plot spectra in list

    if "channelx" in params:
        channelx = params["channelx"]
    else:
        channelx = specs[0].channels[0]

    if "channely" in params:
        channely = params["channely"]
    else:
        channely = specs[0].channels[1]

    if "direction" in params:
        direction = params["direction"]
    else:
        direction = list(specs[0].data.keys())[0]

    if "color" in params:
        color = params["color"]
    else:
        color = pl.cm.rainbow(np.linspace(0, 1, len(specs)))
    if "print_legend" in params:
        print_legend = params["print_legend"]
    else:
        print_legend = False
    if "offset" in params:
        offset = params["offset"]
    else:
        offset = 0

    fig = plt.figure(figsize=(6, 4))

    counter = 0
    for s, c in zip(specs, color):
        (x_data, x_unit) = s.get_channel(channelx, direction=direction)
        (y_data, y_unit) = s.get_channel(channely, direction=direction)

        plt.plot(x_data, y_data + counter * offset, color=c, label=s.name)
        counter = counter + 1

    plt.xlabel(f"{channelx} ({x_unit})")
    plt.ylabel(f"{channely} ({y_unit})")

    if print_legend:
        plt.legend()

    plt.show()

    return fig


def ref_spec_plotting(
    ref_file: str, spec_files: list, fname_ref: str, fname_specs: str, **params
):
    """
    Plots several spectra in line plot together with respective spectra locations.

    Parameters
    ----------
    ref_file : str
        Filename of reference sxm image
    spec_files : list
        List of filenames of dat spectra.
    fname_ref : str
        Filename of saved reference image.
    fname_specs : str
        Filename of saved spectra image.
    **params : TYPE
        channelx_plot: channel to plot on x-axis of dat spectra (default: V)
        channely_plot: channel to plot on y-axis of dat spectra (default: dIdV)
        offset: offset between spectra (default: 0)
        annotate_spec: annotate spectra with numbers.

    Returns
    -------
    specs_fig : matplotlib.pyplot.Figure
        Spectra figure with respective spectra locations.

    """

    # Load modules

    # Optional input
    if "channelx_plot" in params:
        channelx_plot = params["channelx_plot"]
    else:
        channelx_plot = "V"

    if "channely_plot" in params:
        channely_plot = params["channely_plot"]
    else:
        channely_plot = "dIdV"

    if "print_legend" not in params:
        params["print_legend"] = True

    if "annotate_spec" in params:
        annotate_spec = params["annotate_spec"]
    else:
        annotate_spec = False

    ## Import
    # Importing reference image
    ref = spm(ref_file)

    # Importing spec
    sp = []
    for s in spec_files:
        sp.append(spm(s))

    # Defining color scale
    col = pl.cm.rainbow(np.linspace(0, 1, len(sp)))

    ## Plotting reference image with spec locations
    ref.plot(show_params=True, offset=False, show=False, channel="z", close_fig=False)

    # plot circle for each location
    for s, c in zip(sp, col):
        (x, y) = relative_position(ref, s)
        plt.plot(x, y, "o", color=c)
        if annotate_spec:
            plt.annotate(s.name, (x, y))

    plt.xlabel("x (nm)")
    plt.ylabel("y (nm)")

    fig_dir = os.path.abspath(os.path.join(fname_ref, os.pardir))
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    plt.savefig(fname_ref + ".png", dpi=500)
    plt.savefig(fname_ref + ".svg", dpi=500)

    # Plotting specs
    specs_fig = specs_plot(
        sp,
        channelx=channelx_plot,
        channely=channely_plot,
        direction="forward",
        color=col,
        **params,
    )

    fig_dir = os.path.abspath(os.path.join(fname_specs, os.pardir))
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    specs_fig.savefig(fname_specs + ".png", dpi=500)
    specs_fig.savefig(fname_specs + ".svg", dpi=500)

    plt.show()

    return specs_fig
