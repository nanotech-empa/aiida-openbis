import os
import re
import warnings

import matplotlib.pyplot as plt
import numpy as np
import spiepy
import yaml
from matplotlib.colors import LogNorm

from . import nanonispy as nap


class Spm:
    DEFAULT_CONFIG = os.path.dirname(os.path.abspath(__file__)) + "/machine_config.yaml"

    # constructor
    def __init__(self, path: str, config_file=DEFAULT_CONFIG):
        """
        Spm constructor

        Parameters
        ----------
        path : str
            SPM filepath.

        Returns
        -------
        None.

        """

        config = self.load_machine_configuration_from_yaml_file(config_file)

        # self.path = path.replace('//', '/')
        abspath = os.path.abspath(path)
        self.path = "/".join(abspath.split("\\")[-4:])
        self.name = self.path.split("/")[-1]
        file_extension = os.path.splitext(path)[1]

        if file_extension == ".sxm":
            self.napImport = nap.read.Scan(path)
            self.type = "scan"
        elif file_extension == ".dat":
            self.napImport = nap.read.Spec(path)
            self.type = "spec"
        elif file_extension == ".3ds":
            self.napImport = nap.read.Grid(path)
            self.type = "grid"
        else:
            print("Datatype not supported.")
            return
        
        self.SignalsList = []
        for key in self.napImport.signals.keys():
            in_config = key in config["Channels"]
            self.SignalsList.append(
                {
                    "ChannelName": key,
                    "ChannelNickname": config["Channels"][key]["nickname"]
                    if in_config
                    else key,
                    "ChannelScaling": config["Channels"][key]["scaling"]
                    if in_config
                    else 1.0,
                    "ChannelUnit": config["Channels"][key]["unit"]
                    if in_config
                    else "N/A",
                }
            )

        self.channels = [c["ChannelNickname"] for c in self.SignalsList]

        self.ParamListReference = []
        for key in config["Parameters"].keys():
            self.ParamListReference.append(
                {
                    "ParamName": key,
                    "ParamNickname": config["Parameters"][key]["nickname"],
                    "ParamScaling": config["Parameters"][key]["scaling"],
                    "ParamUnit": config["Parameters"][key]["unit"],
                }
            )

        self.header = self.napImport.header

    def load_machine_configuration_from_yaml_file(self, filepath):
        """
        Load the machine configuration from a yaml file.

        Inputs:
        -------
        filepath: str
            Path to the yaml file

        Outputs:
        --------
        config: dict
            Dictionary containing the machine configuration
        """

        # Define specific loader for float values
        mode = "r"
        loader = yaml.SafeLoader
        loader.add_implicit_resolver(
            "tag:yaml.org,2002:float",
            re.compile(
                """^(?:
            [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
            |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
            |\\.[0-9_]+(?:[eE][-+][0-9]+)?
            |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
            |[-+]?\\.(?:inf|Inf|INF)
            |\\.(?:nan|NaN|NAN))$""",
                re.X,
            ),
            list("-+0123456789."),
        )

        # Load the yaml file
        with open(filepath, mode=mode) as fhandle:
            config = yaml.load(fhandle, Loader=loader)

        # Check if the machine configuration is set to THzSTM
        if config["MachineConfig"] != "THzSTM":
            warnings.warn("MachineConfig not set to THzSTM")

        return config

    def __repr__(self):
        return self.path

    # get channel
    def get_channel(
        self,
        channel: str,
        direction: str = "forward",
        flatten: bool = False,
        offset: bool = False,
        zero: bool = False,
    ):
        """
        Returns the measurement values and the channel unit.

        Parameters
        ----------
        channel : str
            Name of the channel stored in the reference config file or in the NANONIS file.
        direction : str, optional
            Measurement direction. It can only be forward or backward. The default is 'forward'.
        flatten : bool, optional
            Substract background plane from image data. The default is False.
        offset : bool, optional
            Substract constant offset = mean value. The default is False.
        zero : bool, optional
            Substract constant offset = minimum value. The default is False.

        Returns
        -------
        tuple
            Measurement values and channel unit

        """

        if self.type == "scan":
            chNum = [d["ChannelNickname"] for d in self.SignalsList].index(channel)
            im = self.napImport.signals[self.SignalsList[chNum]["ChannelName"]][
                direction
            ]
            im = im * self.SignalsList[chNum]["ChannelScaling"]

            if flatten:
                if ~np.isnan(np.sum(im)):
                    im, _ = spiepy.flatten_xy(im)
                else:
                    m, n = np.shape(im)
                    i = np.argwhere(np.isnan(im))[0, 0]
                    im_cut = im[: i - 1, :]
                    im, _ = spiepy.flatten_xy(im_cut)
                    empty = np.full((m - i, n), np.nan)
                    im = np.vstack((im, empty))

                    warnings.warn("NaN values in image")
                    # im-np.nanmean(im)

            if offset:
                im = im - np.mean(im)
            if zero:
                im = im + abs(np.min(im))

            unit = self.SignalsList[chNum]["ChannelUnit"]

            return (im, unit)

        elif self.type == "spec":
            if direction == "backward":
                channel = channel + "_bw"
                # print(channel)

            chNum = [d["ChannelNickname"] for d in self.SignalsList].index(channel)
            data = self.napImport.signals[self.SignalsList[chNum]["ChannelName"]]
            data = data * self.SignalsList[chNum]["ChannelScaling"]
            unit = self.SignalsList[chNum]["ChannelUnit"]

            return (data, unit)

        elif self.type == "grid":
            if direction == "backward":
                channel = channel + "_bw"
                # print(channel)

            chNum = [d["ChannelNickname"] for d in self.SignalsList].index(channel)
            data = self.napImport.signals[self.SignalsList[chNum]["ChannelName"]]
            data = data * self.SignalsList[chNum]["ChannelScaling"]
            unit = self.SignalsList[chNum]["ChannelUnit"]

            return (data, unit)

        return

    # get parameter
    def get_param(self, param: str):
        """
        Returns the parameter value and the parameter unit.

        Parameters
        ----------
        param : str
            Name of the parameter stored in the reference config file or in the NANONIS file.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if any(d["ParamNickname"] == param for d in self.ParamListReference):
            paNum = [d["ParamNickname"] for d in self.ParamListReference].index(param)

            if self.ParamListReference[paNum]["ParamScaling"] is None:
                return self.napImport.header[
                    self.ParamListReference[paNum]["ParamName"]
                ]
            else:
                return (
                    float(
                        self.napImport.header[
                            self.ParamListReference[paNum]["ParamName"]
                        ]
                    )
                    * self.ParamListReference[paNum]["ParamScaling"],
                    self.ParamListReference[paNum]["ParamUnit"],
                )

        elif param == "width" or param == "height":
            # height, width = self.get_param('scan_range')
            # height, width = [height*10**9, width*10**9]
            scanfield = self.get_param("scan>scanfield")
            width, height = (
                float(scanfield.split(";")[2]) * 10**9,
                float(scanfield.split(";")[3]) * 10**9,
            )

            return (eval(param), "nm")

        else:
            if param in self.napImport.header.keys():
                return self.napImport.header[param]
            else:
                return

    # print essential parameters for plotting
    def print_params(self, show: bool = True):
        """
        Print essential parameters for plotting

        Parameters
        ----------
        show : bool, optional
            Whether label if printed. The default is True.

        Returns
        -------
        label : list
            Set of important parameters.

        """

        label = []

        if self.type == "scan":
            parameters = {
                "fb_enable": "z-controller>controller status",
                "bias": "V",
                "set_point": "setpoint",
                "z_offset": "z_offset",
                "comment": "comments",
                "height": "height",
                "width": "width",
                "angle": "angle",
                "controller": "z-controller>controller name",
                "date": "rec_date",
                "time": "rec_time",
            }
            for key, value in parameters.items():
                try:
                    parameters[key] = self.get_param(value)
                except Exception as err:
                    if type(err) == type(KeyError()):
                        print(f"Missing Header Key, {err}")
                    parameters[key] = "N/A"

            fb_enable = parameters["fb_enable"]

            bias = parameters["bias"]
            set_point = parameters["set_point"]
            z_offset = parameters["z_offset"]
            date = parameters["date"]
            time = parameters["time"]
            comment = parameters["comment"]
            parameters["controller"]
            parameters["height"]
            parameters["width"]
            parameters["angle"]
            if fb_enable == "OFF":
                label.append("constant height")
                label.append(f"z-offset: {z_offset[0]:.3f}{z_offset[1]}")

            if np.abs(bias[0]) < 0.1:
                bias = list(bias)
                bias[0] = bias[0] * 1000
                bias[1] = "mV"
                bias = tuple(bias)

            label.append(f"I = {set_point[0]:.0f}{set_point[1]}")
            label.append(f"bias = {bias[0]:.2f}{bias[1]}")
            # label.append(f'size: {width[0]}{width[1]} x {height[0]:.1f}{height[1]} ({angle[0]:.0f}{angle[1]})')
            label.append(f"comment: {comment}")
            label.append(f"Date: {date} - {time}")

        elif self.type == "spec":
            parameters = {
                "fb_enable": "Z-Ctrl hold",
                "bias": "V_spec",
                "set_point": "setpoint_spec",
                "z_offset": "Bias Spectroscopy>Z offset (m)",
                "comment": "comments",
                "lockin_amplitude": "lockin_amplitude",
                "lockin_phase": "lockin_phase",
                "lockin_frequency": "lockin_frequency",
                "lockin_status": "Lock-in>Lock-in status",
                "date": "Saved Date",
            }
            for key, value in parameters.items():
                try:
                    parameters[key] = self.get_param(value)
                except Exception as err:
                    if type(err) == type(KeyError()):
                        print(f"Missing Header Key, {err}")
                    parameters[key] = "N/A"

            fb_enable = parameters["fb_enable"]
            bias = parameters["bias"]
            set_point = parameters["set_point"]
            z_offset = parameters["z_offset"]
            date = parameters["date"]
            parameters["lockin_status"]
            lockin_amplitude = parameters["lockin_amplitude"]
            lockin_phase = parameters["lockin_phase"]
            lockin_frequency = parameters["lockin_frequency"]
            comment = parameters["comment"]

            # if lockin_status == 'ON':
            label.append(
                f"lockin: A = {lockin_amplitude[0]:.0f}{lockin_amplitude[1]} (θ = {lockin_phase[0]:.0f}{lockin_phase[1]}, f = {lockin_frequency[0]:.0f}{lockin_frequency[1]})"
            )

            if fb_enable == "FALSE":
                label.append("feedback on")

            elif fb_enable == "TRUE":
                label.append("feedback off")

            label.append(
                f"setpoint: I = {set_point[0]:.0f}{set_point[1]}, V = {bias[0]:.1f}{bias[1]}"
            )
            label.append(f"comment: {comment}")
            label.append(f"Date: {date}")

        label.append(f"path: {self.path}")
        label = "\n".join(label)

        if show:
            print(label)

        return label

    # plot
    def plot(self, **params):
        """
        Plots a channel of the Spm objects

        Parameters
        ----------
        **params : TYPE
            Optional parameters for .sxm files:
                channel: Name of the channel that is plotted
                direction: Measurement direction (forward or backward)
                flatten: Boolean whether the plane subtraction
                offset: Boolean whether the constant offset (mean value) should be subtracted
                cmap: Name of color map
                clim: Lower bound and upper bound of the color limits
                log: Boolean whether the color scale is logarithmic
                show_params: Whether in the title essential parameters are shown
                show: Whether the plot is shown
                save: Whether plot is saved
                save_name: Name of the saved plot
                close_fig: Whether figure is closed at the end
                save_format: Save format
                zero: Boolean whether the constant offset (minimum value) should be subtracted

            Optional parameters for .dat files:
                channelx: Name of the X channel
                channely: Name of the Y channel
                direction: Measurement direction (forward or backward)
                log: Boolean whether the color scale is logarithmic
                loglog: Log-log plot
                show_params: Whether in the title essential parameters are shown
                show: Whether the plot is shown
                save: Whether plot is saved
                save_name: Name of the saved plot
                close_fig: Whether figure is closed at the end
                save_format: Save format

        Returns
        -------
        fig : matplotlib.pyplot.Figure
            Figure of a channel of the Spm object

        """

        # cmaps = sorted(m for m in plt.cm.datad)
        # ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'jet', 'jet_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spectral', 'spectral_r', 'spring', 'spring_r', 'summer', 'summer_r', 'terrain', 'terrain_r', 'winter', 'winter_r']

        # ml.interactive(0)

        if self.type == "scan":
            if "channel" in params:
                channel = params["channel"]
            else:
                channel = self.channels[0]

            if "direction" in params:
                direction = params["direction"]
            else:
                direction = "forward"

            if "flatten" in params:
                flatten = params["flatten"]
            else:
                flatten = False

            if "offset" in params:
                offset = params["offset"]
            else:
                offset = False

            if "cmap" in params:
                cmap = params["cmap"]
            else:
                cmap = "gray"

            if "clim" in params:
                clim = params["clim"]
            else:
                clim = False

            if "log" in params:
                log = params["log"]
            else:
                log = False

            if "show_params" in params:
                show_params = params["show_params"]
            else:
                show_params = False

            if "show" in params:
                show = params["show"]
            else:
                show = True

            if "save" in params:
                save = params["save"]
            else:
                save = False

            if "save_name" in params:
                save_name = params["save_name"]
            else:
                save_name = False

            if "close_fig" in params:
                close_fig = params["close_fig"]
            else:
                close_fig = True

            if "save_format" in params:
                save_format = params["save_format"]
            else:
                save_format = "png"

            if "zero" in params:
                zero = params["zero"]
            else:
                zero = False

            (chData, chUnit) = self.get_channel(
                channel, direction=direction, flatten=flatten, offset=offset, zero=zero
            )

            if "vmin" in params:
                vmin = params["vmin"]
            else:
                vmin = np.min(chData)

            if "vmax" in params:
                vmax = params["vmax"]
            else:
                vmax = np.max(chData)

            if direction == "backward":
                chData = np.fliplr(chData)

            width = self.get_param("width")
            height = self.get_param("height")
            pix_y, pix_x = np.shape(chData)

            fig = plt.figure(figsize=(7, 8))

            ImgOrigin = "lower"
            if self.get_param("scan_dir") == "down":
                ImgOrigin = "upper"

            if log:
                im = plt.imshow(
                    np.abs(chData),
                    aspect="equal",
                    extent=[0, width[0], 0, height[0]],
                    cmap=cmap,
                    norm=LogNorm(),
                    origin=ImgOrigin,
                )
            else:
                im = plt.imshow(
                    chData,
                    aspect="equal",
                    extent=[0, width[0], 0, height[0]],
                    vmin=vmin,
                    vmax=vmax,
                    cmap=cmap,
                    origin=ImgOrigin,
                )

            if clim:
                plt.clim(clim)

            im.axes.set_xticks([0, width[0]])
            im.axes.set_xticklabels([0, np.round(width[0], 2)])
            im.axes.set_yticks([0, height[0]])
            im.axes.set_yticklabels([0, np.round(height[0], 2)])

            if show_params:
                title = self.print_params(show=False)
            else:
                title = self.path

            plt.title(title + "\n", loc="left")
            plt.xlabel(f"x ({width[1]})")
            plt.ylabel(f"y ({height[1]})")

            cbar = plt.colorbar(
                im, fraction=0.046, pad=0.02, format="%.2g", shrink=0.5, aspect=10
            )
            cbar.set_label(f"{channel} ({chUnit}")

            if show:
                plt.show()
            # else:
            #     plt.close(fig)

            if close_fig:
                plt.close(fig)

            if save:
                if save_name:
                    fname = save_name
                else:
                    fname = self.path.split("/")[-1]
                    fname = fname.split(".")[-2]
                fig.savefig(fname + "." + save_format, dpi=500)

            return fig

        elif self.type == "spec":
            if "channelx" in params:
                channelx = params["channelx"]
            else:
                channelx = self.channels[0]

            if "channely" in params:
                channely = params["channely"]
            else:
                channely = self.channels[1]

            if "direction" in params:
                direction = params["direction"]
            else:
                direction = direction = "forward"

            if "log" in params:
                log = params["log"]
            else:
                log = False

            if "loglog" in params:
                loglog = params["loglog"]
            else:
                loglog = False

            if "show_params" in params:
                show_params = params["show_params"]
            else:
                show_params = False

            if "show" in params:
                show = params["show"]
            else:
                show = True

            if "save" in params:
                save = params["save"]
            else:
                save = False

            if "save_name" in params:
                save_name = params["save_name"]
            else:
                save_name = False

            if "close_fig" in params:
                close_fig = params["close_fig"]
            else:
                close_fig = True

            if "save_format" in params:
                save_format = params["save_format"]
            else:
                save_format = "png"

            (x_data, x_unit) = self.get_channel(channelx, direction=direction)
            (y_data, y_unit) = self.get_channel(channely, direction=direction)

            fig = plt.figure(figsize=(6, 4))

            if log:
                plt.semilogy(x_data, np.abs(y_data))

            elif loglog:
                plt.loglog(np.abs(x_data), np.abs(y_data))

            else:
                plt.plot(x_data, y_data)

            if show_params:
                title = self.print_params(show=False)
            else:
                title = self.path

            plt.title(title + "\n", loc="left")

            plt.xlabel(f"{channelx} ({x_unit})")
            plt.ylabel(f"{channely} ({y_unit})")

            if show:
                plt.show()
            # else:
            #     plt.close(fig)

            if close_fig:
                plt.close(fig)

            if save:
                if save_name:
                    fname = save_name
                else:
                    fname = self.path.split("/")[-1]
                    fname = fname.split(".")[-2]
                fig.savefig(fname + "." + save_format, dpi=500)

            return fig

    def importall(FilePath: str, FilePrefix="", ImportOnly=""):
        """
        Returns a list of Spm objects from a specific folder

        Parameters
        ----------
        FilePath : str
            File path of local directory containing Nanonis files
        FilePrefix : str
            Optional name of file base to be imported e.g. 'img_'
        ImportOnly : 'scan' or 'spec'
            Optional specifier for file type to be imported

        Returns
        -------
        list of Spm objects

        """

        from os import walk

        files = []
        NumSpec = 0
        NumScan = 0
        NumGrid = 0

        for root, dirs, filenames in walk(FilePath):
            for file in filenames:
                if file.endswith(".sxm") and file.startswith(FilePrefix):
                    if not ((ImportOnly == "spec") or (ImportOnly == "grid")):
                        files.append(Spm(root + "/" + file))
                        NumScan = NumScan + 1
                elif file.endswith(".dat") and file.startswith(FilePrefix):
                    if not ((ImportOnly == "scan") or (ImportOnly == "grid")):
                        files.append(Spm(root + "/" + file))
                        NumSpec = NumSpec + 1
                elif file.endswith(".3ds") and file.startswith(FilePrefix):
                    if not ((ImportOnly == "scan") or (ImportOnly == "spec")):
                        files.append(Spm(root + "/" + file))
                        NumGrid = NumGrid + 1

        print(
            str(len(files))
            + " files imported; "
            + str(NumScan)
            + " scan(s), "
            + str(NumSpec)
            + " spectra, and "
            + str(NumGrid)
            + " grids"
        )

        return files
