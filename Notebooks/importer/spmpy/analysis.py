import numpy as np


def fit_parabola(
    bias,
    df,
    single_spectrum: bool = False,
    fit_min: bool = False,
    fit_max: bool = False,
):
    """
    Fits parabolas to a list of df vs Bias spectrum

    Parameters
    ----------
    bias : numpy.array
        1D or 2D array of bias voltage values
    df : numpy.array
        1D or 2D array of corresponding frequency shift values.
    single_spectrum : bool
        If bias and df is only one spectrum
    fitMin : TODO
        TODO
    fitMax : TODO
        TODO
    **params : TODO

    Returns
    -------
    tuple
        p_fit: List of second degree polyfits (p_fit = [a,b,c], where y = ax**2 + bx + c)
        err_p: List of uncertainties of the polyfit
        df_fit: List of Y values of the fitted parabola
        bias_max: List of X coordinate of the extrema of each parabolas
        err_bias_max: List of uncertainties of the X coordinate of the extrema of each parabolas
        df_max: List of Y coordinate of the extrema of each parabolas
        err_df_max: List of uncertainties of the Y coordinate of the extrema of each parabolas

    """

    if single_spectrum:
        bias = np.array([bias])
        df = np.array([df])
    p_fit = np.zeros((len(bias), 3))
    df_fit = np.zeros((len(bias), len(bias[0])))
    err_p = np.zeros((len(bias), 3))
    bias_max = np.zeros(len(bias))
    err_bias_max = np.zeros(len(bias))
    df_max = np.zeros(len(bias))
    err_df_max = np.zeros(len(bias))

    for i in range(len(bias)):
        if not fit_min and not fit_max:
            (p, cov) = np.polyfit(bias[i], df[i], 2, cov=True)
        elif fit_min and not fit_max:
            fit_i = np.where(bias[i] > fit_min[i])[0]
            (p, cov) = np.polyfit(bias[i, fit_i], df[i, fit_i], 2, cov=True)
        p_fit[i] = p
        err1 = np.sqrt(abs(cov))
        err_p[i] = np.array([err1[0, 0], err1[1, 1], err1[2, 2]])

        df_fit[i] = np.polyval(p, bias[i])

        bias_max[i] = -p[1] / (2 * p[0])

        # if errors independent

        # err_bias_max[i] = abs(bias_max[i]) * np.sqrt((err1[1, 1]
        #        / p[1]) ** 2 + (err1[0, 0] / p[0]) ** 2)
        err_bias_max[i] = np.sqrt(
            (err1[1, 1] / (2 * p[0])) ** 2 + (err1[0, 0] * p[1] / (2 * p[0] ** 2)) ** 2
        )

        # maximum error in any case (John Tayler Error analysis book)

        # err_bias_max[i] = abs(bias_max[i]) * (err1[1, 1] + err1[0, 0])

        df_max[i] = np.polyval(p, bias_max[i])
        err_df_max[i] = (
            p[0]
            * bias_max[i] ** 2
            * (err1[0, 0] + 2 * err_bias_max[i] / abs(bias_max[i]))
            + p[1] * bias_max[i] * (err1[1, 1] + err_bias_max[i] / abs(bias_max[i]))
            + p[2] * err1[2, 2]
        )
    return (
        p_fit,
        err_p,
        df_fit,
        bias_max,
        err_bias_max,
        df_max,
        err_df_max,
    )


def relative_position(img, spec, **params):
    """
    Returns the (X,Y) position of the spectrum with respect to the position of the image origin.

    Parameters
    ----------
    img : Spm
        Spm object of the SXM file.
    spec : Spm
        Spm object of the DAT file.
    **params : TODO

    Returns
    -------
    list
        Coordinates of the position of the spectrum with respect to the position of the image origin.

    """

    # width = ref.get_param('width')
    # height = ref.get_param('height')
    # [px_x,px_y] = ref.get_param('scan_pixels')

    [o_x, o_y] = img.get_param("scan_offset")
    width = img.get_param("width")[0]
    height = img.get_param("height")[0]
    [o_x, o_y] = [o_x * 10**9, o_y * 10**9]

    angle = float(img.get_param("scan_angle")) * -1 * np.pi / 180

    x_spec = spec.get_param("x")[0]
    y_spec = spec.get_param("y")[0]

    if angle != 0:
        # Transforming to relative coordinates with angle
        x_rel = (
            (x_spec - o_x) * np.cos(angle) + (y_spec - o_y) * np.sin(angle) + width / 2
        )
        y_rel = (
            -(x_spec - o_x) * np.sin(angle)
            + (y_spec - o_y) * np.cos(angle)
            + height / 2
        )
    else:
        x_rel = x_spec - o_x + width / 2
        y_rel = y_spec - o_y + height / 2
    return [x_rel, y_rel]


def kpfm(files: list, **params):
    """
    Returns dict of fitted kpfm parabolas for list of spm objects

    Parameters
    ----------
    files : list
        List of Spm objects.
    **params : TODO

    Returns
    -------
    data: dict
        Dictionary of the parabola fit parameters returned by fit_parabola function.

    """

    if "range" in params:
        range = params["range"]
    else:
        range = [0, len(files[0].get_channel("V")[0]) - 1]

    data = {
        "V": [],
        "df": [],
        "V_max": [],
        "df_max": [],
        "p_fit": [],
        "err_p": [],
        "V_fit": [],
        "df_fit": [],
        "err_V_max": [],
        "err_df_max": [],
        "position": [],
    }

    for f in files:
        if f.type == "spec":
            if any(d["ChannelNickname"] == "df" for d in f.SignalsList):
                data["V"].append(f.get_channel("V")[0])
                data["df"].append(f.get_channel("df")[0])

            if any(d["ChannelNickname"] == "df_bw" for d in f.SignalsList):
                if not ("df_bw" in data.keys()):
                    data["V_bw"] = []
                    data["df_bw"] = []

                else:
                    data["V_bw"].append(f.get_channel("V_bw")[0])
                    data["df_bw"].append(f.get_channel("df_bw")[0])

            x = f.get_param("x")[0]
            y = f.get_param("y")[0]

        data["position"].append([x, y])

    (p_fit, err_p, df_fit, V_max, err_V_max, df_max, err_df_max) = fit_parabola(
        np.array(data["V"])[:, range[0] : range[1]].tolist(),
        np.array(data["df"])[:, range[0] : range[1]].tolist(),
        single_spectrum=False,
        fit_min=False,
        fit_max=False,
    )

    data["p_fit"] = p_fit
    data["err_p"] = err_p
    data["V_fit"] = np.array(data["V"])[:, range[0] : range[1]].tolist()
    data["df_fit"] = df_fit
    data["V_max"] = V_max
    data["err_V_max"] = err_V_max
    data["df_max"] = df_max
    data["err_df_max"] = err_df_max

    if "df_bw" in data.keys():
        (
            p_fit_bw,
            err_p_bw,
            df_fit_bw,
            v_max_bw,
            err_v_max_bw,
            df_max_bw,
            err_df_max_bw,
        ) = fit_parabola(
            np.array(data["V_bw"])[:, range[0] : range[1]].tolist(),
            np.array(data["df_bw"])[:, range[0] : range[1]].tolist(),
            single_spectrum=False,
            fit_min=False,
            fit_max=False,
        )

        data["p_fit_bw"] = p_fit_bw
        data["err_p_bw"] = err_p_bw
        data["V_fit_bw"] = np.array(data["V_bw"])[:, range[0] : range[1]].tolist()
        data["df_fit_bw"] = df_fit_bw
        data["V_max_bw"] = v_max_bw
        data["err_V_max_bw"] = err_v_max_bw
        data["df_max_bw"] = df_max_bw
        data["err_df_max_bw"] = err_df_max_bw

    return data
