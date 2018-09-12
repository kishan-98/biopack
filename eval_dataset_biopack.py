# Usage: python eval_dataset_biopack.py <path to text file to be analysed> [<location to store graphs> or <locations to store graphs of EDA & PPG respectively at different places>]
import sys, os, re, importlib, statistics, math, operator
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, integrate, interpolate
from numpy import NaN, Inf, arange, isscalar, asarray, array

def find_peaks(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []

    if x is None:
        x = arange(len(v))

    v = asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN

    lookformax = False

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append(mxpos)
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append(mnpos)
                mx = this
                mxpos = x[i]
                lookformax = True
    return maxtab, mintab
    # return array(maxtab + mintab)

def slope_threshold(v, i1, i2, threshold=0.093, dx=1):
    '''
    v is entire array of data; i1 and i2 are range of indices of min peak and max peak respectively; threshold is cutoff for the first order derivative for required peaks
    '''
    for i in range(i1+1, i2+1):
        if (v[i] - v[i-1])/dx > threshold:
            return True
    return False

def filter_peaks(v_y, v_min, v_max, threshold_slope=0.093, threshold_time=2, dx=1):
    '''
    v_y is entire array of data; v_max and v_min are list of indices of max peak and min peak respectively; threshold is cutoff for the first order derivative for required peaks
    '''
    filtered_v_max, filtered_v_min = [], []
    for i1, i2 in zip(v_min, v_max):
        # Slope based filtering - threshold 0.093 microSeimens per second
        if slope_threshold(v_y, i1, i2, threshold_slope, dx):
            filtered_v_min.append(i1)
            filtered_v_max.append(i2)
    # return filtered_v_max, filtered_v_min
    v_max, v_min = filtered_v_max, filtered_v_min
    v_len = len(v_max)
    if not v_len:
        return filtered_v_max, filtered_v_min
    filtered_v_max, filtered_v_min = [v_max[0]], [v_min[0]]
    for i in range(1, v_len):
        # Proximity based filtering - threshold 2 second(0.5 Hz)
        if v_min[i] - filtered_v_max[-1] > threshold_time/dx:# and v_min[i] - filtered_v_min[-1] > threshold_time/dx:
            filtered_v_min.append(v_min[i])
            filtered_v_max.append(v_max[i])
        elif v_y[filtered_v_max[-1]] < v_y[v_max[i]]:
            filtered_v_max[-1] = v_max[i]
    return filtered_v_max, filtered_v_min

def get_baseline(v, v_min, v_max):
    '''
    v is entire array of data; v_max and v_min are list of indices of max peak and min peak respectively;
    '''
    baseline = []
    k = min(len(v_min), len(v_max)) - 1
    # Before first minimum, the baseline value is equal to the first minimum
    for i in range(v_min[0]):
        baseline.append(v[v_min[0]])
    # For all the [minimum, maximum] value pair, the temporary baseline is considered to be the linear line between the ith minimum and (i+1)th minimum
    # The actual baseline is the ratio of area under curve of the EDA data wrt to temporary baseline
    # This sequence can be repeated n number of times to get accurate baseline
    # Currently it is repeated 1 time only
    for i in range(0, k):
        diff_min_values = v[v_min[i+1]] - v[v_min[i]]
        temp_baseline = [v[v_min[i]] + ((j - v_min[i])/((v_min[i+1] - v_min[i])))*(diff_min_values) for j in range(v_min[i], v_min[i+1] + 1)]
        # total_area = integrate.simps([v[j+v_min[j]] - temp_baseline[j] for j in range(v_min[i], v_min[i+1] + 1)])
        total_area = np.trapz(list(map(abs, [v[j] - temp_baseline[j-v_min[i]] for j in range(v_min[i], v_min[i+1] + 1)])))
        curr_area = 0
        for j in range(v_min[i], v_min[i+1]):
            # baseline.append(temp_baseline[j-v_min[i]])
            baseline.append(v[v_min[i]] + (curr_area/total_area)*(diff_min_values))
            curr_area = curr_area + abs(np.trapz([v[j] - temp_baseline[j-v_min[i]], v[j+1] - temp_baseline[j+1-v_min[i]]]))
        # print("for i =", i," curr_area =", curr_area, " and total_area =", total_area)

    # For all the last [minimum, maximum] value pair, it is handled little differently, but with the same concept
    temp_min_index, temp_min = min(enumerate(v[v_max[k]:]), key = operator.itemgetter(1))
    temp_min_index = temp_min_index + v_max[k]
    diff_min_values = temp_min - v[v_min[k]]
    temp_baseline = [v[v_min[k]] + ((j - v_min[k])/((temp_min_index - v_min[k])))*(diff_min_values) for j in range(v_min[k], temp_min_index + 1)]
    # total_area = integrate.simps([v[j+v_min[j]] - temp_baseline[j] for j in range(v_min[i], v_min[i+1] + 1)])
    total_area = np.trapz(list(map(abs, [v[j] - temp_baseline[j-v_min[k]] for j in range(v_min[k], temp_min_index + 1)])))
    curr_area = 0
    for j in range(v_min[k], temp_min_index):
        # baseline.append(temp_baseline[j-v_min[k]])
        baseline.append(v[v_min[k]] + (curr_area/total_area)*(diff_min_values))
        curr_area = curr_area + abs(np.trapz([v[j] - temp_baseline[j-v_min[k]], v[j+1] - temp_baseline[j+1-v_min[k]]]))
    # print("for i =", k," curr_area =", curr_area, " and total_area =", total_area)
    i = temp_min_index
    while i < len(v):
        baseline.append(temp_min)
        i = i + 1
    # print("len(v) == len(baseline):", len(v) == len(baseline))
    return baseline

def get_details_EDA(eda_y, eda_x, sampling_rate_eda, display_details=True, display_graph=False, graph_title="", store_graph=False, store_location=""):
    EDA = {"x" : eda_x, "y" : eda_y, "tonic" : [], "mean" : 0, "median" : 0, "stdev" : 0, "peak_count" : 0, "max_peak_indices" : [], "max_peak_position" : [], "max_peak_values" : [], "min_peak_indices" : [], "min_peak_position" : [], "min_peak_values" : [], "peak_values" : [], "area" : 0}
    # print("length:", len(eda_y), len(eda_x))
    if len(eda_y) < 1 or len(eda_x) < 1:
        return EDA
    # EDA
    EDA["mean"]   = statistics.mean(EDA["y"])
    EDA["median"] = statistics.median(EDA["y"])
    EDA["stdev"]  = statistics.stdev(EDA["y"]) if len(EDA["y"]) > 1 else EDA["y"]
    # window   = signal.general_gaussian(51, p=0.5, sig=20)
    # filtered = signal.fftconvolve(window, EDA["y"])
    # filtered = (np.average(EDA["y"]) / np.average(filtered)) * filtered
    # filtered = np.roll(filtered, -25)
    # EDA["peak_indices"] = signal.find_peaks_cwt(filtered, np.arange(10, 15), noise_perc=0.1).tolist()
    # EDA["peak_indices"] = signal.find_peaks_cwt(EDA["y"], np.arange(1, 2)).tolist()
    EDA["max_peak_indices"], EDA["min_peak_indices"] = find_peaks(EDA["y"], 0.01)
    EDA["peak_count"]    = min(len(EDA["max_peak_indices"]), len(EDA["min_peak_indices"]))
    EDA["max_peak_indices"] = EDA["max_peak_indices"][:EDA["peak_count"]]
    EDA["min_peak_indices"] = EDA["min_peak_indices"][:EDA["peak_count"]]
    EDA["max_peak_indices"], EDA["min_peak_indices"] = filter_peaks(EDA["y"], EDA["min_peak_indices"], EDA["max_peak_indices"], threshold_slope=0.093, threshold_time=2, dx=1/sampling_rate_eda)
    EDA["peak_count"]    = min(len(EDA["max_peak_indices"]), len(EDA["min_peak_indices"]))
    EDA["max_peak_position"] = [EDA["x"][x] for x in EDA["max_peak_indices"]]
    EDA["min_peak_position"] = [EDA["x"][x] for x in EDA["min_peak_indices"]]
    EDA["max_peak_values"] = [EDA["y"][x] for x in EDA["max_peak_indices"]]
    EDA["min_peak_values"] = [EDA["y"][x] for x in EDA["min_peak_indices"]]
    EDA["peak_values"]   = [EDA["y"][EDA["max_peak_indices"][i]] - EDA["y"][EDA["min_peak_indices"][i]] for i in range(EDA["peak_count"])]
    # EDA["tonic"] = interpolate.InterpolatedUnivariateSpline(EDA["min_peak_position"], EDA["min_peak_values"])(EDA["x"]) if EDA["peak_count"] > 3 else EDA["y"]
    EDA["tonic"] = get_baseline(EDA["y"], EDA["min_peak_indices"], EDA["max_peak_indices"])
    for i in range(EDA["peak_count"]):
        temp_y = [EDA["y"][x] - EDA["y"][EDA["min_peak_indices"][i]] for x in range(EDA["min_peak_indices"][i], EDA["max_peak_indices"][i] + 1)]
        EDA["area"] = EDA["area"] + 2*integrate.simps(temp_y, dx=1/sampling_rate_eda)
    # print("EDA[\"peak_values\"]:", EDA["peak_values"])
    if display_details:
        print("mean,", EDA["mean"])
        print("median,", EDA["median"])
        print("stdev,", EDA["stdev"])
        print("area,", EDA["area"])
        print("peak_count,", EDA["peak_count"])
        # print(", ".join(map(str, EDA["peak_values"])))
        print("average_peak, %s," % statistics.mean(EDA["peak_values"]), "max_peak, %s," % max(EDA["peak_values"]), "min_peak, %s" % min(EDA["peak_values"])) if EDA["peak_count"] > 0 else print("average_peak, ,", "max_peak, ,", "min_peak, ")
    if display_graph:
        plt.plot(EDA["x"], EDA["y"], 'r')
        plt.scatter(EDA["max_peak_position"], EDA["max_peak_values"], c='b')
        plt.scatter(EDA["min_peak_position"], EDA["min_peak_values"], c='g')
        plt.plot(EDA["x"], EDA["tonic"], "k-", label="baseline")
        plt.title(graph_title)
        plt.xlabel("Time(minutes)")
        plt.ylabel("Conductance(µSiemens)")
        plt.legend()
        plt.show()
        plt.gcf().clear()
    if store_graph:
        plt.plot(EDA["x"], EDA["y"], 'r')
        plt.scatter(EDA["max_peak_position"], EDA["max_peak_values"], c='b')
        plt.scatter(EDA["min_peak_position"], EDA["min_peak_values"], c='g')
        plt.title(graph_title)
        plt.xlabel("Time(minutes)")
        plt.ylabel("Conductance(µSiemens)")
        plt.savefig(store_location, format="png", dpi=1000)
        plt.gcf().clear()

    return EDA

def get_details_PPG(ppg_y, ppg_x, sampling_rate_ppg, display_details=True, display_graph=False, graph_title="", store_graph=False, store_location=""):
    PPG = {"x" : ppg_x, "y" : ppg_y, "mean" : 0, "median" : 0, "stdev" : 0, "peak_count" : 0, "max_peak_indices" : [], "max_peak_position" : [], "min_peak_indices" : [], "min_peak_position" : [], "peak_values" : [], "power_spectrum" : {"freq" : [], "density": []}, "total_energy" : 0}
    # print("length:", len(ppg_y), len(ppg_x))
    if len(ppg_y) < 1 or len(ppg_x) < 1:
        return PPG
    # PPG
    PPG["mean"]   = statistics.mean(PPG["y"])
    PPG["median"] = statistics.median(PPG["y"])
    PPG["stdev"]  = statistics.stdev(PPG["y"]) if len(PPG["y"]) > 1 else PPG["y"]
    # window   = signal.general_gaussian(51, p=0.5, sig=20)
    # filtered = signal.fftconvolve(window, PPG["y"])
    # filtered = (np.average(PPG["y"]) / np.average(filtered)) * filtered
    # filtered = np.roll(filtered, -25)
    # PPG["peak_indices"] = signal.find_peaks_cwt(filtered, np.arange(10, 15), noise_perc=0.1).tolist()
    # EDA["peak_indices"] = signal.find_peaks_cwt(PPG["y"], np.arange(1, 2)).tolist()
    PPG["max_peak_indices"], PPG["min_peak_indices"] = find_peaks(PPG["y"], 0.1)
    PPG["peak_count"]    = min(len(PPG["max_peak_indices"]), len(PPG["min_peak_indices"]))
    PPG["max_peak_indices"] = PPG["max_peak_indices"][:PPG["peak_count"]]
    PPG["min_peak_indices"] = PPG["min_peak_indices"][:PPG["peak_count"]]
    PPG["max_peak_position"] = [PPG["x"][x] for x in PPG["max_peak_indices"]]
    PPG["min_peak_position"] = [PPG["x"][x] for x in PPG["min_peak_indices"]]
    PPG["max_peak_values"] = [PPG["y"][x] for x in PPG["max_peak_indices"]]
    PPG["min_peak_values"] = [PPG["y"][x] for x in PPG["min_peak_indices"]]
    PPG["peak_values"]   = [PPG["y"][PPG["max_peak_indices"][i]] - PPG["y"][PPG["min_peak_indices"][i]] for i in range(PPG["peak_count"])]
    # print("PPG[\"peak_values\"]:", PPG["peak_values"])
    PPG["power_spectrum"]["freq"], PPG["power_spectrum"]["density"] = signal.periodogram(PPG["y"], sampling_rate_ppg)
    max_density_index, max_density = max(enumerate(PPG["power_spectrum"]["density"]), key = operator.itemgetter(1))
    max_frequency = PPG["power_spectrum"]["freq"][max_density_index]
    filtered_density = list(filter(lambda x: x[1] > 0.8 * max_density, enumerate(PPG["power_spectrum"]["density"])))
    filtered_frequency = list(map(lambda x: PPG["power_spectrum"]["freq"][x[0]], filtered_density))
    PPG["total_energy"] = sum(list(map(lambda x: x**2, PPG["y"])))/sampling_rate_ppg
    if display_details:
        print("mean,", PPG["mean"])
        print("median,", PPG["median"])
        print("stdev,", PPG["stdev"])
        print("peak_count,", PPG["peak_count"])
        # print(", ".join(map(str, PPG["peak_values"])))
        print("average_peak, %s," % statistics.mean(PPG["peak_values"]), "max_peak, %s," % max(PPG["peak_values"]), "min_peak, %s" % min(PPG["peak_values"])) if PPG["peak_count"] > 0 else print("average_peak, ,", "max_peak, ,", "min_peak, ")
        print("power_spectrum max frequency,", max_frequency)
        print("power_spectrum frequencies,",", ".join(map(str, filtered_frequency)))
        print("total_energy,", PPG["total_energy"])
    if display_graph:
        plt.plot(PPG["x"], PPG["y"], 'r')
        plt.scatter(PPG["max_peak_position"], PPG["max_peak_values"], c='b')
        plt.scatter(PPG["min_peak_position"], PPG["min_peak_values"], c='g')
        plt.title(graph_title)
        plt.xlabel("Time(minutes)")
        plt.ylabel("volts")
        plt.show()
        plt.gcf().clear()
    if store_graph:
        plt.plot(PPG["x"], PPG["y"], 'r')
        plt.scatter(PPG["max_peak_position"], PPG["max_peak_values"], c='b')
        plt.scatter(PPG["min_peak_position"], PPG["min_peak_values"], c='g')
        plt.title(graph_title)
        plt.xlabel("Time(minutes)")
        plt.ylabel("volts")
        plt.savefig(store_location, format="png", dpi=1000)
        plt.gcf().clear()

    return PPG

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except(ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - np.asarray(y[0]))
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - np.asarray(y[-1]))
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def compress_data(y, compr_factor):
    '''
    y is the data to be compressed
    compr_factor is defined as number of samples out of which one sample is to be picked
    '''
    start = 0
    end = min(compr_factor, len(y))
    y_compr = []
    while start < len(y):
        y_compr.append(statistics.mean(y[start:end]))
        start = end
        end = min(end + compr_factor, len(y))
    return y_compr

# print("command line arguments:", sys.argv[:])
program_name = sys.argv[0]
filepath = sys.argv[1]
# storage_location_PPG = storage_location_EDA = os.path.dirname(filepath)
storage_location_PPG = storage_location_EDA = ""
if len(sys.argv) == 3:
    storage_location_PPG = storage_location_EDA = sys.argv[2]
elif len(sys.argv) > 3:
    storage_location_PPG = sys.argv[2]
    storage_location_EDA = sys.argv[3]
if not re.search("/$", storage_location_PPG):
    storage_location_PPG = storage_location_PPG + "/"
if not re.search("/$", storage_location_EDA):
    storage_location_EDA = storage_location_EDA + "/"
# print("program_name:", program_name)
# print("filepath:", filepath)
if os.path.exists(filepath):
    # print()
    # print("filepath:", filepath)
    filename  = os.path.splitext(os.path.basename(filepath))[0] # filename without extension
    extension = os.path.splitext(os.path.basename(filepath))[1] # extension of the file
    file = open(filepath, "r")
    # print("for filepath ", filepath, ":")
    default_split = True
    first_line = next(file).strip()
    second_line = next(file).strip()
    if first_line.split(":")[0] == "question":
        default_split = False
    else:
        file.seek(0, 0)
    details = [next(file) for x in range(9)]
    # print("details:", details)
    sampling_rate_EDA = 1000 # Try to get it from the file
    sampling_rate_PPG = 1000 # Try to get it from the file
    compression_factor_PPG = 1000 # Adjust manually
    compression_factor_EDA = 1000 # Adjust manually
    ppg_x, ppg_y, eda_x, eda_y = ([] for x in range(4))
    for line in file:
        line = line.strip()
        # print("line:", line)
        line = line.split("\t")
        # print("line:", line)
        ppg_x.append(float(line[0]))
        ppg_y.append(float(line[1]))
        eda_x.append(float(line[0]))
        eda_y.append(float(line[2]))

    # Plotting graphs and making it smooth
    # ppg_smooth = savitzky_golay(ppg_y, 101, 4)
    # plt.plot(ppg_x, ppg_y, 'b')
    # plt.plot(ppg_x, ppg_smooth, 'r')
    # plt.title("PPG v/s PPG smoothened Data")
    # plt.show()
    #
    # eda_smooth = savitzky_golay(eda_y, 1001, 4)
    # plt.plot(eda_x, eda_y, 'b')
    # plt.plot(eda_x, eda_smooth, 'r')
    # plt.title("EDA v/s EDA smoothened Data")
    # plt.show()

    # manager = plt.get_current_fig_manager()
    # manager.window.showMaximized()
    # manager.resize(*manager.window.maxsize())
    # manager.frame.Maximize(True)

    fig = plt.gcf()
    fig.set_size_inches(10, 7.5, forward=True)

    # Smoothing of graph
    # plt.plot(ppg_x, ppg_y, 'r', label="raw")
    # ppg_y = savitzky_golay(ppg_y, 101, 5)
    ppg_y = signal.wiener(ppg_y, mysize=101)
    # plt.plot(ppg_x, ppg_y, 'g', label="smooth")
    # plt.title("PPG data before v/s after smoothening")
    # plt.legend()
    # plt.xlabel("Time(minutes)")
    # plt.ylabel("volts")
    # plt.show()
    # plt.savefig(storage_location + filename + "_PPG(overall).png", format="png", dpi=1000)
    # plt.gcf().clear()
    # plt.plot(eda_x, eda_y, 'r', label="raw")
    # eda_y_s = savitzky_golay(eda_y, 2001, 5)
    # eda_y_s = savitzky_golay(eda_y_s, 2001, 5)
    # eda_y_s = savitzky_golay(eda_y_s, 2001, 5)
    # eda_y_s = savitzky_golay(eda_y_s, 2001, 5)
    eda_y = signal.wiener(eda_y, mysize=3001)
    # eda_y = signal.wiener(eda_y, mysize=3001)
    # eda_y = signal.wiener(eda_y, mysize=3001)
    # plt.plot(eda_x, eda_y, 'b', label="wiener smooth")
    eda_y = savitzky_golay(eda_y, 1001, 5)
    # plt.plot(eda_x, eda_y, 'g', label="weiner + golay")
    # plt.title("EDA data before v/s after smoothening")
    # plt.legend()
    # plt.xlabel("Time(minutes)")
    # plt.ylabel("Conductance(µSeimens)")
    # plt.show()
    # plt.savefig(storage_location + filename + "_EDA(overall).png", format="png", dpi=1000)
    # plt.gcf().clear()

    # sys.exit()

    if default_split:
        '''
        if no splitting instructions are provided, split into default(entire, one minute, half minute)
        '''

        orig_stdout = sys.stdout

        try:
            f = open(storage_location_EDA + filename + "_analysis_EDA" + extension, 'w+')
            sys.stdout = f
        except:
            print("Error in opening", storage_location_EDA + filename + "_analysis_EDA" + extension)
            print("Displaying on terminal")

        # Details
        print("Details:", "".join(details))

        # EDA Channel Data
        print("EDA Channel Data")
        # print("")
        # Overall Data
        print("Overall Data:")
        # print("")
        EDA = get_details_EDA(eda_y, eda_x, sampling_rate_EDA, display_details=True, display_graph=False, graph_title="Overall Data EDA", store_graph=False, store_location=storage_location_EDA + filename + "_EDA(overall).png")
        print("")
        # sys.exit()

        # Overall compressed data
        # eda_x = compress_data(eda_x, compression_factor_EDA)
        # eda_y = compress_data(eda_y, compression_factor_EDA)
        # sampling_rate_EDA = sampling_rate_EDA/compression_factor_EDA

        EDA_one, EDA_half = ([] for x in range(2))

        time_unit = 60
        # One Minute Data for EDA
        start_pointer = 0
        end_pointer = time_unit * sampling_rate_EDA
        print("One Minute Data:")
        # print("")
        while start_pointer < len(eda_x):
            print(start_pointer/sampling_rate_EDA, "-", min(end_pointer/sampling_rate_EDA, len(eda_x)/sampling_rate_EDA), "seconds")
            temp_eda = get_details_EDA(eda_y[int(start_pointer):int(min(end_pointer, len(eda_x)))], eda_x[int(start_pointer):int(min(end_pointer, len(eda_x)))], sampling_rate_EDA, display_details=True, display_graph=False, graph_title=str(start_pointer/sampling_rate_EDA) + "-" + str(end_pointer/sampling_rate_EDA) + "seconds(EDA)", store_graph=False, store_location=storage_location_EDA + filename + "_(" + str(start_pointer) + "-" + str(end_pointer) + ").png")
            EDA_one.append(temp_eda)
            start_pointer = end_pointer
            end_pointer = end_pointer + time_unit * sampling_rate_EDA
            # print("")
            # break
            print("")

        time_unit = 30
        # Half Minute Data for EDA
        start_pointer = 0
        end_pointer = time_unit * sampling_rate_EDA
        print("Half Minute Data:")
        # print("")
        while start_pointer < len(eda_x):
            print(start_pointer/sampling_rate_EDA, "-", min(end_pointer/sampling_rate_EDA, len(eda_x)/sampling_rate_EDA), "seconds")
            temp_eda = get_details_EDA(eda_y[int(start_pointer):int(min(end_pointer, len(eda_x)))], eda_x[int(start_pointer):int(min(end_pointer, len(eda_x)))], sampling_rate_EDA, display_details=True, display_graph=False, graph_title=str(start_pointer/sampling_rate_EDA) + "-" + str(end_pointer/sampling_rate_EDA) + "seconds(EDA)", store_graph=False, store_location=storage_location_EDA + filename + "_(" + str(start_pointer) + "-" + str(end_pointer) + ").png")
            EDA_half.append(temp_eda)
            start_pointer = end_pointer
            end_pointer = end_pointer + time_unit * sampling_rate_EDA
            # print("")
            # break
            print("")

        try:
            f = open(storage_location_PPG + filename + "_analysis_PPG" + extension, 'w+')
            sys.stdout = f
        except:
            print("Error in opening", storage_location_PPG + filename + "_analysis_PPG" + extension)
            print("Displaying on terminal")

        # Details
        print("Details:", "".join(details))

        # PPG Channel Data
        print("PPG Channel Data")
        # print("")
        # Overall Data
        print("Overall Data:")
        # print("")
        PPG = get_details_PPG(ppg_y, ppg_x, sampling_rate_PPG, display_details=True, display_graph=False, graph_title="Overall Data PPG", store_graph=False, store_location=storage_location_PPG + filename + "_PPG(overall).png")
        print("")
        # Plotting and saving overall spectral power data
        # plt.plot(PPG["power_spectrum"]["freq"], PPG["power_spectrum"]["density"], c='r')
        # plt.title("Power Spectrum of Overall Data PPG")
        # plt.xlabel("Frequency")
        # plt.ylabel("density")
        # # plt.show()
        # plt.savefig(storage_location_PPG + filename + "_PPG_spectral(overall).png")
        # plt.gcf().clear()

        # Overall compressed data
        # ppg_x = compress_data(ppg_x, compression_factor_PPG)
        # ppg_y = compress_data(ppg_y, compression_factor_PPG)
        # sampling_rate_PPG = sampling_rate_PPG/compression_factor_PPG

        PPG_one, PPG_half = ([] for x in range(2))

        time_unit = 60
        # One Minute Data for PPG
        start_pointer = 0
        end_pointer = time_unit * sampling_rate_PPG
        print("One Minute Data:")
        # print("")
        while start_pointer < len(ppg_x):
            print(start_pointer/sampling_rate_PPG, "-", min(end_pointer/sampling_rate_PPG, len(ppg_x)/sampling_rate_PPG), "seconds")
            temp_ppg = get_details_PPG(ppg_y[int(start_pointer):int(min(end_pointer, len(ppg_x)))], ppg_x[int(start_pointer):int(min(end_pointer, len(ppg_x)))], sampling_rate_PPG, display_details=True, display_graph=False, graph_title=str(start_pointer/sampling_rate_PPG) + "-" + str(end_pointer/sampling_rate_PPG) + "seconds", store_graph=False, store_location=storage_location_PPG + filename + "_(" + str(start_pointer) + "-" + str(end_pointer) + ").png")
            PPG_one.append(temp_ppg)
            start_pointer = end_pointer
            end_pointer = end_pointer + time_unit * sampling_rate_PPG
            # print("")
            # break
            print("")

        time_unit = 30
        # Half Minute Data for PPG
        start_pointer = 0
        end_pointer = time_unit * sampling_rate_PPG
        print("Half Minute Data:")
        # print("")
        while start_pointer < len(ppg_x):
            print(start_pointer/sampling_rate_PPG, "-", min(end_pointer/sampling_rate_PPG, len(ppg_x)/sampling_rate_PPG), "seconds")
            temp_ppg = get_details_PPG(ppg_y[int(start_pointer):int(min(end_pointer, len(ppg_x)))], ppg_x[int(start_pointer):int(min(end_pointer, len(ppg_x)))], sampling_rate_PPG, display_details=True, display_graph=False, graph_title=str(start_pointer/sampling_rate_PPG) + "-" + str(end_pointer/sampling_rate_PPG) + "seconds", store_graph=False, store_location=storage_location_PPG + filename + "_(" + str(start_pointer) + "-" + str(end_pointer) + ").png")
            PPG_half.append(temp_ppg)
            start_pointer = end_pointer
            end_pointer = end_pointer + time_unit * sampling_rate_PPG
            # print("")
            # break
            print("")
    else:
        '''
        splitting instructions are provided, split into given question and its range
        '''
        questions = ((first_line.split(":")[1]).strip()).split(",")
        # print("questions:", questions)
        time_range = list(map(float, ((second_line.split(":")[1]).strip()).split(",")))
        # print("time_range:", time_range)
        question_range = [min(time_range), max(time_range)]
        # print("question_range:", question_range)
        drawing = False
        if questions[-1] == "drawing":
            del questions[-1]
            drawing = True

        try:
            f = open(storage_location_EDA + filename + "_analysis_EDA" + extension, 'w+')
            sys.stdout = f
        except:
            print("Error in opening", storage_location_EDA + filename + "_analysis_EDA" + extension)
            print("Displaying on terminal")

        # Details
        print("Details:", "".join(details))

        # EDA Channel Data
        print("EDA Channel Data")
        # print("")

        if drawing:
            # Overall Data with Drawing
            print("Overall Data(for questions):")
            # print("")
            EDA_questions = get_details_EDA(eda_y[int(math.floor(question_range[0]*sampling_rate_EDA)):int(math.floor(question_range[1]*sampling_rate_EDA+1))], eda_x[int(math.floor(question_range[0]*sampling_rate_EDA)):int(math.floor(question_range[1]*sampling_rate_EDA+1))], sampling_rate_EDA, display_details=True, display_graph=True, graph_title="Overall Data EDA(questions)", store_graph=False, store_location=storage_location_EDA + filename + "_EDA(questions).png")
            # sys.exit()
            print("")
            print("Overall Data(for drawing):")
            # print("")
            EDA_drawing = get_details_EDA(eda_y[int(math.floor(question_range[1]*sampling_rate_EDA+1)):], eda_x[int(math.floor(question_range[1]*sampling_rate_EDA+1)):], sampling_rate_EDA, display_details=True, display_graph=True, graph_title="Overall Data EDA(drawing)", store_graph=False, store_location=storage_location_EDA + filename + "_EDA(drawing).png")
            print("")
        else:
            # Overall Data without Drawing
            print("Overall Data:")
            # print("")
            EDA = get_details_EDA(eda_y, eda_x, sampling_rate_EDA, display_details=True, display_graph=True, graph_title="Overall Data EDA(questions)", store_graph=False, store_location=storage_location_EDA + filename + "_EDA(questions).png")
            print("")
        # sys.exit()

        # Overall compressed data
        # eda_x = compress_data(eda_x, compression_factor_EDA)
        # eda_y = compress_data(eda_y, compression_factor_EDA)
        # sampling_rate_EDA = int(sampling_rate_EDA/compression_factor_EDA)

        EDA_que = []

        # Question-wise data
        print("Question-wise data:")
        # print("")
        for i, val in enumerate(questions):
            start_pointer = min(int(math.floor(float(time_range[2*i])*sampling_rate_EDA)), len(eda_y))
            end_pointer   = min(int(math.floor(float(time_range[2*i + 1])*sampling_rate_EDA)), len(eda_y))
            print(start_pointer/sampling_rate_EDA, "-", min(end_pointer/sampling_rate_EDA, len(eda_x)/sampling_rate_EDA), "seconds", "(question %s)" % val)
            temp_eda = get_details_EDA(eda_y[int(start_pointer):int(min(end_pointer, len(eda_x)))], eda_x[int(start_pointer):int(min(end_pointer, len(eda_x)))], sampling_rate_EDA, display_details=True, display_graph=False, graph_title=str(start_pointer/sampling_rate_EDA) + "-" + str(end_pointer/sampling_rate_EDA) + "seconds", store_graph=False, store_location=storage_location_EDA + filename + "_(" + "question" + "-" + str(val) + ").png")
            EDA_que.append(temp_eda)
            print("")

        try:
            f = open(storage_location_PPG + filename + "_analysis_PPG" + extension, 'w+')
            sys.stdout = f
        except:
            print("Error in opening", storage_location_PPG + filename + "_analysis_PPG" + extension)
            print("Displaying on terminal")

        # Details
        print("Details:", "".join(details))

        # PPG Channel Data
        print("PPG Channel Data")
        # print("")

        if drawing:
            # Overall Data with Drawing
            print("Overall Data(for question):")
            # print("")
            PPG_questions = get_details_PPG(ppg_y[int(math.floor(question_range[0]*sampling_rate_PPG)):int(math.floor(question_range[1]*sampling_rate_PPG+1))], ppg_x[int(math.floor(question_range[0]*sampling_rate_PPG)):int(math.floor(question_range[1]*sampling_rate_PPG+1))], sampling_rate_PPG, display_details=True, display_graph=False, graph_title="Overall Data PPG(questions)", store_graph=False, store_location=storage_location_PPG + filename + "_PPG(questions).png")
            print("")
            # Plotting and saving questions spectral power data
            # plt.plot(PPG_questions["power_spectrum"]["freq"], PPG_questions["power_spectrum"]["density"], c='r')
            # plt.title("Power Spectrum of Overall Data PPG")
            # plt.xlabel("Frequency")
            # plt.ylabel("density")
            # # plt.show()
            # plt.savefig(storage_location_PPG + filename + "_PPG_spectral(questions).png")
            # plt.gcf().clear()

            print("Overall Data(for drawing):")
            # print("")
            PPG_drawing = get_details_PPG(ppg_y[int(math.floor(question_range[1]*sampling_rate_PPG+1)):], ppg_x[int(math.floor(question_range[1]*sampling_rate_PPG+1)):], sampling_rate_PPG, display_details=True, display_graph=False, graph_title="Overall Data PPG(drawing)", store_graph=False, store_location=storage_location_PPG + filename + "_PPG(drawing).png")
            print("")
            # Plotting and saving drawing spectral power data
            # plt.plot(PPG_drawing["power_spectrum"]["freq"], PPG_drawing["power_spectrum"]["density"], c='r')
            # plt.title("Power Spectrum of Overall Data PPG")
            # plt.xlabel("Frequency")
            # plt.ylabel("density")
            # # plt.show()
            # plt.savefig(storage_location_PPG + filename + "_PPG_spectral(drawing).png")
            # plt.gcf().clear()
        else:
            # Overall Data without Drawing
            print("Overall Data:")
            # print("")
            PPG = get_details_PPG(ppg_y, ppg_x, sampling_rate_PPG,  display_details=True, display_graph=False, graph_title="Overall Data PPG(questions)", store_graph=False, store_location=storage_location_PPG + filename + "_PPG(questions).png")
            print("")
            # Plotting and saving questions spectral power data
            # plt.plot(PPG["power_spectrum"]["freq"], PPG["power_spectrum"]["density"], c='r')
            # plt.title("Power Spectrum of Overall Data PPG")
            # plt.xlabel("Frequency")
            # plt.ylabel("density")
            # # plt.show()
            # plt.savefig(storage_location_PPG + filename + "_PPG_spectral(overall).png")
            # plt.gcf().clear()

        # sys.exit()

        # Overall compressed data
        # ppg_x = compress_data(ppg_x, compression_factor_PPG)
        # ppg_y = compress_data(ppg_y, compression_factor_PPG)
        # sampling_rate_PPG = int(sampling_rate_PPG/compression_factor_PPG)

        PPG_que = []

        # Question-wise data
        print("Question-wise data:")
        # print("")
        for i, val in enumerate(questions):
            start_pointer = min(int(math.floor(float(time_range[2*i])*sampling_rate_PPG)), len(ppg_y))
            end_pointer   = min(int(math.floor(float(time_range[2*i + 1])*sampling_rate_PPG)), len(ppg_y))
            print(start_pointer/sampling_rate_PPG, "-", min(end_pointer/sampling_rate_PPG, len(ppg_x)/sampling_rate_PPG), "seconds", "(question %s)" % val)
            temp_ppg = get_details_PPG(ppg_y[int(start_pointer):int(min(end_pointer, len(ppg_x)))], ppg_x[int(start_pointer):int(min(end_pointer, len(ppg_x)))], sampling_rate_PPG, display_details=True, display_graph=False, graph_title=str(start_pointer/sampling_rate_PPG) + "-" + str(end_pointer/sampling_rate_PPG) + "seconds", store_graph=False, store_location=storage_location_PPG + filename + "_(" + "question" + "-" + str(val) + ").png")
            PPG_que.append(temp_ppg)
            print("")
else:
    print("No such file exists")
