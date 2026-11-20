from pyteomics import achrom 

def predict_lc_retention_times(peptides):
    """
    This function predicts the liquid chromatography retention times for a list of peptides.

    Parameters
    ----------
    peptides : list
        A list of peptide sequences.

    Returns
    -------
    list
        The predicted relative retention times of each peptide with only two decimal places.
    """
    retention_map = {}
    for peptide in peptides:
        relative_RT = achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0)
        retention_map[peptide] = round(relative_RT, 2)
    return retention_map


import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    This function plots the distribution of retention times.

    Parameters
    ----------
    retention_times : list
        A list of retention times.
    resolution : int, optional
        The number of bins for the histogram (default is 30).

    Returns
    -------
    Histogram
        The function displays the plot.
    """
    plt.hist(retention_times, bins=resolution) 
    plt.xlabel('Retention time (minutes)')
    plt.ylabel('Intensity')
    plt.title('Histogram of Predicted LC Retention Times')
    plt.show()


def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    This function selects peptides within a specified retention time window.

    Parameters
    ----------
    peptide_rt_map : dict
        A dictionary mapping peptide names to their retention times.
    lower_ret_time : float
        The lower limit of the retention time window.
    upper_ret_time : float
        The upper limit of the retention time window.

    Returns
    -------
    list
        A list of peptides that fall within the specified retention time window.
    """
    peptides_in_window = [
        peptides 
        for peptides, rt in peptide_rt_map.items() 
        if lower_ret_time <= rt <= upper_ret_time
        ]

    return peptides_in_window

