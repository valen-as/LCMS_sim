from LCMS_sim.liquid_chromatography import (predict_lc_retention_times, plot_retention_time, select_retention_time_window)

def test_predict_lc_retention_times():
    peptides = ['PEPTIDE', 'AAA', 'ACDEFG', 'TEST']
    expected = {'PEPTIDE': 7.80, 'AAA': 6.60, 'ACDEFG': 9.70, 'TEST': -1.20}

    actual = predict_lc_retention_times(peptides)
    assert actual == expected

def test_select_retention_time_window():
    peptide_rt_map = {
        "peptide1": 8.5,
        "peptide2": 15.2,
        "peptide3": 20.8
    }
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=8.5, upper_ret_time=15.2)

    assert selected == ["peptide1", "peptide2"]