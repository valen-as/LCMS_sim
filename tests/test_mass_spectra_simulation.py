from LCMS_sim.mass_spectra_simulation import (amino_acid_mass_dalton, calculate_mol_mass, calculate_mol_mass_collection, calculate_mz_collection, plot_spectrum, fragment_peptide)

def test_calculate_mol_mass():
    aa_mass_dict = {
        "A" : 1.0,
        "B" : 2.0,
        "C" : 3.0
    }
    test_peptide = 'AB'
    expected_mass = {'AB':3.0}
    calculated_mass = calculate_mol_mass(test_peptide, aa_mass_dict)
    assert calculated_mass == expected_mass

def test_calculate_mol_mass_collection():
    aa_mass_dict = {
        'A':10,
        'B':20,
        'C':30
    }
    peptides = ['AAA', 'BBB', 'CCC', 'ABC']
    expected = {'AAA':30, 'BBB':60, 'CCC':90, 'ABC':60}
    actual = calculate_mol_mass_collection(peptides,aa_mass_dict)

    assert actual == expected

def test_calculate_mz_collection():
    peptide_mass_map = {
        'pep1':100,
        'pep2':200,
        'pep3':300
    }
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {'pep1': (100 + 2 * 1.007) / 2,
                'pep2': (200 + 2 * 1.007) / 2,
                'pep3': (300 + 2 * 1.007) / 2
                }

    assert actual == expected

def test_fragment_peptide():
    peptide = 'ABCD'
    expected = ['A', 'AB', 'ABC', 'D', 'CD', 'BCD']
    actual = fragment_peptide(peptide)

    assert set(actual) == set(expected)
    assert 'ABCD' not in expected