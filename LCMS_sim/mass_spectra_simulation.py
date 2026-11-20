amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    This function calculates the molar mass of a peptide.

    Parameters
    ----------
    peptide_seq: str
        Sequence of the peptide
    amino_acid_mass_dict: dict
        Dictionary that assigns a mass to each aminoacid

    Returns
    -------
    dict
        Dictionary with peptide as key and mass as value
    """
    mass_dict = {}
    for peptide in peptide_seq:
        total_mass = 0
        for aa in peptide_seq:
            if aa in amino_acid_mass_dict:
                total_mass += amino_acid_mass_dict[aa]
            else:
                raise ValueError(f"Unknown amino acid: {aa}")
        mass_dict[peptide_seq] = total_mass

    return mass_dict

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    This function calculates the molar mass of a list of peptides.
    
    Parameters
    ----------
    peptides: list
        list of peptides
    amino_acid_mass_dict: dict
        dictionary that assigns a mass to each aminoacid
    Returns
    -------
    dict
        dictionary with peptide as key and mass as value
    """
    mass_dict = {}
    total_mass = []
    
    for peptide in peptides:
        pep_mass = 0
        for aa in peptide: 
            if aa in amino_acid_mass_dict:
                pep_mass += amino_acid_mass_dict[aa]
            else:
                raise ValueError(f"Unknown amino acid: {aa}")
        total_mass.append(pep_mass)       
        
    mass_dict = dict(zip(peptides, total_mass))
    
    return mass_dict

def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    This function takes a peptide-mass dictionary and gives another dictionary mapping peptides to their m/z value.

    Parameters
    ----------
    peptide_mass_map: dict
        peptide-mass dictionary
    charge: int
        charge (default is 2)
    proton_mass: float
        mass of a proton (default is 1.007)

    Returns
    -------
    dict
        dictionary mapping peptides to their m/z value

    """
    mz_dict={}
    mz_mass=[]

    for mass in peptide_mass_map.values():
        mz_value = 0
        mz_value = (mass + charge * proton_mass) / charge
        mz_mass.append(mz_value)
    mz_dict=dict(zip(peptide_mass_map.keys(), mz_mass))
    return mz_dict

import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    This function plots a simulated mass spectrum by assigning random intensities to each m/z values.

    Parameters
    ----------
    mz_values: dict
        dictionary with the m/z values of each peptide
    random_count_range: any
        range for the intensisties that are going to be generated (default is (0, 30000))
    seed: int
        seed for the random count range (default is 42)

    Returns
    -------
    bar plot
        plots a simulated mass spectrum by assigning random intensities to each m/z values

    """
                                           
    mz_list = []
    intensity = []
    for mz in mz_values.values():
        mz_list.append(mz)
        intensity.append(np.random.randint(random_count_range[0], random_count_range[1]+1))
    
    plt.bar(mz_list, intensity, width=1.5)
    plt.xlabel('m/z value')
    plt.ylabel('Intensity')
    plt.title('Simulated MS Spectrum')

def fragment_peptide(peptide):
    """
    The function returns the combined list of b- and y-ion sequences for any peptide input.

    Parameters
    ----------
    peptide: str
        the peptide to be fragmented

    Returns
    -------
    list
        list of b- and y-ion fragments

    """
    b_ions = []
    y_ions = []

    for i in range(1, len(peptide)):
        b_ions.append(peptide[:i])
        y_ions.append(peptide[-i:])
    return b_ions + y_ions