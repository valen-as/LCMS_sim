import re

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)', 
    'LysN': r'(?=K)', 
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)'
}


def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    This function takes a protein sequence and a cleavage pattern, and returns the digested peptides if they are bigger than 5 and smaller than 30 amino acids.

    Parameters
    ----------
    protein_seq : str
        The protein sequence to be digested.
    cleave_pattern : str
        The cleavage pattern to use for digestion.
    min_pep_len : int, optional
        The minimum length of the peptides to include (default is 5).
    max_pep_len : int, optional
        The maximum length of the peptides to include (default is 30).

    Returns
    -------
    list
        A list of the digested peptides that meet the length criteria.
    """

    peptides = re.split(cleave_pattern, protein_seq)
    
    digested_peptides = []
    for peptide in peptides: 
        if min_pep_len <= len(peptide) <= max_pep_len:
            digested_peptides.append(peptide)
    return(digested_peptides)

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digest a collection of proteins and return a dictionary of the resulting peptides.

    Parameters
    ----------
    protein_map : dict
        A dictionary mapping protein IDs to their sequences.
    cleave_pattern : str
        The cleavage pattern for the protease.
    min_pep_len : int, optional
        The minimum peptide length (default is 5).
    max_pep_len : int, optional
        The maximum peptide length (default is 30).

    Returns
    -------
    dict
        A dictionary mapping protein IDs to lists of digested peptides.

    """
    digested_proteins = {}
    for protein_id, protein_seq in protein_map.items():
        peptides = digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len, max_pep_len)
        digested_proteins[protein_id] = peptides
        
    return digested_proteins

def compute_sequence_coverage(protein_seq, peptides):
    """
    The function gives how much coverage the digested peptides provide for the protein sequence.

    Parameters
    ----------
    protein_seq : str
        The amino acid sequence of the protein.
    peptides : list
        A list of digested peptides.

    Returns
    -------
    float
        The percentage of the protein sequence covered by the peptides.
    """
    
    covered_positions = set()  
    for peptide in peptides:
        start_index = protein_seq.find(peptide) 
        if start_index != -1: 
           covered_positions.update(range(start_index, start_index + len(peptide))) 
    return coverage


