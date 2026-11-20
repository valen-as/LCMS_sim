from LCMS_sim.protein_digestion import (digest_protein_sequence, digest_protein_collection, compute_sequence_coverage, enzyme_cleavage_patterns)

def test_digest_protein_sequence():
    dummy_seq_1 = 'DUMMYTESTSEQUENCEONE'
    dummy_seq_2 = 'DUMMYTESTSEQUENCETWO'
    
    pep_1 = digest_protein_sequence(dummy_seq_1, 'S', min_pep_len=5, max_pep_len=30)
    pep_2 = digest_protein_sequence(dummy_seq_2, 'S', min_pep_len=5, max_pep_len=30)

    assert pep_1 == ['DUMMYTE', 'EQUENCEONE']
    assert pep_2 == ['DUMMYTE', 'EQUENCETWO']
    #print("All tests passed.")


def test_digest_protein_collection():
    dummy_proteins = {
        "Protein1": "ACDEFGHIJKLMNOPQRSTUVWXYZ",
        "Protein2": "AAABCDEFGHIJKL"
    }

    assert digest_protein_collection(dummy_proteins, enzyme_cleavage_patterns['Trypsin']) == {
        "Protein1": ["ACDEFGHIJK", "LMNOPQR", "STUVWXYZ"],
        "Protein2": ["AAABCDEFGHIJK"]
    }
    #print("All tests passed.")

def test_compute_sequence_coverage():
    dummy_prot_seq = 'THISISATESTFORCOVERAGE'
    dummy_peps = ['THISIS', 'A', 'TEST', 'FOR', 'COVERAGE']

    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 100.0
    #print("All tests passed.")

