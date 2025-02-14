# Importing Dependencies
import pytest
from align.align import NeedlemanWunsch, read_fasta
import numpy as np

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from align import NeedlemanWunsch, read_fasta

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    #in the read_fasta, the return order is seq, header. so header is our throwaway variable, we only care for the sequence
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open= -10, gap_extend= -1)
    nw.align(seq1, seq2)

    #check if matrix is correctly initialized 
    assert nw._align_matrix.shape == (len(seq1) + 1, len(seq2) + 1)


def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open= -10, gap_extend= -1)
    score, aligned_seqA, aligned_seqB = nw.align(seq3, seq4)

    expected_aligned_seqA = "MAVHQLIRRP"
    expected_aligned_seqB = "M---QLIRHP"

    assert aligned_seqA == expected_aligned_seqA
    assert aligned_seqB == expected_aligned_seqB






