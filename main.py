# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")
    
    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    #BLOSUM matrix is a substitution matrix containing precomputed scores that indicate how likely one aa is to be substituted by another due to evolution
    #assigns positive scores for similar aa, and neg for dissimilar subst. 
    #if a pair seqA[i-1], seqB[j-1] is found is BLOSUM, it is retreived, else a default mismatch score of -1 is used
    nw = NeedlemanWunsch("BLOSUM.txt", gap_open= -10, gap_extend= -1)

    species = {
        #header is used as a dictionary key, will see e.g. Gallus_gallus: 95.0 for the NW score
        gg_header: nw.align(hs_seq, gg_seq)[0],
        mm_header: nw.align(hs_seq, mm_seq)[0],
        br_header: nw.align(hs_seq, br_seq)[0],
        tt_header: nw.align(hs_seq, tt_seq)[0]

    }
    
    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    #use sorted() to sort species by similarity to human BRD2
    #species.item() returns a list of tuples, so gg_header, score_gg. using the .items method you can convert a dict into a list of tuples.
    #the key=lambda x: x[1] tells sorted() to sort based on the second item in each tuple (the alignment score)
    #reverse=True ensures that the higher score appears first
    #sorted() is sorting the dictionary items. x represents each item in species.items(), x[1] extracts the second element of each tuple, which is the align score
    #the key parameter defines the sorting criteria. without it, python would sort based on first element, so the species name. but we want to sort by alignment score
    sorted_species = sorted(species.items(), key = lambda x: x[1], reverse=True)
    for name, score in sorted_species:
        print(f'{name}: {score}')


if __name__ == "__main__":
    main()
