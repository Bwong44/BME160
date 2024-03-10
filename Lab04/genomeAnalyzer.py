import sequenceAnalysis as sa

def main (fileName=None):
    """Takes an imput of sequences and returns the relative codon usage of the sequences."""
    myReader = sa.FastAreader(fileName) 
    myNuc = sa.NucParams()

    #Reads the given fasta file
    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)   

    #Prints sequence length
    baseLen = myNuc.nucCount() / 1000000
    print(f"sequence length = {baseLen:.2f} Mb\n")

    #Prints GC content by adding Gs and Cs and dividing them by the total count of nuc.
    gcContent = (myNuc.nucComp["G"] + myNuc.nucComp["C"]) / myNuc.nucCount() * 100
    print(f"GC content = {gcContent:.1f}%\n")

    #Prints the output of codon freq.
    matchingKeys = [] #Stores the codons
    for aa in "-ACDEFGHIKLMNPQRSTVWY": #sorts codons in alpha order since it will iterate through a string in that order.
        match = [key for key, value in myNuc.rnaCodonTable.items() if value == aa] #Creates a list of the codons that match to a specific aa.
        match.sort()
        matchingKeys += match #Merges each sorted list 
    for codon in matchingKeys: 
        val = myNuc.codonComp[codon] / myNuc.aaComp[myNuc.rnaCodonTable[codon]]
        print(f"{codon} : {myNuc.rnaCodonTable.get(codon)} {(val*100):5.1f} ({myNuc.codonComp[codon]:6d})")

if __name__ == "__main__":
    main() # make sure to change this in order to use stdin