import sequenceAnalysis as sa

def main (fileName=None):
    """Takes an imput of sequences and returns the relative codon usage of the sequences."""
    #Declarations
    myReader = sa.FastAreader("testGenome.fa") 
    haloReader = sa.FastAreader("haloVolc1_1-genes.fa")  
    haloNuc = sa.NucParams()
    myNuc = sa.NucParams()
    #Names of the sequences that are being compared so that it could be used in the future for other sequences
    seq1 = "pyroOgun"
    seq2 = "haloVolc1"

    for head, seq in myReader.readFasta(): #Reads through testGenome.fa
        myNuc.addSequence(seq)   
    for head, seq in haloReader.readFasta(): #Reads through haloVolc1_1-genes.fa
        haloNuc.addSequence(seq) 

    #Compares sequence length    
    baseLen = myNuc.nucCount() / 1000000
    haloBaseLen = haloNuc.nucCount() / 1000000
    diff = abs(haloBaseLen - baseLen)
    if haloBaseLen > baseLen:
        print(f"The {seq2} sequence is longer by {diff:.2f} Mb")
    elif haloBaseLen < baseLen:
        print(f"The {seq2} sequence is shorter by {diff:.2f} Mb")
    else:
        print(f"The {seq2} sequence is the same lenght as {seq1} with a sequence length of {baseLen:.2f} Mb\n")
    print(f"{seq1} sequence length = {baseLen:.2f} Mb\n{seq2} sequence length = {haloBaseLen:.2f} Mb\n")

    #Compares GC content
    seq1gcContent = (myNuc.nucComp["G"] + myNuc.nucComp["C"]) / myNuc.nucCount() * 100 #GC content of seq1
    seq2gcContent = (haloNuc.nucComp["G"] + haloNuc.nucComp["C"]) / haloNuc.nucCount() * 100 #GC content of seq2
    compare = abs(seq2gcContent - seq1gcContent)
    if seq2gcContent > seq1gcContent:
        print(f"The GC content of {seq2} is {compare:.1f}% greater than the GC content of {seq1}")
    elif seq2gcContent < seq1gcContent:
        print(f"The GC content of {seq2} is {compare:.1f}% less than the GC content of {seq1}")
    else:
        print(f"The GC content of {seq2} is the same as {seq1} at {seq1gcContent}%")
    print(f"GC content of {seq1} = {seq1gcContent:.1f}%\nGC content of {seq2} = {seq2gcContent:.1f}%\n")

    #Compares amino acid composition and returns the % difference of each amino acid
    seq1aaComp = myNuc.aaComp
    seq2aaComp = haloNuc.aaComp
    print(f"The overall amino acid composition comparison of {seq2} to {seq1} is as followed:")
    for aa in "-ACDEFGHIKLMNPQRSTVWY":
        seq1Freq = seq1aaComp[aa] / (sum(myNuc.aaComp.values()))
        seq2Freq = seq2aaComp[aa] / (sum(haloNuc.aaComp.values()))
        aaDiff = (seq1Freq - seq2Freq) * 100
        if aaDiff < 0: #If the diff is negative
            result = f"{abs(aaDiff):.2f}% higher "
        elif aaDiff > 0: #If the diff is positive
            result = f"{aaDiff:.2f}% lower  "
        else:
            result = "same"
        print(f"{aa} : {seq2}: ({seq2aaComp[aa]:6d}) is {result} than {seq1}: ({seq1aaComp[aa]:6d})")

    #Compares relative codon bias of the 2 genomes
    matchingKeys = [] #Stores the codons
    print(f"\nThe relative codon bias of {seq2} and {seq1} is as followed:")
    for aa in "-ACDEFGHIKLMNPQRSTVWY": #sorts codons in alpha order since it will iterate through a string in that order.
        match = [key for key, value in myNuc.rnaCodonTable.items() if value == aa] #Creates a list of the codons that match to a specific aa.
        match.sort()
        matchingKeys += match #Merges each sorted list 
    for codon in matchingKeys: 
        val = myNuc.codonComp[codon] / myNuc.aaComp[myNuc.rnaCodonTable[codon]]
        val2 = haloNuc.codonComp[codon] / haloNuc.aaComp[haloNuc.rnaCodonTable[codon]]
        codonDiff = (val2-val) * 100 #halo - pyro
        if codonDiff < 0: #If the difference is negative
            result = f"{abs(codonDiff):5.2f}% lower "
        elif codonDiff > 0:
            result = f"{codonDiff:5.2f}% higher"
        else:
            result = "the same %   "
        print(f"{codon} : {haloNuc.rnaCodonTable.get(codon)} {seq2}: {(val2*100):5.1f} ({haloNuc.codonComp[codon]:6d}) is {result} {seq1}: {(val*100):5.1f} ({myNuc.codonComp[codon]:6d})")


if __name__ == "__main__":
    main()