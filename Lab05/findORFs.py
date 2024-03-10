#!/usr/bin/env python3
# Name: Brandon Wong (bwong44)
# Group Members: Tim Lee, Vibha Rao, Penetha Jayakumar

"""Finds ORFs in a given sequence given parameters from the command line."""

"""Design Method:
Create a class in sequenceAnalysis to handle all the operations reading the sequence frame by frame.

Iterates 3 times for each frame
    Iterates for each codon stepping by 3 characters
        If has start codon:
            Append to start codon list
        If has stop codon:
            For items in start codon list:
                Appends to geneFragList (frame,start,stop,length)
            Clears start codon list
        If there has been no start codons but found a stop codon we have found a dangling stop:
            Appends to geneFragList(frame,start,stop,length)
            Clears start codon list
    If has start but no stop we have found a dangling start:
        Appends to geneFragList(frame,start,stop,length)
    If doesn't have start, stop we have found a dangling start stop:
        Append to geneFragList(frame,start,stop,length)
    adjustFrame += 1
"""

import sequenceAnalysis

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   
def main(inFile = None, options = None):
    """Finds all valid gene fragments using the terminal to specify options."""
    thisCommandLine = sequenceAnalysis.CommandLine(options)
    reader = sequenceAnalysis.FastAreader(inFile)
    finder = sequenceAnalysis.OrfFiner()
    finder.startSet = set(thisCommandLine.args.start) #Converts to set
    finder.stopSet = set(thisCommandLine.args.stop) #Converts to set
    finder.minOrf = thisCommandLine.args.minGene
    finder.longestGene = thisCommandLine.args.longestGene
    for head, seq in reader.readFasta():
        print(head) #Returns name of sequence
        finder.readFrame(seq)  
        finder.readReverseComp(seq) 
        finder.minimumOrf()
        finder.sortBySize() 
        finder.largestOrf()
        finder.reverseOrder()
        for gene in finder.geneFragList:
            print(f"{gene[0]:+d} {gene[1]:>5d}..{gene[2]:>5d} {gene[3]:>5d}")
        finder.clearGeneFragList()

    
if __name__ == "__main__":
    main()