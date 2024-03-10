#!/usr/bin/env python3
# Name: Brandon Wong (bwong44)
# Group Members: Tim Lee, Vibha Rao, Penetha Jayakumar

import sys

class OrfFiner:
    """
    Defines a class that reads the sequence and finds the open reading frames.

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
    
    def __init__(self):
        """Initializes the OrfFinder class and creates dictionaries that stores the ORF data."""
        self.startSet = set()
        self.stopSet = set()
        self.orfList = [] #Holds start positions
        self.geneFragList = [] #[reading frame, start, end, length]
        self.minOrf = 100
        self.longestGene = False

    def readFrame(self, seq):
        """Accepts the sequence as a string from fasta and parses it for the ORFs, looking for gene fragments."""
        adjustFrame = 0 #Stores value of which reading frame it is in 'TAG','TGA','TAA
        danglingStartStop = []
        hasStart = 0
        hasCodon = 0
        for orf in range(3): #Iterates for the 3 reading frames
            frameOne = range(adjustFrame,len(seq),3)
            hasStart = 0
            hasStop  = 0
            self.orfList.clear()
            for frame in frameOne: #Iterates through the sequence 3 codons at a time
                if seq[frame:frame+3] in self.startSet: #Checks for start
                    self.orfList.append(frame + 1)
                    hasStart = 1

                elif seq[frame:frame+3] in self.stopSet and hasStart: #Checkts for stop
                    for pos in self.orfList:
                        self.geneFragList.append([adjustFrame + 1,pos,frame+3,frame+3 - pos + 1])
                    hasStop = 1
                    self.orfList.clear()

                elif hasStart == 0 and seq[frame:frame+3] in self.stopSet: #Checks for dangling stop (we find stop but no start)
                    self.geneFragList.append([adjustFrame+1,1,frame+3,frame+3])
                    hasStop = 1
                    self.orfList.clear()

            if hasStart == 1 and hasStop == 0: #Checks for dangling start (we find start but no stop)
                self.geneFragList.append([adjustFrame+1,1,len(seq),len(seq)])
            elif hasStart == 0 and hasStop == 0 and self.orfList == []: #Checks for dangling start and stop 
                self.geneFragList.append([adjustFrame+1,1,len(seq),len(seq)])
            adjustFrame +=1
                            
    def readReverseComp(self, seq):
        """Takes the reverse complement of the sequuence and searches for valid gene fragments"""
        complementDict = {"A":"T","T":"A","G":"C","C":"G",}
        seqReverse = seq[::-1]
        reverseComplement = "".join(complementDict[nucleotide] for nucleotide in seqReverse)
        adjustFrame = 0 #Stores value of which reading frame it is in 'TAG','TGA','TAA
        negativeFrame = 0 #Stores negative frame value
        hasStart = 0
        hasCodon = 0
        for orf in range(3): #Iterates for the 3 reading frames
            frameOne = range(adjustFrame,len(reverseComplement),3)
            hasStart = 0
            hasStop  = 0
            self.orfList.clear()
            for frame in frameOne: #Iterates through the sequence 3 codons at a time
                if reverseComplement[frame:frame+3] in self.startSet: #Checks for start
                    self.orfList.append(frame + 1)
                    hasStart = 1

                elif reverseComplement[frame:frame+3] in self.stopSet and hasStart: #Checkts for stop
                    for pos in self.orfList:
                        self.geneFragList.append([negativeFrame-1,len(reverseComplement)-pos+1,len(reverseComplement)-frame-2,frame+3 - pos + 1])
                    hasStop = 1
                    self.orfList.clear()

                elif hasStart == 0 and reverseComplement[frame:frame+3] in self.stopSet: #Checks for dangling stop (we find stop but no start)
                    self.geneFragList.append([negativeFrame-1,len(reverseComplement),len(reverseComplement)-frame-2,frame+3])
                    hasStop = 1
                    self.orfList.clear()

            if hasStart == 1 and hasStop == 0: #Checks for dangling start (we find start but no stop)
                self.geneFragList.append([negativeFrame-1,1,len(reverseComplement),len(reverseComplement)])
            elif hasStart == 0 and hasStop == 0 and self.orfList == []: #Checks for dangling start and stop 
                self.geneFragList.append([negativeFrame-1,1,len(reverseComplement),len(reverseComplement)])
            adjustFrame +=1
            negativeFrame -= 1
        
    def minimumOrf(self):
        "Excludes the genes with the length < minOrf"
        minLenList = []
        for gene in range(len(self.geneFragList)): #Iterates through the list of gene fragments
            if self.geneFragList[gene][3] > self.minOrf: #Checks if the gene fragment is greater than the min
                minLenList.append(self.geneFragList[gene])
        self.geneFragList = minLenList

    def sortBySize(self):
        """Takes the list of gene fragments and sorts by the largest gene size then starting number."""
        sortedGeneList = sorted(self.geneFragList, key = lambda x: (-x[-1],x[1])) #Sorts by gene length
        self.geneFragList = sortedGeneList

    def largestOrf(self):
        "Returns only the largest genes if the option is set."
        if self.longestGene == True: #Checks if longest gene is written in terminal
            largestList = []
            compareSet = set()
            for gene in self.geneFragList: #Compares the ending position since it is already sorted by largest gene
                if gene[2] not in compareSet:
                    largestList.append(gene)
                    compareSet.add(gene[2])
            self.geneFragList = largestList
        
    def reverseOrder(self):
        """Reverses the order so that the start and stop match the top strand."""
        for gene in self.geneFragList:
            if gene[0] < 0:
                startHolder = gene[1]
                stopHolder = gene[2]
                gene[1] = stopHolder
                gene[2] = startHolder
            
    def clearGeneFragList(self):
        """Clears out the geneFragList object so that it doesn't store the data of previous sequences if there are any."""
        self.geneFragList = []


class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ["ATG"],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


class NucParams:
    """Defines a class that is used to read sequences and return information about genetic composition."""
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        """Constructs dictionaries to handle an arbitrarily large number of sequences."""
        #Creates a dict for each composition 
        self.aaComp = {aa:0 for aa in self.rnaCodonTable.values()}
        self.nucComp = {nucc:0 for nucc in "ACGTNU"}
        self.codonComp = {codon:0 for codon in self.rnaCodonTable.keys()}
        self.addSequence(inString)

    def addSequence (self, inSeq):
        """Accepets new sequences and adds them to the their necessary dictionaries."""
        seq = inSeq.replace(" ","").upper() #Cleans up the input to handle lower case characters and white spaces
        
        #nucleotide dict
        for nucc in seq:
            if nucc in "ACGTUN": #Searchs for valid necleotides
                self.nucComp[nucc] = self.nucComp.get(nucc) + 1
        
        #codon dict and amino acid dict
        rnaSeq = seq.replace("T","U")
        codonRange = range(0,len(rnaSeq),3)
        for codon in codonRange: #codon is used for indexing
            if rnaSeq[codon:codon+3] in self.codonComp: #Checks if the slice is a valid key and if it is it will add its count
                self.codonComp[rnaSeq[codon:codon+3]] = self.codonComp.get(rnaSeq[codon:codon+3]) + 1
                if rnaSeq[codon:codon+3] in self.rnaCodonTable.keys(): #Decodes the codon to aa and if it's valid it will add to the count
                    self.aaComp[self.rnaCodonTable[rnaSeq[codon:codon+3]]] = self.aaComp.get(self.rnaCodonTable[rnaSeq[codon:codon+3]]) + 1


    def aaComposition(self):
        """Returns the dictionary of the count of amino acids in the sequences."""
        return self.aaComp
    
    def nucComposition(self):
        """Returns the dictionary of the count of nucleotides in the sequences."""
        return self.nucComp
    
    def codonComposition(self):
        """Returns the dictionary of the count of codons in the sequences."""
        return self.codonComp
    
    def nucCount(self):
        "Returns the count of nucleotides in the sequences."
        return sum(self.nucComp.values())
    

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence


class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O
    """ Ex input:  VLSPADKTNVKAAW

    Ex Output:
    Number of Amino Acids: 14
    Molecular Weight: 1499.7
    molar Extinction coefficient: 5500.00
    mass Extinction coefficient: 3.67
    Theoretical pI: 9.88
    Amino acid composition:
    A = 21.43%
    C = 0.00%
    D = 7.14%
    E = 0.00%
    F = 0.00%
    G = 0.00%
    H = 0.00%
    I = 0.00%
    K = 14.29%
    L = 7.14%
    M = 0.00%
    N = 7.14%
    P = 7.14%
    Q = 0.00%
    R = 0.00%
    S = 7.14%
    T = 7.14%
    V = 14.29%
    W = 7.14%
    Y = 0.00%
    """

    def __init__ (self, protein):
        """Constructs an amino acid composition dictionary, removing anything that isn't a valid AA."""
        #Declarations
        self.protein = protein.upper()
        self.aaCompositionDict = {aa:0 for aa in "ACDEFGHIKLMNPQRSTVWY"} #Constructs a dictionary of all valid amino acids.
        self.count = 0 #Holds the count of valid amino acids that can be used by the dictionary and the aaCount method.
            
        for aa in self.protein: #Iterates through the amino acid sequence, looking for valid amino acids.
            if aa in "ACDEFGHIKLMNPQRSTVWY":
                self.count += 1
                self.aaCompositionDict[aa] = self.aaCompositionDict.get(aa) + 1 #Adds to value each time it goes over one of the amino acids in the input.
        
    def aaCount (self):
        """Returns the count of valid amino acids in the """
        return self.count

    def pI (self, precision = 2): #pH that is closest to a charge of 0.
        """Returns the theoretical isolelectric point by using binary search over the range of pHs."""
        #Declarations
        low = 0 #Lowest boundary 
        high = 14 #Highest boundary
        p = 1 / (10 ** precision) #Calculates the precision 

        #Binary Search over the pH range.
        while p <= high - low: 
            mid = (high + low) / 2
            thisCharge = self._charge_(mid)
            if thisCharge > 0: #When the given charge is too low it moves the lowest pH boundry to "mid" and searchs through the new range.
                low = mid
            elif thisCharge < 0: #When the given charge is too high it moves the highest pH boundry to "mid" and searchs through the new range.
                high = mid 
        return mid
            
    def aaComposition (self) :
        """Returns the full dictionary of the amino acids and their composition in the protein."""
        return self.aaCompositionDict

    def _charge_ (self, pH):
        """Returns the net charge of the protien for the given pH"""
        #Declarations
        aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
        aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
        aaNterm = 9.69
        aaCterm = 2.34
        sumPos = 0 #stores the sum of the positive charges
        sumNeg = 0 #stores the sum of the negative charges
        
        for aa in self.aaCompositionDict:
            if aa in "KRH" and self.aaCompositionDict[aa] != 0: #Positive charge conditional
                sumPos += self.aaCompositionDict[aa] * (10 ** aa2chargePos[aa]) / (10 ** aa2chargePos[aa] + 10 ** pH)
            elif aa in "DECY" and self.aaCompositionDict[aa] != 0: #Negative charge conditional
                sumNeg += self.aaCompositionDict[aa] * (10 ** pH) / (10 ** aa2chargeNeg[aa] + 10 ** pH)
        sumPos += (10 ** aaNterm) / (10 ** aaNterm + 10 ** pH)
        sumNeg += (10 ** pH) / (10 ** aaCterm + 10 ** pH)
        netCharge = sumPos - sumNeg
        return netCharge
    
    def molarExtinction (self, cystine=True):
        """Returns the molar extinction coefficient with cystine as an optional parameter.
        
        cystine = True: oxidizing conditions 
        cystine = False: reducing conditions 
        """
        #Declarations
        aa2abs280= {'Y':1490, 'W': 5500, 'C': 125} #absorbance
        yCoefficient = self.aaCompositionDict["Y"] * aa2abs280["Y"]
        WCoefficient = self.aaCompositionDict["W"] * aa2abs280["W"]
        cCoefficient = self.aaCompositionDict["C"] * aa2abs280["C"]
        
        #Condition flow for calculating the coefficient depending on oxidizing and reducing conditions. 
        if cystine == True:
            extinctionCoefficient = yCoefficient + WCoefficient + cCoefficient
        else:
            extinctionCoefficient = yCoefficient + WCoefficient
        return extinctionCoefficient

    def massExtinction (self, cystine=True):
        """Returns the mass extinction coefficient as a single number."""
        myMW =  self.molecularWeight()
        return self.molarExtinction(cystine=True) / myMW if myMW else 0.0

    def molecularWeight (self):
        """Returns the molecular weight of the protein by comparing the dictionary to the aa2mw dictionary."""
        weight = 0
        mwH2O = 18.015
        aa2mw = { #molecular weights dict
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
        
        #Iterates through the aa composition and compares them to mw dict to get the mw of each aa. 
        for aa in self.aaCompositionDict.keys():
            if aa in aa2mw.keys():
                weight += ((aa2mw.get(aa) - mwH2O) * self.aaCompositionDict[aa])
        weight += mwH2O
        return weight

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()