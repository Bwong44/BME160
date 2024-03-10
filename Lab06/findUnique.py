#!/usr/bin/env python3
# Name: Brandon Wong (bwong44)
# Group Members: Tim Lee, Vibha Rao, Penetha Jayakumar

import sequenceAnalysis as sa

########################################################################
# findUnique Class
# 
# 
########################################################################

class findUnique:
    """Finds unique subsequences across the different tRNAs."""
    InstanceList = [] #Holds all the tRNA objects
    
    def __init__(self):
        """Initializes the findUnique class."""
        self.powerset = set() #Holds powerset of the sequence
        self.head = "" #Holds the header line of the sequence
        self.essential = set()  #Holds essential set of sequence
        self.sequence = "" #Initial sequence without alignment characters
        findUnique.InstanceList.append(self) #Appends the tRNA object

    def cleanSeq(self, seq):
        """Cleans the sequence of any alignment characters from initial sequences."""
        cleanedSeq = seq.replace(".","").replace("_","").replace("-","")
        self.sequence = cleanedSeq
        return cleanedSeq
    
    def createPowerset(self, seq):
        """Creates powerset for each sequence."""
        for i in range(len(seq)): #Nested loop to get all subsequences
            for j in range(i+1, len(seq)+1): #Accouts for slicing exclusion
                self.powerset.add(seq[i:j])

    def findUnique(self):
        """Finds the unique subsequences in the 22 tRNAs."""
        unionSet = set()
        copyInstanceList = findUnique.InstanceList.copy()
        copyInstanceList.remove(self)
        for tRNA in copyInstanceList: #Creates union of all other sets
            unionSet = unionSet.union(tRNA.powerset)
        uniqueSet = self.powerset.difference(unionSet) #Differnce of self vs union of all the others
        #Find Essentials
        self.essential = uniqueSet.copy() #Copy set so we don't run into iteration error
        for unique in uniqueSet:
            for i in uniqueSet:
                try:
                    if unique != i and unique in i:
                        self.essential.remove(i)
                except KeyError: #Handles errors where unique is bigger than i
                    continue    

    def alignSeq(self):
        """Aligns the subsequences based on where they are in the sequence."""
        alignList = [] #Holds a list of list, [position,essential]

        for essentialSeq in self.essential: #Finds the position of where each essential segment starts
            pos = self.sequence.find(essentialSeq)
            while pos != -1:
                alignList.append([pos,essentialSeq])
                pos = self.sequence.find(essentialSeq,pos+1)
        alignList.sort()
        for item in alignList: #Multiples . by the position
            print("."*item[0] + item[1])


########################################################################
# Main
# Here is the main program
# 
########################################################################

def main(inCL=None):
    """Runs the findUnique class to find the unique subsequences of each tRNA."""
    reader = sa.FastAreader()
    for head, seq in reader.readFasta(): #Iterates through the fasta file
            tRNA = findUnique()
            tRNA.head = head.strip()
            tRNA.createPowerset(tRNA.cleanSeq(seq))

    sortedHead = sorted(findUnique.InstanceList, key=lambda x: x.head) #Sorts by header
    for object in sortedHead: #Iterates tRNA objects to find uniques
        object.findUnique()
        print(object.head)
        print(object.sequence)
        object.alignSeq()

if __name__ == "__main__":
    main()  
