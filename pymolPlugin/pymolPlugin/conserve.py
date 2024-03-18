#!/usr/bin/env python3
#Name: Brandon Wong (bwong44) 
#Group Members: Tim Lee

# This program calculates the conservation of residues across multiple protein sequences within pymol.
# It also visualizes the conservation data.
#
# Dependencies: Install matpoltlib and Bio python packages within the PyMOL environment
#
# Input: The user must have the protein sequences loaded into PyMOL
#
# Output: File containing the protein sequences(fasta), conservation scores(txt file),
#         heatmap of the residue frequencies(png), and alligned sequences(png)

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from Bio import AlignIO
from collections import Counter
from pymol import cmd, stored

####################################################################################
# conserve class
# finds the conservation of residues across multiple protein sequences within pymol
# 
####################################################################################

class conserve:
    """Creates labels for the alpha carbons of all the protein."""
    InstanceList = []

    def __init__(self):
        """Initializes the conserve class."""
        self.header = "" #Holds the header line of the sequence
        self.sequence = "" #Holds the sequence of the protein
        conserve.InstanceList.append(self)

    def outputData(fileName):
        """Outputs the data of the proteins to a file."""
        with open(fileName, "w") as f:
            for instance in conserve.InstanceList:
                f.write(f"{instance.header}\n")
                f.write(f"{instance.sequence}\n")

    def calculateConservation(self, proteins):
        """Calculates the conservation of residues across multiple protein sequences."""
        # Create an alignment from the protein sequences
        alignment = [protein.sequence for protein in proteins]

        #Calculate conservation scores for each position in the alignment
        conservationScores = []
        for i in range(len(alignment[0])): # Count the residues at this position
            residueCount = Counter(seq[i] for seq in alignment)
            #Calculate the conservation score as the frequency of the most common residue
            conservationScore = max(residueCount.values()) / len(alignment)
            conservationScores.append(conservationScore)
        return conservationScores

    def calculateFrequencies(self, sequences):
        """Calculates the frequency of each residue at each position."""
        #Initialize a dictionary to store the frequencies
        frequencies = {residue: [0]*len(sequences[0]) for residue in 'ACDEFGHIKLMNPQRSTVWY'}

        #Calculate the frequencies
        for sequence in sequences:
            for i, residue in enumerate(sequence):
                if residue in frequencies:
                    frequencies[residue][i] += 1
        return frequencies 
    
    def frequenciesToArray(self, frequencies):
        """Converts the frequencies dictionary to a 2D array."""
        residues, freqs = zip(*frequencies.items()) #Get the residues and their frequencies
        freqsArray = np.array(freqs) #Convert the frequencies to a 2D array
        return residues, freqsArray

########################################################################
# class dataVisualization
# visulaizes the conservation data
# 
########################################################################

class dataVisualization:
    """Visualizes the conservation data."""
    def __init__(self):
        """Initializes the dataVisualization class."""

    def plotHeatmap(self, residues, freqsArray):
        """Creates a heatmap of the residue frequencies."""
        cmap = LinearSegmentedColormap.from_list("myCmap", ["white", "red"])
        plt.imshow(freqsArray, aspect="auto", cmap=cmap, interpolation="nearest")
        plt.colorbar(label="Frequency")
        plt.yticks(range(len(residues)), residues)
        plt.xlabel("Residue ID")
        plt.xticks(range(1, freqsArray.shape[1] + 1, 5))
        plt.ylabel("Residue")
        plt.grid(True, which="both", color="black", linewidth=1)
        plt.show()
    
    def plotAlignmentWithScores(self, alignmentFile):
        """Creates a PNG image of a sequence alignment."""
        alignment = AlignIO.read(alignmentFile, "fasta")
        fig, ax = plt.subplots(figsize=(10, len(alignment)*0.5)) #Create a figure and axes

        for i, record in enumerate(alignment):
            aminoAcids = list(str(record.seq)) #AA list
            x = range(len(aminoAcids)) #x positions list
            y = [i*0.5 for _ in aminoAcids] #y positions list so they aren't as spread out
            ax.scatter(x, y, marker="s", s=1000) #Create a scatter plot of the amino acids
            for a, xPos in zip(aminoAcids, x): #Iterate through AA and x positions
                ax.text(xPos, i*0.5, a, ha="center", va="center", color="white", fontsize=10) #Add the amino acid letter as text
        
        for position, xPos in enumerate(x, start=1): #Mark every 10th position in the alignment
            if position % 10 == 0 or position == 1:
                ax.text(xPos, len(alignment)*0.5, str(position), ha="center", va="center", fontsize=10)  # Adjust the font size here
        ax.axis("off")
        plt.show()

####################################################################################
# Main
# Here is the main program
# 
####################################################################################

def main():
    """Runs the createLabels class to create objects for each protein."""
    proteinNames = cmd.get_names("objects")
    proteinObjects = [] #stores projtein objects

    for proteinName in proteinNames:
        if "_" in proteinName: #Checks if the object is a sequence or a grouping of sequences
            continue
        else:
            proteinObject = conserve() #Creates an class object for each protein
            sequenceObject = cmd.get_fastastr(proteinName)
            fastaLines = sequenceObject.split("\n") #Split the sequence into lines
            for line in fastaLines: #Stores the header and sequence of the protein
                if line.startswith(">"):
                    proteinObject.header = line
                else:
                    if proteinObject.sequence:
                        proteinObject.sequence += line
                    else:
                        proteinObject.sequence = line
            proteinObjects.append(proteinObject) #Adds to proteinObjects list
    
    outputFile = "out.fasta" #Sets name of the output file that will store the protein sequences
    conserve.outputData(outputFile)

    conservationScores = proteinObject.calculateConservation(proteinObjects) #Calculate conservation scores for the proteins
    residues = range(1, len(conservationScores) + 1)
    with open("conservationScores.txt", "w") as file: #Write the residue IDs and conservation scores to the file
        for residueId, score in zip(residues, conservationScores):
            file.write(f"Residue ID: {residueId}, Conservation Score: {score}\n")

    #Iterate through the conservation scores and set the B-factors of the residues in the protein
    for proteinName in proteinNames:
        for i, score in enumerate(conservationScores):
            cmd.alter(f"{proteinName} and resi {i+1}", f"b={score}")
        #Get the current B-factors
        stored.b = []
        cmd.iterate(f"{proteinName} and name CA", "stored.b.append(b)")
        minB = min(stored.b) #Find the minimum and maximum B-factor
        maxB = max(stored.b)
        cmd.alter(f"{proteinName} and name CA", f"b=(b-{minB})/{maxB - minB}") #Normalize the B-factors

        #Adds residue to alpha carbon labels
        cmd.select("alpha_carbons", f"{proteinName} and name CA") #Selects all alpha carbons
        stored.alpha_carbons = [] # Initialize an empty list to store alpha carbon identifiers
        #Iterate over the selected alpha carbons and store their identifiers in the list
        cmd.iterate("alpha_carbons", "stored.alpha_carbons.append((model, segi, chain, resi, resn))")
        for ac in stored.alpha_carbons: # Labels all the alpha carbons with their residue information
            #Get the conservation score for this residue
            score = conservationScores[int(ac[3])-1]
            #Add the residue ID and conservation score to the label
            cmd.label(f"{ac[0]} and chain {ac[2]} and resi {ac[3]} and name CA", f"'{ac[4]}\\nResidue id:{ac[3]}\\nConservation Score:{score}\\n'")
        cmd.spectrum("b", "white_red", f"{proteinName} and name CA") #Colorize the protein, white = low conservation, red = high conservation

    sequences = [protein.sequence for protein in proteinObjects]
    frequencies = proteinObjects[0].calculateFrequencies(sequences) #Calculate the frequencies
    residues, freqsArray = proteinObjects[0].frequenciesToArray(frequencies) #Convert the frequencies to a 2D array
    plotter = dataVisualization()
    plotter.plotHeatmap(residues, freqsArray) #Create a heatmap of the frequencies
    plotter2 = dataVisualization()
    plotter2.plotAlignmentWithScores(outputFile)

main()
