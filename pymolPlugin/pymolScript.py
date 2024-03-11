#run E:\pymol\pymolPlugin\pymolScript.py
from pymol import cmd

def saveFasta():
    """Saves the sequence in FASTA format, redirecting to a file."""
    sequence = cmd.get_fastastr(selection="all")
    with open('test.txt', 'w') as f: #Saves to same directory as where the code is 
        f.write(sequence)



def main():
    """Runs the code in this file."""
    saveFasta()

main()