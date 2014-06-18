import os,sys
import string
import numpy as np


class LoadPDB:
    """ PDB has the following format:
        HETATM    1  O1  UNK     1    -146.252  54.702  44.478  1.00  0.00           O
        HETATM    2  O2  UNK     1    -147.422  56.482  43.850  1.00  0.00           O
        HETATM    3  O3  UNK     1    -150.638  53.942  46.578  1.00  0.00           O
        HETATM    4  O4  UNK     1    -147.541  55.434  50.459  1.00  0.00           O
         ......
    """
    def __init__(self,fn):
        pdb = open(fn,'r')
        lines = pdb.readlines()
        pdb.close()
        newlines = []
        for line in lines:
            if string.strip(line[0:7])=="HETATM":
                newlines.append(line)
        self.atomid = [string.strip(line[7:12]) for line in newlines]
        self.atomname = [string.strip(line[12:16]) for line in newlines]
        self.resname = [string.strip(line[16:21]) for line in newlines]
        self.resid = [string.strip(line[21:29]) for line in newlines]
        self.x = [string.strip(line[29:38]) for line in newlines]
        self.y = [string.strip(line[38:46]) for line in newlines]
        self.z = [string.strip(line[46:55]) for line in newlines]

    def Dict(self):    
        XYZ=np.array([self.x,self.y,self.z]) #angstrom
        AtomID = np.array(self.atomid,'int')
        AtomNames = np.array(self.atomname,'str') 
        ResidueNames = np.array(self.resname,'str')
        ResidueID = np.array(self.resid,'str')
        Dict = {"XYZ":XYZ,"AtomID":AtomID,"AtomNames":AtomNames,"ResidueNames":ResidueNames,"ResidueID":ResidueID}
        return Dict        

def main():
   fn = '/Users/tud51931/projects/QM/albocycline/albo-cineromycinB/minimized_pdbs_for_QM/output/3.pdb'
   pdb = LoadPDB(fn) 
   dic = pdb.Dict()
   print dic
if __name__=="__main__": main()
