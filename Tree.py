# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:59:34 2017

@author: cathy
"""
import pickle
from cStringIO import StringIO
from Bio import Phylo
import sys


#a class to construct the newick string of the complete binary transmission(genealoty) tree
epifile, write_to=sys.argv[1:3]

class Tree():
    def __init__(self, epi):
        self.epi=epi
    def Run(self):
        for n in range(1,len(self.epi.Incidence)+1):
            ppl=self.epi.Incidence[-n]
            if ppl.infectby !=None:
                if ppl.infector==False:
                    s1='No'+repr(ppl.ID)+':'+repr(ppl.remove_t-ppl.infecttime)
                else:
                    s1=ppl.newickstring
                if bool(ppl==ppl.infectby.infectee[-1])==True:
                    s2='No'+repr(ppl.infectby.ID)+':'+repr(ppl.infectby.remove_t-ppl.infecttime)
                else:
                    s2=ppl.infectby.newickstring
                s3='No'+repr(ppl.infectby.ID)+':'+repr(ppl.branchlength)

                ppl.infectby.newickstring='('+s1+','+s2+')'+s3
            
        return self.epi.Seed.newickstring

def main():
    with open(epifile, 'rb') as epidata:
        epi = pickle.load(epidata)
        tree=Tree(epi)
        newick_string=tree.Run()
        tree = Phylo.read(StringIO(newick_string), "newick")
        Phylo.write(tree,write_to,'newick')

if __name__ == '__main__': main()
        
        