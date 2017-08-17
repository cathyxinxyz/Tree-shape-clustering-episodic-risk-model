# -*- coding: utf-8 -*-
import dendropy
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy import stats
import numpy as np

       
def High_acute_cluster(tree,cutoffs):
    clusters=dict()
    high_acute_transmission=dict()
    high_to_high=dict()
    cluster_number=dict()
    cluster_size=dict()
    for cf in cutoffs:
        clusters[cf]=list()
    node_not_labeled=['I' in n.label for n in tree.internal_nodes()].count(False)

    if node_not_labeled<5:
        for n in tree.internal_nodes():
            for cf in cutoffs:            
                #print n, n.child_nodes(),n.distance_from_root(),len(n.leaf_nodes())
                if n.distance_from_root()>cf and n.parent_node !=None and n.parent_node.distance_from_root()<cf:
                    clusters[cf].append((n.label, len(n.leaf_nodes())))

        for cf in cutoffs:
            #print [l[-7:-4] for l in clusters[cf]].count('I1h')
            if clusters[cf]!=[]:
                high_acute_transmission[cf]=float([l[0][-7:-4] for l in clusters[cf]].count('I1h'))/len(clusters[cf])
                high_to_high[cf]=float([l[0][-7:] for l in clusters[cf]].count('I1h I1h'))/len(clusters[cf])
                cluster_number[cf]=len(clusters[cf])
                cluster_size[cf]=np.mean([l[1] for l in clusters[cf]])
                #acute_transmission[cf]=float([l[-7:-5] for l in clusters[cf]].count('I1'))/len(clusters[cf])
        #print high_acute_transmission, high_to_high
        return high_acute_transmission, high_to_high, cluster_number,cluster_size
    else:
        return None, None, None,None


####function that read in tree file simulated using treesim function and return the fraction of tree branching events that are caused by transmissions 
#from high risk acute HIV infection to high risk susceptible individuals, frac_H1H
#from high risk acute HIV infection, frac_H1
#from acute HIV infection, frac_1
def Fraction_Acute_transmission(tree):
    node_not_labeled=['I' in n.label for n in tree.internal_nodes()].count(False)
    #print "{} nodes are not labeled".format(node_not_labeled)
    if node_not_labeled<5:
        labels=list()
        for n in tree.internal_nodes():
            labels.append(n.label[-7:])
        frac_H1H=float(labels.count('I1h I1h'))/len(labels)
        frac_H1=(float(labels.count('I1h I1h'))+labels.count('I1h I1l'))/len(labels)
        frac_1=(float(labels.count('I1h I1h'))+labels.count('I1h I1l')+labels.count('I1l I1h')+labels.count('I1l I1l'))/len(labels)
                
        return frac_H1H, 
    else:
        return None,None,None


        
