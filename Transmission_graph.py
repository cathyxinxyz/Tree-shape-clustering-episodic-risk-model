# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 15:47:02 2017

@author: cathy
"""
import networkx as nx
import sys
import cPickle as pickle

epifile, pickled_as=sys.argv[1:3]

G=nx.DiGraph()

with open(epifile, 'rb') as epidata:
    epi = pickle.load(epidata)
    for incid in epi.Incidence:
        ## add an edge between infective and infectee on the graph 
        if incid.infectby !=None:
            from_node=str(incid.infectby.ID)
            to_node=str(incid.ID)
            G.add_edge(from_node, to_node, length=incid.branchlength)
            G.node[from_node]['state']=incid.infector_state
            G.node[to_node]['state']=incid.state_at_infection

nx.write_gpickle(G, pickled_as)

