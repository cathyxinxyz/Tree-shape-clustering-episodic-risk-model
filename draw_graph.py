# -*- coding: utf-8 -*-

import networkx as nx
import sys

pickled=sys.argv[1:]
G=nx.read_gpickle(pickled)
print G.nodes()
colors = {'I1h': 'r','I1l': 'b','I2h': 'm','I2l':'purple'}
states=nx.get_node_attributes(G,'state')
color_map=[colors[states[node]] for node in G]
nx.draw_graphviz(G, node_color=color_map)