# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 11:29:30 2016

@author: shong
"""

import networkx as nx
import pygraphviz as pgv
import matplotlib.pyplot as plt

plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman') 
plt.rcParams.update({'font.size': 18})
plt.figure(figsize=(10,5))

c3 = nx.complete_graph(3)
c4 = nx.complete_graph(4)
c5 = nx.complete_graph(5)

###
#A = pgv.AGraph(strict=True,directed=False)
#A = nx.to_agraph(G)

#A.node_attr['shape']='circle'

#A.draw('gschem.png',prog='dot')
#A.write('6clique.dot')
####

fig = plt.figure(figsize=(6,2))

plt.subplot(131)
plt.title('3-clique',fontsize=18)
pos=nx.graphviz_layout(c3,prog="circo",root=0)
nx.draw(c3,pos,node_size=100,alpha=0.5,node_color="blue", with_labels=False)
plt.axis('equal')

plt.subplot(132)
plt.title('4-clique',fontsize=18)
pos=nx.graphviz_layout(c4,prog="circo",root=0)
nx.draw(c4,pos,node_size=100,alpha=0.5,node_color="blue", with_labels=False)
plt.axis('equal')

plt.subplot(133)
plt.title('5-clique',fontsize=18)
pos=nx.graphviz_layout(c5,prog="circo",root=0)
nx.draw(c5,pos,node_size=100,alpha=0.5,node_color="blue", with_labels=False)
plt.axis('equal')

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig('schemaclique.eps')

plt.show()


s6 = nx.star_graph(6)
r6 = nx.cycle_graph(7)
c7 = nx.complete_graph(7)

fig = plt.figure(figsize=(9,3))

labels={}
labels[0]='2'
labels[1]='2'
labels[2]='2'
labels[3]='2'
labels[4]='2'
labels[5]='2'
labels[6]='2'

plt.subplot(131)
plt.title('Ring',fontsize=16)
pos=nx.graphviz_layout(r6,prog="circo",root=0)
nx.draw(r6,pos,node_size=500,node_color="white")
nx.draw_networkx_labels(r6,pos,labels,font_size=16)
plt.text(25,-85,'Centralization = 0',fontsize=16)
plt.text(25,-130,'Transitivity     = 0',fontsize=16)
plt.axis('equal')

print pos

labels={}
labels[0]='6'
labels[1]='1'
labels[2]='1'
labels[3]='1'
labels[4]='1'
labels[5]='1'
labels[6]='1'

plt.subplot(132)
plt.title('Star',fontsize=16)
pos=nx.graphviz_layout(s6,prog="twopi",root=0)
nx.draw(s6,pos,node_size=500,node_color="white")
nx.draw_networkx_labels(s6,pos,labels,font_size=16)
plt.text(30,-15,'Centralization = 1',fontsize=16)
plt.text(30,-35,'Transitivity     = 0',fontsize=16)
plt.axis('equal')

print pos

labels={}
labels[0]='6'
labels[1]='6'
labels[2]='6'
labels[3]='6'
labels[4]='6'
labels[5]='6'
labels[6]='6'

plt.subplot(133)
plt.title('Clique',fontsize=16)
pos=nx.graphviz_layout(c7,prog="twopi",root=0)
nx.draw(c7,pos,node_size=500,node_color="white")
nx.draw_networkx_labels(c7,pos,labels,font_size=16)
plt.text(10,-15,'Centralization = 0',fontsize=16)
plt.text(10,-35,'Transitivity     = 1',fontsize=16)
plt.axis('equal')

plt.tight_layout(pad=0.8, w_pad=0.0, h_pad=0.0)
plt.savefig('schemacentralization.eps')

print pos

plt.show()


# make a triangle for demonstrating transitivity
tr = nx.Graph()
tr.add_nodes_from([1,2,3])
tr.add_edges_from([(1,2),(2,3),(3,1)])



fig = plt.figure(figsize=(5,5))

plt.subplot(211)
plt.title('Transitivity',fontsize=18)
plt.axis('off')
pos=nx.graphviz_layout(tr,prog="circo",root=0)
#nx.draw(tr,pos,node_size=100,alpha=0.5,node_color="blue", with_labels=False)
nx.draw_networkx_nodes(tr,pos,node_size=100,alpha=0.5,node_color="blue", with_labels=False)
nx.draw_networkx_edges(tr,pos,edgelist=[(1,2),(1,3)],edge_color='k')
nx.draw_networkx_edges(tr,pos,edgelist=[(2,3)],edge_color='k',style='dotted')
plt.text(70,100,'?',fontsize=18)
#plt.text(150,70,r'$\frac{a}{b}$',fontsize=18)
plt.axis('equal')


# make a triangle for demonstrating transitivity
bc = nx.Graph()
bc.add_nodes_from([0,1,2,3,4,5])
bc.add_edges_from([(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(2,5)])

pos=nx.graphviz_layout(bc,prog="circo")
pos[0]=(-150,100)
pos[1]=(0,20.0)
pos[2]=(200,20.0)
pos[3]=(350,100)
pos[4]=(200,180)
pos[5]=(0,180)

labels={}
labels[0]='i'
labels[3]='j'
labels[1]=''
labels[2]=''
labels[5]=''
labels[4]=''


plt.subplot(212)
plt.title('Diameter and Betweenness Centrality',fontsize=18)
plt.axis('off')
plt.axis([-250,450,-50,250])
#nx.draw(tr,pos,node_size=100,alpha=0.5,node_color="blue", with_labels=False)
nx.draw_networkx_nodes(bc,pos,nodelist=[0,3],node_size=500,node_color="white",with_labels=True)
nx.draw_networkx_nodes(bc,pos,nodelist=[2,5],node_size=100,alpha=0.5,node_color="red", with_labels=False)
nx.draw_networkx_nodes(bc,pos,nodelist=[1,4],node_size=100,alpha=0.5,node_color="blue", with_labels=False)
nx.draw_networkx_labels(bc,pos,labels,font_size=16)
nx.draw_networkx_edges(bc,pos,edgelist=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(2,5)],edge_color='k')
plt.text(-50,200,r'$\frac{1}{3} + \frac{1}{3}$',fontsize=20)
plt.text(150,-25,r'$\frac{1}{3} + \frac{1}{3}$',fontsize=20)
plt.text(225,200,r'$\frac{1}{3}$',fontsize=20)
plt.text(-40,-25,r'$\frac{1}{3}$',fontsize=20)
#plt.axis('equal')



plt.tight_layout(pad=0.8, w_pad=0.0, h_pad=0.0)
plt.savefig('schematrbc.eps')
plt.show()

















