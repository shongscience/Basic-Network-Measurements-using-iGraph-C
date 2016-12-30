# Basic-Network-Measurements-using-iGraph-C
1. Some basic (and ad-hoc) codes to get network quantities using iGraph/C/C++

2. iGraph is available at http://igraph.org/

3. Once you have installed iGraph/C properly, you can compile 
  - netscalar.c : print basic network quantities from input "edge" files
  - netprintmembership.c : print a membership id for each vertex 
 
    ** network example files : ./misc/ simple.edge, simpleplus.edge, ring20.edge, and star20.edge

  3.1 netscalar.c
    - usage : netscalar.bin in.edge out.network
    - ex) netscalar.bin ring20.edge dummyresult.network
      
  3.2 netprintmembership.c
    - usage : netprintmembership.bin in.edge out.network
    - ex) netprintmembership.bin simpleplus.edge dummyresult.network
    
    ** These are simple modifed suites from igraph example files. Most jobs are likely done by shell- and python-scripts utilizing these simple c-binaries. As an example, you can check "generateScripts.py".
  
    ** When the numbers of vertices and edges are not large, you'd better use igraph/python in an ipython environment.

4. Misc folder
  - Some utility files related to astrometry, network edges, etc.
