# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 11:03:24 2016

@author: shong
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:16:44 2016

@author: shong
"""

import matplotlib.pyplot as plt
import numpy as np
import os

#from scipy.integrate import simps
from scipy.interpolate import interp1d
from halomod.integrate_corr import AngularCF, angular_corr_gal, flat_z_dist, dxdz
from hmf.cosmo import Cosmology
#from mpmath import gamma as Gamma
from astropy import units as u
from astropy.coordinates import SkyCoord


#plot the filter selections 
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman') 
plt.rcParams.update({'font.size': 16})


# read ra,dec from plaeposit.dat
ora, odec = np.loadtxt("plaeposit.dat",unpack=True)
rora, rodec = np.loadtxt("randgal.dat",unpack=True)

# ra,dec box rabox = [218.9,217.55] decbox = [32.95,33.85]
ramax, ramin = [218.9,217.55]
decmin, decmax = [32.95,33.85]

#%%
""" Make "catbox.dat" ==========================================
details :
    for ./geach12_model1/ and ./geach12_model2/,  generate dirs cat0,cat1, ... , cat59 and box-cut catalog files "catbox.dat"
"""
dirheader = "./geach12_model1/"
catheader = "cat"

for imodel in range(0,60):
    # make a dir for each catalog
    curdir = dirheader+catheader+str(imodel)
    print curdir
    if not os.path.exists(curdir):
        os.makedirs(curdir)
    
    curcatFileName = dirheader+catheader+'_'+str(imodel)
    #print curcatFileName
    curra, curdec = np.loadtxt(curcatFileName,unpack=True)
    idxBox = np.where((curra < ramax) & (curra > ramin) & (curdec > decmin) & (curdec < decmax))
    #print idxBox
    #print len(curra)
    #print len(curra[idxBox])
    np.savetxt(curdir+"/catbox.dat",np.transpose([curra[idxBox],curdec[idxBox]]),fmt="%.6f %.6f")

dirheader = "./geach12_model2/"
catheader = "cat"

for imodel in range(0,60):
    # make a dir for each catalog
    curdir = dirheader+catheader+str(imodel)
    print curdir
    if not os.path.exists(curdir):
        os.makedirs(curdir)
    
    curcatFileName = dirheader+catheader+'_'+str(imodel)
    #print curcatFileName
    curra, curdec = np.loadtxt(curcatFileName,unpack=True)
    idxBox = np.where((curra < ramax) & (curra > ramin) & (curdec > decmin) & (curdec < decmax))
    #print idxBox
    #print len(curra)
    #print len(curra[idxBox])
    np.savetxt(curdir+"/catbox.dat",np.transpose([curra[idxBox],curdec[idxBox]]),fmt="%.6f %.6f")


""" Make "catbox.lsout" =====================================
Now using "catbox.dat", measure two-point correlation function
"""
dirheader = "./geach12_model1/"
catheader = "cat"

cmdheader = "~/cosmologybin/landyang.out "
cmdtail = " ~/cosmologybin/laebootesrandom1milcut.dat 1 10000 20"
cmdls =""
cmdcopy =""

shellscript = open("genLSmodel1.sh",'w')
lslistfile = open("LSmodel1.list",'w')

for imodel in range(0,60):
    # make a dir for each catalog
    curcatbox = dirheader+catheader+str(imodel)+"/catbox.dat"
    cmdls = cmdheader+curcatbox+cmdtail
    print cmdls
    cmdcopy = "cp ./_lsout.data "+dirheader+catheader+str(imodel)+"/catbox.lsout"
    print cmdcopy
    shellscript.write(cmdls+"\n")
    shellscript.write(cmdcopy+"\n")
    print dirheader+catheader+str(imodel)+"/catbox.lsout"
    lslistfile.write(dirheader+catheader+str(imodel)+"/catbox.lsout"+"\n")

shellscript.close()
lslistfile.close()



dirheader = "./geach12_model2/"
catheader = "cat"

cmdheader = "~/cosmologybin/landyang.out "
cmdtail = " ~/cosmologybin/laebootesrandom1milcut.dat 1 10000 20"
cmdls =""
cmdcopy =""

shellscript = open("genLSmodel2.sh",'w')
lslistfile = open("LSmodel2.list",'w')

for imodel in range(0,60):
    # make a dir for each catalog
    curcatbox = dirheader+catheader+str(imodel)+"/catbox.dat"
    cmdls = cmdheader+curcatbox+cmdtail
    print cmdls
    cmdcopy = "cp ./_lsout.data "+dirheader+catheader+str(imodel)+"/catbox.lsout"
    print cmdcopy
    shellscript.write(cmdls+"\n")
    shellscript.write(cmdcopy+"\n")
    print dirheader+catheader+str(imodel)+"/catbox.lsout"
    lslistfile.write(dirheader+catheader+str(imodel)+"/catbox.lsout"+"\n")

shellscript.close()
lslistfile.close()


"""
Now cat genLSmodel1.sh genLSmodel2.sh > genLSall.sh 
genLSall.sh will do every calculation.

    - then, compare these with observed LAE clustering 
"""
logang, wang = np.loadtxt("../plaelsoutcut.data",skiprows=2 ,usecols =(0,4), unpack=True)

ang = np.power(10.0,logang)
wangerr = np.loadtxt("../plaebooterror.data",skiprows=1 ,usecols = (2,), unpack=True)

print logang
print ang
print wang
print wangerr
print len(ang),len(wang),len(wangerr)


fig = plt.figure(figsize=(16,9))
plt.subplot(231)

#plt.plot(logang,wang,'s')
plt.errorbar(ang,wang,yerr=wangerr,fmt='o')
plt.axis([1,1000,0.001,50])
plt.title("Least Chi-Square Model")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$\theta$ (arcsec)')
plt.ylabel(r'$\omega_{LS}(\theta)$')
plt.scatter(ang[:14],wang[:14],marker='x',color='r',s=100)



dirheader = "./geach12_model1/"
catheader = "cat"
curlsname = ""

chivalues = np.zeros(60)
ngals = np.zeros(60)
smallchivalues = np.zeros(60)

for imodel in range(0,60):
    curlsname = dirheader+catheader+str(imodel)+"/catbox.lsout"
    curlogang, curwang = np.loadtxt(curlsname,skiprows=2 ,usecols =(0,4), unpack=True)
    curang = np.power(10.0,curlogang)
    plt.plot(curang[:14],curwang[:14],color='grey',linestyle='dotted')
    chisquare = np.sum(((wang[:14] - curwang[:14])/wangerr[:14])**2)
    chivalues[imodel] = chisquare

    smallchisquare = np.sum(((wang[:4] - curwang[:4])/wangerr[:4])**2)
    smallchivalues[imodel] = smallchisquare


    curlsname = dirheader+catheader+str(imodel)+"/catbox.dat"
    curra, curdec = np.loadtxt(curlsname, unpack=True)
    ngals[imodel] = len(curra)
    print imodel,ngals[imodel],chivalues[imodel],smallchivalues[imodel]


ibest = np.argsort(chivalues)
for idx in range(0,5):
    bestlsname = dirheader+catheader+str(ibest[idx])+"/catbox.lsout"
    bestlogang, bestwang = np.loadtxt(bestlsname,skiprows=2 ,usecols =(0,4), unpack=True)
    if idx == 0 :
        plt.plot(curang[:14],bestwang[:14],color='red',linestyle='solid',linewidth=2)
    else :
        plt.plot(curang[:14],bestwang[:14],color='red',linestyle='dotted',linewidth=1.5)
      


plt.subplot(232)
plt.yscale('log')
plt.axis([1000,1700,1,1000])
plt.title("Least Chi-Square Model")
plt.xlabel(r'$N_{galaxy}$')
plt.ylabel(r'$\chi^2$')
plt.scatter(ngals,chivalues,color='grey')
plt.scatter(ngals[ibest[:5]],chivalues[ibest[:5]],marker='x',color='r',s=100)


plt.subplot(233)
plt.yscale('log')
plt.axis([1000,1700,1,100])
plt.title("Least Chi-Square Model")
plt.xlabel(r'$N_{galaxy}$')
plt.ylabel(r'$\chi^2$ ($\theta < 10$ arcsec)')
plt.scatter(ngals,smallchivalues,color='grey')
plt.scatter(ngals[ibest[:5]],smallchivalues[ibest[:5]],marker='x',color='r',s=100)


#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#plt.show()


# model 2

plt.subplot(234)

#plt.plot(logang,wang,'s')
plt.errorbar(ang,wang,yerr=wangerr,fmt='o')
plt.axis([1,1000,0.001,50])
plt.title("MCMC Model")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$\theta$ (arcsec)')
plt.ylabel(r'$\omega_{LS}(\theta)$')
plt.scatter(ang[:14],wang[:14],marker='x',color='r',s=100)


dirheader = "./geach12_model2/"
catheader = "cat"
curlsname = ""

chivalues = np.zeros(60)
ngals = np.zeros(60)
smallchivalues = np.zeros(60)

for imodel in range(0,60):
    curlsname = dirheader+catheader+str(imodel)+"/catbox.lsout"
    curlogang, curwang = np.loadtxt(curlsname,skiprows=2 ,usecols =(0,4), unpack=True)
    curang = np.power(10.0,curlogang)
    plt.plot(curang[:14],curwang[:14],color='grey',linestyle='dotted')
    chisquare = np.sum(((wang[:14] - curwang[:14])/wangerr[:14])**2)
    chivalues[imodel] = chisquare

    smallchisquare = np.sum(((wang[:4] - curwang[:4])/wangerr[:4])**2)
    smallchivalues[imodel] = smallchisquare


    curlsname = dirheader+catheader+str(imodel)+"/catbox.dat"
    curra, curdec = np.loadtxt(curlsname, unpack=True)
    ngals[imodel] = len(curra)
    print imodel,ngals[imodel],chivalues[imodel],smallchivalues[imodel]

ibest = np.argsort(chivalues)
for idx in range(0,5):
    bestlsname = dirheader+catheader+str(ibest[idx])+"/catbox.lsout"
    bestlogang, bestwang = np.loadtxt(bestlsname,skiprows=2 ,usecols =(0,4), unpack=True)
    if idx == 0 :
        plt.plot(curang[:14],bestwang[:14],color='red',linestyle='solid',linewidth=2)
    else :
        plt.plot(curang[:14],bestwang[:14],color='red',linestyle='dotted',linewidth=1.5)
      

plt.subplot(235)
plt.yscale('log')
plt.axis([1000,1700,1,1000])
plt.title("MCMC Model")
plt.xlabel(r'$N_{galaxy}$')
plt.ylabel(r'$\chi^2$')
plt.scatter(ngals,chivalues,color='grey')
plt.scatter(ngals[ibest[:5]],chivalues[ibest[:5]],marker='x',color='r',s=100)


plt.subplot(236)
plt.yscale('log')
plt.axis([1000,1700,1,100])
plt.title("MCMC Model")
plt.xlabel(r'$N_{galaxy}$')
plt.ylabel(r'$\chi^2$ ($\theta < 10$ arcsec)')
plt.scatter(ngals,smallchivalues,color='grey')
plt.scatter(ngals[ibest[:5]],smallchivalues[ibest[:5]],marker='x',color='r',s=100)


plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.show()



""" Make a script to generate network edges "dij.dat" for each "catbox.dat"

1. basic command
    ~/networkbin/getpairarcsecdistInputRADECDegree.bin catbox.dat dij.dat

2. running script "genDij.sh"

"""

shellscript = open("genDij.sh",'w')

cmdheader = "~/networkbin/getpairarcsecdistInputRADECDegree.bin "
cmdtail = ""

catheader = "cat"
cmd =""


dirheader = "./geach12_model1/"
for imodel in range(0,60):
    # make a dir for each catalog
    curcatbox = dirheader+catheader+str(imodel)+"/catbox.dat "
    curdij = dirheader+catheader+str(imodel)+"/dij.dat"
    cmd = cmdheader+curcatbox+curdij
    print cmd
    shellscript.write(cmd+"\n")

dirheader = "./geach12_model2/"
for imodel in range(0,60):
    # make a dir for each catalog
    curcatbox = dirheader+catheader+str(imodel)+"/catbox.dat "
    curdij = dirheader+catheader+str(imodel)+"/dij.dat"
    cmd = cmdheader+curcatbox+curdij
    print cmd
    shellscript.write(cmd+"\n")

shellscript.close()




"""
Generating edge files : R??arcsec??.edge from "dij.dat"


"""
#shellscript = open("getNetwork.sh",'w')


cmdheader = "~/networkbin/netscalar.out "
cmdtail = ""

catheader = "cat"
cmd =""
fedgename=""
curdir=""
curdijfile=""

arcseclist = np.linspace(1,200,num=200,endpoint=True)
#arcseclist = np.linspace(1,5,num=5,endpoint=True)



dirheader = "./geach12_model1/"
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    fedgename = curdir+"edge.list"
    #print fedgename
    fedgelist = open(fedgename,'w')
    
    
    tmpstr = ''
    tmpfloat = 0.0

    curdijfile = curdir+"dij.dat"
    oi,oj,odij = np.loadtxt(curdijfile,unpack=True)
    numtotaledge = len(oi)

    for ridx, curarcsec in enumerate(arcseclist):
        curfname=str("R%04d.edge"%curarcsec)
        curnetname=str("R%04d.network"%curarcsec)
 
        ## add this line to the edgelist file
        fedgelist.write(curfname+'\n')
        
        #make cmd 
        cmd = cmdheader+curdir+curfname+" "+curdir+curnetname
        print curdir+curfname
        print curdir+curnetname
        print cmd        
#        shellscript.write(cmd+'\n')
        
        curfile=open(curdir+curfname,'w')

        for iedge in range(numtotaledge):
            if odij[iedge] <= curarcsec:
                #print "%d %d %f\n"%(oi[iedge], oj[iedge], odij[iedge])
                curfile.write("%d %d\n"%(oi[iedge], oj[iedge]))
        curfile.close()

    fedgelist.close()

dirheader = "./geach12_model2/"
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    fedgename = curdir+"edge.list"
    #print fedgename
    fedgelist = open(fedgename,'w')
    
    
    tmpstr = ''
    tmpfloat = 0.0

    curdijfile = curdir+"dij.dat"
    oi,oj,odij = np.loadtxt(curdijfile,unpack=True)
    numtotaledge = len(oi)

    for ridx, curarcsec in enumerate(arcseclist):
        curfname=str("R%04d.edge"%curarcsec)
        curnetname=str("R%04d.network"%curarcsec)
 
        ## add this line to the edgelist file
        fedgelist.write(curfname+'\n')
        
        #make cmd 
        cmd = cmdheader+curdir+curfname+" "+curdir+curnetname
        print curdir+curfname
        print curdir+curnetname
        print cmd        
#        shellscript.write(cmd+'\n')
        
        curfile=open(curdir+curfname,'w')

        for iedge in range(numtotaledge):
            if odij[iedge] <= curarcsec:
                #print "%d %d %f\n"%(oi[iedge], oj[iedge], odij[iedge])
                curfile.write("%d %d\n"%(oi[iedge], oj[iedge]))
        curfile.close()

    fedgelist.close()


#shellscript.close()






"""
Now parsing the network result files 


"""



catheader = "cat"
fedgename=""
curdir=""
curcatfile=""
curnetworkfile=""

arcseclist = np.linspace(1,200,num=200,endpoint=True)
#arcseclist = np.linspace(1,5,num=5,endpoint=True)


## Model 1

dirheader = "./geach12_model1/"
catheader = "cat"
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    fedgename = curdir+"edge.list"
    curcatfile= curdir+"catbox.dat"
    #print fedgename
    
    #read a networkfile list
    filelist = open(fedgename,"r")
    eachfile = filelist.readlines()
    numlist = len(eachfile)
    filelist.close()
    print 'Num of network files = ', numlist

    linkLength = np.zeros(numlist)
    diameter = np.zeros(numlist)
    transitivity = np.zeros(numlist)
    avgclustering = np.zeros(numlist)
    gcomponent = np.zeros(numlist)
    numvert = np.zeros(numlist)
    density = np.zeros(numlist)
    cendg = np.zeros(numlist)
    cenbc = np.zeros(numlist)
    cencl = np.zeros(numlist)
    maxclique = np.zeros(numlist)



    tmpfile = open(curcatfile,"r")
    tmplines = tmpfile.readlines()
    numsize = np.float_(len(tmplines))
    tmpfile.close()
    numsize = np.float_(len(tmplines))
    #print eachfile

    for idx, curname in enumerate(eachfile, start=0):
        #parse the filename to load it and to extract Radius
        curname = curname.rstrip()
        curname = curname.replace("edge","network")
        curRadius = np.float_(curname.split("R")[1].split(".n")[0])
        curname = curdir+curname
        print 'current file = ', curname
        print "curent radius = ", curRadius
    
        #parse the content of each result file
        curfile = open(curname,"r")
        curresult = curfile.readlines()
        print "current result = ", curresult
        curDia = np.float_(curresult[0].rstrip().split(':')[1])
        print "current diameter = ", curDia
        curClu = np.float_(curresult[1].rstrip().split(':')[1])
        print "current transitivity = ", curClu
        curAvgClu = np.float_(curresult[2].rstrip().split(':')[1])
        print "current AvgClustering = ", curAvgClu
        curGcomp = np.float_(curresult[3].rstrip().split(':')[1])
        print "current Giant Component = ", curGcomp
        curNumvert = np.float_(curresult[4].rstrip().split(':')[1])
        print "current Number of Vertices = ", curNumvert
        curDensity = np.float_(curresult[5].rstrip().split(':')[1])
        print "current Edge Density = ", curDensity
        curCenDG = np.float_(curresult[6].rstrip().split(':')[1])
        print "current CenDG = ", curCenDG
        curCenBC = np.float_(curresult[7].rstrip().split(':')[1])
        print "current CenDG = ", curCenBC
        curCenCL = np.float_(curresult[8].rstrip().split(':')[1])
        print "current CenDG = ", curCenCL
        curMaxClique = np.float_(curresult[9].rstrip().split(':')[1])
        print "current macClique = ", curMaxClique
        curfile.close()
    
        linkLength[idx] = curRadius
        diameter[idx] = curDia
        transitivity[idx] = curClu
        avgclustering[idx] = curAvgClu
        gcomponent[idx] = curGcomp/numsize
        numvert[idx] = curNumvert/numsize
        density[idx] = curDensity
        cendg[idx] = curCenDG
        cenbc[idx] = curCenBC
        cencl[idx] = curCenCL
        maxclique[idx] = curMaxClique

    print linkLength
    isort = np.argsort(linkLength)
    print linkLength[isort]

    linkLength = linkLength[isort]
    diameter = diameter[isort]
    transitivity = transitivity[isort]
    avgclustering = avgclustering[isort]
    gcomponent = gcomponent[isort]
    numvert = numvert[isort]
    density = density[isort]
    cendg = cendg[isort]
    cenbc = np.absolute(cenbc[isort])
    cencl = cencl[isort]
    maxclique = maxclique[isort]


    # write out the results
    curnetworkfile=curdir+"network.results"
    ofile = open(curnetworkfile,"w")
    data = np.array([linkLength,diameter,transitivity,avgclustering,gcomponent,numvert,density,cendg,cenbc,cencl,maxclique])
    data = data.T
    np.savetxt(ofile, data, fmt=['%f','%f','%f','%f','%f','%f','%f','%f','%f','%f','%f'])
    ofile.close()



### Model 2

dirheader = "./geach12_model2/"
catheader = "cat"
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    fedgename = curdir+"edge.list"
    curcatfile= curdir+"catbox.dat"
    #print fedgename
    
    #read a networkfile list
    filelist = open(fedgename,"r")
    eachfile = filelist.readlines()
    numlist = len(eachfile)
    filelist.close()
    print 'Num of network files = ', numlist

    linkLength = np.zeros(numlist)
    diameter = np.zeros(numlist)
    transitivity = np.zeros(numlist)
    avgclustering = np.zeros(numlist)
    gcomponent = np.zeros(numlist)
    numvert = np.zeros(numlist)
    density = np.zeros(numlist)
    cendg = np.zeros(numlist)
    cenbc = np.zeros(numlist)
    cencl = np.zeros(numlist)
    maxclique = np.zeros(numlist)



    tmpfile = open(curcatfile,"r")
    tmplines = tmpfile.readlines()
    numsize = np.float_(len(tmplines))
    tmpfile.close()
    numsize = np.float_(len(tmplines))
    #print eachfile

    for idx, curname in enumerate(eachfile, start=0):
        #parse the filename to load it and to extract Radius
        curname = curname.rstrip()
        curname = curname.replace("edge","network")
        curRadius = np.float_(curname.split("R")[1].split(".n")[0])
        curname = curdir+curname
        print 'current file = ', curname
        print "curent radius = ", curRadius
    
        #parse the content of each result file
        curfile = open(curname,"r")
        curresult = curfile.readlines()
        print "current result = ", curresult
        curDia = np.float_(curresult[0].rstrip().split(':')[1])
        print "current diameter = ", curDia
        curClu = np.float_(curresult[1].rstrip().split(':')[1])
        print "current transitivity = ", curClu
        curAvgClu = np.float_(curresult[2].rstrip().split(':')[1])
        print "current AvgClustering = ", curAvgClu
        curGcomp = np.float_(curresult[3].rstrip().split(':')[1])
        print "current Giant Component = ", curGcomp
        curNumvert = np.float_(curresult[4].rstrip().split(':')[1])
        print "current Number of Vertices = ", curNumvert
        curDensity = np.float_(curresult[5].rstrip().split(':')[1])
        print "current Edge Density = ", curDensity
        curCenDG = np.float_(curresult[6].rstrip().split(':')[1])
        print "current CenDG = ", curCenDG
        curCenBC = np.float_(curresult[7].rstrip().split(':')[1])
        print "current CenDG = ", curCenBC
        curCenCL = np.float_(curresult[8].rstrip().split(':')[1])
        print "current CenDG = ", curCenCL
        curMaxClique = np.float_(curresult[9].rstrip().split(':')[1])
        print "current macClique = ", curMaxClique
        curfile.close()
    
        linkLength[idx] = curRadius
        diameter[idx] = curDia
        transitivity[idx] = curClu
        avgclustering[idx] = curAvgClu
        gcomponent[idx] = curGcomp/numsize
        numvert[idx] = curNumvert/numsize
        density[idx] = curDensity
        cendg[idx] = curCenDG
        cenbc[idx] = curCenBC
        cencl[idx] = curCenCL
        maxclique[idx] = curMaxClique

    print linkLength
    isort = np.argsort(linkLength)
    print linkLength[isort]

    linkLength = linkLength[isort]
    diameter = diameter[isort]
    transitivity = transitivity[isort]
    avgclustering = avgclustering[isort]
    gcomponent = gcomponent[isort]
    numvert = numvert[isort]
    density = density[isort]
    cendg = cendg[isort]
    cenbc = np.absolute(cenbc[isort])
    cencl = cencl[isort]
    maxclique = maxclique[isort]


    # write out the results
    curnetworkfile=curdir+"network.results"
    ofile = open(curnetworkfile,"w")
    data = np.array([linkLength,diameter,transitivity,avgclustering,gcomponent,numvert,density,cendg,cenbc,cencl,maxclique])
    data = data.T
    np.savetxt(ofile, data, fmt=['%f','%f','%f','%f','%f','%f','%f','%f','%f','%f','%f'])
    ofile.close()




#%%
"""
Plot network results 
"""


# "Network.result" ./randomdist/Nework.reslt
rang, rdia, rtr, ravgcc, rgcomp  = np.loadtxt("../network/randomdist/oneinstance/Network.results",unpack=True)

oang, odia, otr, oavgcc, ogcomp  = np.loadtxt("../network/Network.results",unpack=True)



fig = plt.figure(figsize=(16,13))



plt.subplot(331)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Diameter')
plt.axis([1,200,0,150])
#plt.yscale('log')
#plt.xscale('log')

dirheader = "./geach12_model1/"
catheader = "cat"
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,dia,color='b')
plt.plot(rang,rdia,color='k',linewidth=2.0)
plt.plot(oang,odia,color='r',linewidth=2.0)




plt.subplot(332)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Giant Component Fraction')
plt.axis([1,200,0.0,1.1])
#plt.xscale('log')
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,gcomp,color='b')
plt.plot(rang,rgcomp,color='k',linewidth=2.0)
plt.plot(oang,ogcomp,color='r',linewidth=2.0)



plt.subplot(333)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Transitivity')
plt.axis([1,200,0.5,0.85])
#plt.xscale('log')
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,tr,color='b')
plt.plot(rang,rtr,color='k',linewidth=2.0)
plt.plot(oang,otr,color='r',linewidth=2.0)


plt.subplot(334)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Average CC')
plt.axis([1,200,0.5,0.85])
#plt.xscale('log')
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,avgcc,color='b',linewidth=0.5)
plt.plot(rang,ravgcc,color='black',linewidth=2.0)
plt.plot(oang,oavgcc,color='r',linewidth=2.0)



plt.subplot(335)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Edge Density')
plt.axis([1,200,0.0,0.01])
#plt.xscale('log')
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,den,color='b',linewidth=0.5)


plt.subplot(336)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Degree')
plt.axis([1,200,0.0,0.022])
#plt.xscale('log')
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,cendg,color='b',linewidth=0.5)


plt.subplot(337)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'BC')
plt.axis([1,200,0.0,0.5])
#plt.xscale('log')
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,cenbc,color='b',linewidth=0.5)


plt.subplot(338)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'CL')
plt.axis([1,200,0.00001,0.08])
plt.yscale('log')
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,cencl,color='b',linewidth=0.5)

plt.subplot(339)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Max Clique')
plt.axis([1,200,0,25])
#plt.yscale('log')
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,maxclique,color='b',linewidth=0.5)




plt.show()







fig = plt.figure(figsize=(16,13))
ncat = 60

plt.subplot(331)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Diameter')
plt.axis([1,200,0,150])
#plt.yscale('log')
#plt.xscale('log')

dirheader = "./geach12_model2/"
catheader = "cat"
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,dia,color='b')
plt.plot(rang,rdia,color='k',linewidth=2.0)
plt.plot(oang,odia,color='r',linewidth=2.0)




plt.subplot(332)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Giant Component Fraction')
plt.axis([1,200,0.0,1.1])
#plt.xscale('log')
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,gcomp,color='b')
plt.plot(rang,rgcomp,color='k',linewidth=2.0)
plt.plot(oang,ogcomp,color='r',linewidth=2.0)



plt.subplot(333)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Transitivity')
plt.axis([1,200,0.5,0.85])
#plt.xscale('log')
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,tr,color='b')
plt.plot(rang,rtr,color='k',linewidth=2.0)
plt.plot(oang,otr,color='r',linewidth=2.0)


plt.subplot(334)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Average CC')
plt.axis([1,200,0.5,0.85])
#plt.xscale('log')
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,avgcc,color='b',linewidth=0.5)
plt.plot(rang,ravgcc,color='black',linewidth=2.0)
plt.plot(oang,oavgcc,color='r',linewidth=2.0)



plt.subplot(335)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Edge Density')
plt.axis([1,200,0.0,0.01])
#plt.xscale('log')
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,den,color='b',linewidth=0.5)


plt.subplot(336)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Degree')
plt.axis([1,200,0.0,0.022])
#plt.xscale('log')
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,cendg,color='b',linewidth=0.5)


plt.subplot(337)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'BC')
plt.axis([1,200,0.0,0.5])
#plt.xscale('log')
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,cenbc,color='b',linewidth=0.5)


plt.subplot(338)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'CL')
plt.axis([1,200,0.00001,0.08])
plt.yscale('log')
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,cencl,color='b',linewidth=0.5)

plt.subplot(339)
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Max Clique')
plt.axis([1,200,0,25])
#plt.yscale('log')
for imodel in range(0,ncat):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    plt.plot(ang,maxclique,color='b',linewidth=0.5)


plt.show()


#%%




##########################
# Now geting 90% confidence lines (5% cuts at the bottom and top layers)
#
# 60x200 data set for each network measurement 
rang, rdia, rtr, ravgcc, rgcomp  = np.loadtxt("../network/randomdist/oneinstance/Network.results",unpack=True)
oang, odia, otr, oavgcc, ogcomp  = np.loadtxt("../network/Network.results",unpack=True)



dirheader = "./geach12_model1/"
catheader = "cat"

diabox1 = np.zeros((60,200)) 
trbox1 = np.zeros((60,200))  
avgccbox1 = np.zeros((60,200))  
gcompbox1 = np.zeros((60,200))
nvertbox1 = np.zeros((60,200))  
denbox1 = np.zeros((60,200))  
cendgbox1 = np.zeros((60,200))  
cenbcbox1 = np.zeros((60,200))  
cenclbox1 = np.zeros((60,200))  
maxcliquebox1 = np.zeros((60,200))  



for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    for j in range(200):
        diabox1[imodel,j]=dia[j]
        trbox1[imodel,j]=tr[j]
        avgccbox1[imodel,j]=avgcc[j]
        gcompbox1[imodel,j]=gcomp[j]
        nvertbox1[imodel,j]=nvert[j]
        denbox1[imodel,j]=den[j]                
        cendgbox1[imodel,j]=cendg[j]
        cenbcbox1[imodel,j]=cenbc[j]
        cenclbox1[imodel,j]=cencl[j]
        maxcliquebox1[imodel,j]=maxclique[j]


# sort the results for getting 5% 50% 95% positions
sdiabox1 = np.zeros((60,200)) 
strbox1 = np.zeros((60,200))  
savgccbox1 = np.zeros((60,200))  
sgcompbox1 = np.zeros((60,200))  
snvertbox1 = np.zeros((60,200))  
sdenbox1 = np.zeros((60,200))  
scendgbox1 = np.zeros((60,200))  
scenbcbox1 = np.zeros((60,200))  
scenclbox1 = np.zeros((60,200))  
smaxcliquebox1 = np.zeros((60,200))  

for ir in range(200):
    sdiabox1[:,ir] = np.sort(diabox1[:,ir])
    strbox1[:,ir] = np.sort(trbox1[:,ir])
    savgccbox1[:,ir] = np.sort(avgccbox1[:,ir])
    sgcompbox1[:,ir] = np.sort(gcompbox1[:,ir])
    snvertbox1[:,ir] = np.sort(nvertbox1[:,ir])
    sdenbox1[:,ir] = np.sort(denbox1[:,ir])
    scendgbox1[:,ir] = np.sort(cendgbox1[:,ir])
    scenbcbox1[:,ir] = np.sort(cenbcbox1[:,ir])
    scenclbox1[:,ir] = np.sort(cenclbox1[:,ir])
    smaxcliquebox1[:,ir] = np.sort(maxcliquebox1[:,ir])
   

 
dirheader = "./geach12_model2/"
catheader = "cat"

diabox2 = np.zeros((60,200)) 
trbox2 = np.zeros((60,200))  
avgccbox2 = np.zeros((60,200))  
gcompbox2 = np.zeros((60,200))
nvertbox2 = np.zeros((60,200))  
denbox2 = np.zeros((60,200))  
cendgbox2 = np.zeros((60,200))  
cenbcbox2 = np.zeros((60,200))  
cenclbox2 = np.zeros((60,200))  
maxcliquebox2 = np.zeros((60,200))  



for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    curnetworkfile = curdir+"network.results"
    
    ang, dia, tr, avgcc, gcomp, nvert, den, cendg, cenbc, cencl, maxclique = np.loadtxt(curnetworkfile,unpack=True)
    for j in range(200):
        diabox2[imodel,j]=dia[j]
        trbox2[imodel,j]=tr[j]
        avgccbox2[imodel,j]=avgcc[j]
        gcompbox2[imodel,j]=gcomp[j]
        nvertbox2[imodel,j]=nvert[j]
        denbox2[imodel,j]=den[j]                
        cendgbox2[imodel,j]=cendg[j]
        cenbcbox2[imodel,j]=cenbc[j]
        cenclbox2[imodel,j]=cencl[j]
        maxcliquebox2[imodel,j]=maxclique[j]


# sort the results for getting 5% 50% 95% positions
sdiabox2 = np.zeros((60,200)) 
strbox2 = np.zeros((60,200))  
savgccbox2 = np.zeros((60,200))  
sgcompbox2 = np.zeros((60,200))  
snvertbox2 = np.zeros((60,200))  
sdenbox2 = np.zeros((60,200))  
scendgbox2 = np.zeros((60,200))  
scenbcbox2 = np.zeros((60,200))  
scenclbox2 = np.zeros((60,200))  
smaxcliquebox2 = np.zeros((60,200))  

for ir in range(200):
    sdiabox2[:,ir] = np.sort(diabox2[:,ir])
    strbox2[:,ir] = np.sort(trbox2[:,ir])
    savgccbox2[:,ir] = np.sort(avgccbox2[:,ir])
    sgcompbox2[:,ir] = np.sort(gcompbox2[:,ir])
    snvertbox2[:,ir] = np.sort(nvertbox2[:,ir])
    sdenbox2[:,ir] = np.sort(denbox2[:,ir])
    scendgbox2[:,ir] = np.sort(cendgbox2[:,ir])
    scenbcbox2[:,ir] = np.sort(cenbcbox2[:,ir])
    scenclbox2[:,ir] = np.sort(cenclbox2[:,ir])
    smaxcliquebox2[:,ir] = np.sort(maxcliquebox2[:,ir])
   

import pickle
with open('../obsdist/obsdistNetwork.pickle') as f:
    sdiaboxobs,strboxobs,savgccboxobs,sgcompboxobs,snvertboxobs,sdenboxobs,scendgboxobs,scenbcboxobs,scenclboxobs,smaxcliqueboxobs= pickle.load(f)

with open('../randomdist/randNetwork.pickle') as f:
    sdiaboxrand,strboxrand,savgccboxrand,sgcompboxrand,snvertboxrand,sdenboxrand,scendgboxrand,scenbcboxrand,scenclboxrand,smaxcliqueboxrand = pickle.load(f)



#plot the filter selections 
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman') 
plt.rcParams.update({'font.size': 14,'legend.fontsize': 14})
  
   
import matplotlib.patches as mp
obsp = mp.Patch(color='red', label='Observed LAEs')
m1p = mp.Patch(color='green', label='Model #1')
m2p = mp.Patch(color='blue', label='Model #2')
rp = mp.Patch(color='grey', label='Random')


fig = plt.figure(figsize=(10,12))


plt.subplot(421)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Diameter')
plt.axis([1,200,0,120])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,sdiabox1[30,:],color='green')
plt.fill_between(ang, sdiabox1[3,:],sdiabox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, sdiaboxobs[3,:],sdiaboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,sdiaboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, sdiaboxrand[3,:],sdiaboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(422)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Giant Component Fraction')
plt.axis([1,200,0,1.1])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,sgcompbox1[30,:],color='green')
plt.fill_between(ang, sgcompbox1[3,:],sgcompbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, sgcompboxobs[3,:],sgcompboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,sgcompboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, sgcompboxrand[3,:],sgcompboxrand[56,:], facecolor='grey', alpha=0.2)



plt.subplot(423)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Average CC')
plt.axis([1,200,0.5,0.8])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=1)

#plt.plot(ang,savgccbox1[30,:],color='green')
plt.fill_between(ang, savgccbox1[3,:],savgccbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, savgccboxobs[3,:],savgccboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,savgccboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, savgccboxrand[3,:],savgccboxrand[56,:], facecolor='grey', alpha=0.2)
plt.axvline(170, color='grey', linestyle='--')
plt.text(170, 0.52, r'$ >170^{\prime\prime}$',color='grey', fontsize=13)
plt.text(105, 0.75, '\"CC170\"',color='green', fontsize=25)


plt.subplot(424)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Transitivity')
plt.axis([1,200,0.5,0.8])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=1)
plt.axvline(70, color='grey', linestyle='--')
plt.text(75, 0.52, r'$70^{\prime\prime} = 1.42 h^{-1} $Mpc',color='grey', fontsize=20)
plt.text(75, 0.75, '\"TR70\"',color='green', fontsize=25)


#plt.plot(ang,strbox1[30,:],color='green')
plt.fill_between(ang, strbox1[3,:],strbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, strboxobs[3,:],strboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,strboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, strboxrand[3,:],strboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(425)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Edge Density')
plt.axis([1,200,0.0,0.01])
#plt.yscale('log')
#plt.xscale('log')
plt.axvline(100, color='grey', linestyle='--')
plt.text(105, 0.0015, '\"ED100\"',color='green', fontsize=25)
plt.text(105, 0.0003, r'$ >100^{\prime\prime}$',color='grey', fontsize=13)
plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,sdenbox1[30,:],color='green')
plt.fill_between(ang, sdenbox1[3,:],sdenbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, sdenboxobs[3,:],sdenboxobs[56,:], facecolor='red', alpha=0.3)
#plt.plot(ang,sdenboxobs[30,:],color='red',linewidth=0.5)
plt.fill_between(ang, sdenboxrand[3,:],sdenboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(426)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Size of the Largest Clique')
plt.axis([1,200,0,25])
#plt.yscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,smaxcliquebox1[30,:],color='green')
plt.fill_between(ang, smaxcliquebox1[3,:],smaxcliquebox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, smaxcliqueboxobs[3,:],smaxcliqueboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,smaxcliqueboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, smaxcliqueboxrand[3,:],smaxcliqueboxrand[56,:], facecolor='grey', alpha=0.2)
plt.text(75, 14, r'$70^{\prime\prime}$',color='grey', fontsize=20)
#plt.text(75, 2.3, r'Excess near $70^{\prime\prime}$',color='grey', fontsize=20)
plt.arrow(70, 16, 0, -2, head_width=5.0, head_length=1.0, fc='grey', ec='grey')


plt.subplot(427)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Betweenness Centralization')
plt.axis([1,200,0.0,0.5])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,scenbcbox1[30,:],color='green')
plt.fill_between(ang, scenbcbox1[3,:],scenbcbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, scenbcboxobs[3,:],scenbcboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scenbcboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scenbcboxrand[3,:],scenbcboxrand[56,:], facecolor='grey', alpha=0.2)



plt.subplot(428)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Degree Centralization')
plt.axis([1,200,0.0,0.025])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,scendgbox1[30,:],color='green')
plt.fill_between(ang, scendgbox1[3,:],scendgbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, scendgboxobs[3,:],scendgboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scendgboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scendgboxrand[3,:],scendgboxrand[56,:], facecolor='grey', alpha=0.2)
plt.text(45, 0.011, r'$40^{\prime\prime}$',color='grey', fontsize=20)
#plt.text(75, 2.3, r'Excess near $70^{\prime\prime}$',color='grey', fontsize=20)
plt.arrow(40, 0.013, 0, -0.0025, head_width=5.0, head_length=0.001, fc='grey', ec='grey')



plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1network.pdf', dpi=600)
plt.show()



## model2 2x4


fig = plt.figure(figsize=(10,12))


plt.subplot(421)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Diameter')
plt.axis([1,200,0,120])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,sdiabox2[30,:],color='blue')
plt.fill_between(ang, sdiabox2[3,:],sdiabox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, sdiaboxobs[3,:],sdiaboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,sdiaboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, sdiaboxrand[3,:],sdiaboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(422)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Giant Component Fraction')
plt.axis([1,200,0,1.1])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,sgcompbox2[30,:],color='blue')
plt.fill_between(ang, sgcompbox2[3,:],sgcompbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, sgcompboxobs[3,:],sgcompboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,sgcompboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, sgcompboxrand[3,:],sgcompboxrand[56,:], facecolor='grey', alpha=0.2)



plt.subplot(423)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Average CC')
plt.axis([1,200,0.5,0.8])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=1)

#plt.plot(ang,savgccbox2[30,:],color='blue')
plt.fill_between(ang, savgccbox2[3,:],savgccbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, savgccboxobs[3,:],savgccboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,savgccboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, savgccboxrand[3,:],savgccboxrand[56,:], facecolor='grey', alpha=0.2)
plt.axvline(170, color='grey', linestyle='--')
plt.text(170, 0.52, r'$ >170^{\prime\prime}$',color='grey', fontsize=13)
plt.text(105, 0.75, '\"CC170\"',color='blue', fontsize=25)
#plt.arrow(55, 0.765, 0, -0.02, head_width=10.0, head_length=0.02, fc='blue', ec='blue')
plt.plot([55], [0.73], marker="^",color='blue',markersize=15,fillstyle="full")

plt.subplot(424)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Transitivity')
plt.axis([1,200,0.5,0.8])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=1)
plt.axvline(70, color='grey', linestyle='--')
plt.text(75, 0.52, r'$70^{\prime\prime} = 1.42 h^{-1} $Mpc',color='grey', fontsize=20)
plt.text(75, 0.75, '\"TR70\"',color='blue', fontsize=25)
#plt.arrow(40, 0.795, 0, -0.02, head_width=10.0, head_length=0.02, fc='blue', ec='blue')
plt.plot([42], [0.77], marker="^",color='blue',markersize=15,fillstyle="full")

#plt.plot(ang,strbox2[30,:],color='blue')
plt.fill_between(ang, strbox2[3,:],strbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, strboxobs[3,:],strboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,strboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, strboxrand[3,:],strboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(425)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Edge Density')
plt.axis([1,200,0.0,0.01])
#plt.yscale('log')
#plt.xscale('log')
plt.axvline(100, color='grey', linestyle='--')
plt.text(105, 0.0015, '\"ED100\"',color='blue', fontsize=25)
plt.text(105, 0.0003, r'$ >100^{\prime\prime}$',color='grey', fontsize=13)
plt.legend(handles=[m2p,obsp,rp],loc=2)

#plt.plot(ang,sdenbox2[30,:],color='blue')
plt.fill_between(ang, sdenbox2[3,:],sdenbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, sdenboxobs[3,:],sdenboxobs[56,:], facecolor='red', alpha=0.3)
#plt.plot(ang,sdenboxobs[30,:],color='red',linewidth=0.5)
plt.fill_between(ang, sdenboxrand[3,:],sdenboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(426)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Size of the Largest Clique')
plt.axis([1,200,0,25])
#plt.yscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,smaxcliquebox2[30,:],color='blue')
plt.fill_between(ang, smaxcliquebox2[3,:],smaxcliquebox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, smaxcliqueboxobs[3,:],smaxcliqueboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,smaxcliqueboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, smaxcliqueboxrand[3,:],smaxcliqueboxrand[56,:], facecolor='grey', alpha=0.2)
plt.text(75, 14, r'$70^{\prime\prime}$',color='grey', fontsize=20)
#plt.text(75, 2.3, r'Excess near $70^{\prime\prime}$',color='grey', fontsize=20)
plt.arrow(70, 16, 0, -2, head_width=5.0, head_length=1.0, fc='grey', ec='grey')
#plt.arrow(35, 13, 0, -2, head_width=10.0, head_length=2, fc='blue', ec='blue')
plt.plot([35], [10], marker="^",color='blue',markersize=15,fillstyle="full")

plt.subplot(427)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Betweenness Centralization')
plt.axis([1,200,0.0,0.5])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,scenbcbox2[30,:],color='blue')
plt.fill_between(ang, scenbcbox2[3,:],scenbcbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, scenbcboxobs[3,:],scenbcboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scenbcboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scenbcboxrand[3,:],scenbcboxrand[56,:], facecolor='grey', alpha=0.2)



plt.subplot(428)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Degree Centralization')
plt.axis([1,200,0.0,0.025])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,scendgbox2[30,:],color='blue')
plt.fill_between(ang, scendgbox2[3,:],scendgbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, scendgboxobs[3,:],scendgboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scendgboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scendgboxrand[3,:],scendgboxrand[56,:], facecolor='grey', alpha=0.2)
plt.text(45, 0.011, r'$40^{\prime\prime}$',color='grey', fontsize=20)
#plt.text(75, 2.3, r'Excess near $70^{\prime\prime}$',color='grey', fontsize=20)
plt.arrow(40, 0.013, 0, -0.0025, head_width=5.0, head_length=0.001, fc='grey', ec='grey')
#plt.arrow(20, 0.01, 0, -0.002, head_width=10.0, head_length=0.002, fc='blue', ec='blue')
plt.plot([20], [0.007], marker="^",color='blue',markersize=15,fillstyle="full")



plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model2network.pdf', dpi=600)
plt.show()


#%%

########
# getting KS test p-values
from scipy import stats


print avgccbox1[:,70]
print avgccbox2[:,70]

ksdiarand = np.zeros(200) 
ksdiamone = np.zeros(200) 
ksdiamtwo = np.zeros(200) 

kstrrand = np.zeros(200) 
kstrmone = np.zeros(200) 
kstrmtwo = np.zeros(200) 

ksavgccrand = np.zeros(200) 
ksavgccmone = np.zeros(200) 
ksavgccmtwo = np.zeros(200) 

ksgcomprand = np.zeros(200) 
ksgcompmone = np.zeros(200) 
ksgcompmtwo = np.zeros(200) 


ksdenrand = np.zeros(200) 
ksdenmone = np.zeros(200) 
ksdenmtwo = np.zeros(200) 

kscendgrand = np.zeros(200) 
kscendgmone = np.zeros(200) 
kscendgmtwo = np.zeros(200) 

kscenbcrand = np.zeros(200) 
kscenbcmone = np.zeros(200) 
kscenbcmtwo = np.zeros(200) 

kscqrand = np.zeros(200) 
kscqmone = np.zeros(200) 
kscqmtwo = np.zeros(200) 


for i in range(200):
    dum = stats.ks_2samp(sdiaboxobs[:,i],sdiabox1[:,i])
    ksdiamone[i] = np.double(dum[1])
    dum = stats.ks_2samp(sdiaboxobs[:,i],sdiabox2[:,i])
    ksdiamtwo[i] = np.double(dum[1])
    dum = stats.ks_2samp(sdiaboxobs[:,i],sdiaboxrand[:,i])
    ksdiarand[i] = np.double(dum[1])

    dum = stats.ks_2samp(sgcompboxobs[:,i],sgcompbox1[:,i])
    ksgcompmone[i] = np.double(dum[1])
    dum = stats.ks_2samp(sgcompboxobs[:,i],sgcompbox2[:,i])
    ksgcompmtwo[i] = np.double(dum[1])
    dum = stats.ks_2samp(sgcompboxobs[:,i],sgcompboxrand[:,i])
    ksgcomprand[i] = np.double(dum[1])

    iva = np.where(np.isfinite(sdenboxobs[:,i]))
    ivb = np.where(np.isfinite(sdenbox1[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(sdenboxobs[iva[0],i],sdenbox1[ivb[0],i])
        ksdenmone[i] = np.double(dum[1])
    else:
        ksdenmone[i] = np.nan
    ivb = np.where(np.isfinite(sdenbox2[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(sdenboxobs[iva[0],i],sdenbox2[ivb[0],i])
        ksdenmtwo[i] = np.double(dum[1])
    else:
        ksdenmone[i] = np.nan
    ivb = np.where(np.isfinite(sdenboxrand[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(sdenboxobs[iva[0],i],sdenboxrand[ivb[0],i])
        ksdenrand[i] = np.double(dum[1])
    else:
        ksdenrand[i] = np.nan     

    iva = np.where(np.isfinite(strboxobs[:,i]))
    ivb = np.where(np.isfinite(strbox1[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(strboxobs[iva[0],i],strbox1[ivb[0],i])
        kstrmone[i] = np.double(dum[1])
    else:
        kstrmone[i] = np.nan
    ivb = np.where(np.isfinite(strbox2[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(strboxobs[iva[0],i],strbox2[ivb[0],i])
        kstrmtwo[i] = np.double(dum[1])
    else:
        kstrmone[i] = np.nan
    ivb = np.where(np.isfinite(strboxrand[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(strboxobs[iva[0],i],strboxrand[ivb[0],i])
        kstrrand[i] = np.double(dum[1])
    else:
        kstrrand[i] = np.nan     


    iva = np.where(np.isfinite(savgccboxobs[:,i]))
    ivb = np.where(np.isfinite(savgccbox1[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(savgccboxobs[iva[0],i],savgccbox1[ivb[0],i])
        ksavgccmone[i] = np.double(dum[1])
    else:
        ksavgccmone[i] = np.nan
    ivb = np.where(np.isfinite(savgccbox2[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(savgccboxobs[iva[0],i],savgccbox2[ivb[0],i])
        ksavgccmtwo[i] = np.double(dum[1])
    else:
        ksavgccmone[i] = np.nan
    ivb = np.where(np.isfinite(savgccboxrand[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(savgccboxobs[iva[0],i],savgccboxrand[ivb[0],i])
        ksavgccrand[i] = np.double(dum[1])
    else:
        ksavgccrand[i] = np.nan     



    iva = np.where(np.isfinite(scenbcboxobs[:,i]))
    ivb = np.where(np.isfinite(scenbcbox1[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(scenbcboxobs[iva[0],i],scenbcbox1[ivb[0],i])
        kscenbcmone[i] = np.double(dum[1])
    else:
        kscenbcmone[i] = np.nan
    ivb = np.where(np.isfinite(scenbcbox2[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(scenbcboxobs[iva[0],i],scenbcbox2[ivb[0],i])
        kscenbcmtwo[i] = np.double(dum[1])
    else:
        kscenbcmone[i] = np.nan
    ivb = np.where(np.isfinite(scenbcboxrand[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(scenbcboxobs[iva[0],i],scenbcboxrand[ivb[0],i])
        kscenbcrand[i] = np.double(dum[1])
    else:
        kscenbcrand[i] = np.nan     


    iva = np.where(np.isfinite(scendgboxobs[:,i]))
    ivb = np.where(np.isfinite(scendgbox1[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(scendgboxobs[iva[0],i],scendgbox1[ivb[0],i])
        kscendgmone[i] = np.double(dum[1])
    else:
        kscendgmone[i] = np.nan
    ivb = np.where(np.isfinite(scendgbox2[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(scendgboxobs[iva[0],i],scendgbox2[ivb[0],i])
        kscendgmtwo[i] = np.double(dum[1])
    else:
        kscendgmone[i] = np.nan
    ivb = np.where(np.isfinite(scendgboxrand[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(scendgboxobs[iva[0],i],scendgboxrand[ivb[0],i])
        kscendgrand[i] = np.double(dum[1])
    else:
        kscendgrand[i] = np.nan     


    iva = np.where(np.isfinite(smaxcliqueboxobs[:,i]))
    ivb = np.where(np.isfinite(smaxcliquebox1[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(smaxcliqueboxobs[iva[0],i],smaxcliquebox1[ivb[0],i])
        kscqmone[i] = np.double(dum[1])
    else:
        kscqmone[i] = np.nan
    ivb = np.where(np.isfinite(smaxcliquebox2[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(smaxcliqueboxobs[iva[0],i],smaxcliquebox2[ivb[0],i])
        kscqmtwo[i] = np.double(dum[1])
    else:
        kscqmone[i] = np.nan
    ivb = np.where(np.isfinite(smaxcliqueboxrand[:,i]))
    checkflag = len(iva[0])*len(ivb[0])
    if checkflag > 0 :
        dum = stats.ks_2samp(smaxcliqueboxobs[iva[0],i],smaxcliqueboxrand[ivb[0],i])
        kscqrand[i] = np.double(dum[1])
    else:
        kscqrand[i] = np.nan     


#%%
#plot the filter selections 
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman') 
plt.rcParams.update({'font.size': 14,'legend.fontsize': 14})
  
   
import matplotlib.patches as mp
obsp = mp.Patch(color='red', label='Observed LAEs')
m1p = mp.Patch(color='green', label='Model #1')
m2p = mp.Patch(color='blue', label='Model #2')
rp = mp.Patch(color='grey', label='Random')

thres = np.zeros(len(ang))
thres = 1e-2 # 0.01 p value 

fig = plt.figure(figsize=(10,8))

plt.subplot(421)
plt.title('Diameter')
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.yscale('log')
plt.ylabel(r'$p$ : K-S Test')
plt.axis([1,200,1e-3,1])
plt.plot(ang,ksdiamone,'g')
plt.plot(ang,ksdiamtwo,'b')
plt.plot(ang,ksdiarand,color='grey')
plt.fill_between(ang,thres,ksdiamone,where=ksdiamone>thres ,color='g',alpha=0.2)
plt.fill_between(ang,thres,ksdiamtwo,where=ksdiamtwo>thres ,color='b',alpha=0.2)
plt.fill_between(ang,thres,ksdiarand,where=ksdiarand>thres ,color='grey',alpha=0.2)



plt.subplot(422)
plt.title('Giant Component Fraction')
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.yscale('log')
plt.ylabel(r'$p$ : K-S Test')
plt.axis([1,200,1e-3,1])
plt.plot(ang,ksgcompmone,'g')
plt.plot(ang,ksgcompmtwo,'b')
plt.plot(ang,ksgcomprand,color='grey')
plt.fill_between(ang,thres,ksgcompmone,where=ksgcompmone>thres ,color='g',alpha=0.2)
plt.fill_between(ang,thres,ksgcompmtwo,where=ksgcompmtwo>thres ,color='b',alpha=0.2)
plt.fill_between(ang,thres,ksgcomprand,where=ksgcomprand>thres ,color='grey',alpha=0.2)


plt.subplot(423)
plt.title('Average CC')
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.yscale('log')
plt.ylabel(r'$p$ : K-S Test')
plt.axis([1,200,1e-3,1])
plt.plot(ang,ksavgccmone,'g')
plt.plot(ang,ksavgccmtwo,'b')
plt.plot(ang,ksavgccrand,color='grey')
plt.fill_between(ang,thres,ksavgccmone,where=ksavgccmone>thres ,color='g',alpha=0.2)
plt.fill_between(ang,thres,ksavgccmtwo,where=ksavgccmtwo>thres ,color='b',alpha=0.2)
plt.fill_between(ang,thres,ksavgccrand,where=ksavgccrand>thres ,color='grey',alpha=0.2)
plt.axvline(170, color='grey', linestyle='--')
plt.text(170, 10**(-2.5), r'$ >170^{\prime\prime}$',color='grey', fontsize=13)
plt.text(115, 10**(-1.0), 'CC170',color='grey', fontsize=25)

plt.subplot(424)
plt.title('Transitivity')
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.yscale('log')
plt.ylabel(r'$p$ : K-S Test')
plt.axis([1,200,1e-3,1])
plt.plot(ang,kstrmone,'g')
plt.plot(ang,kstrmtwo,'b')
plt.plot(ang,kstrrand,color='grey')
plt.fill_between(ang,thres,kstrmone,where=kstrmone>thres ,color='g',alpha=0.2)
plt.fill_between(ang,thres,kstrmtwo,where=kstrmtwo>thres ,color='b',alpha=0.2)
plt.fill_between(ang,thres,kstrrand,where=kstrrand>thres ,color='grey',alpha=0.2)
plt.axvline(70, color='grey', linestyle='--')
plt.text(75, 10**(-2.5), r'$70^{\prime\prime}$',color='grey', fontsize=13)
plt.text(75, 10**(-0.8), 'TR70',color='grey', fontsize=25)


plt.subplot(425)
plt.title('Edge Density')
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.yscale('log')
plt.ylabel(r'$p$ : K-S Test')
plt.axis([1,200,1e-3,1])
plt.plot(ang,ksdenmone,'g')
plt.plot(ang,ksdenmtwo,'b')
plt.plot(ang,ksdenrand,color='grey')
plt.fill_between(ang,thres,ksdenmone,where=ksdenmone>thres ,color='g',alpha=0.2)
plt.fill_between(ang,thres,ksdenmtwo,where=ksdenmtwo>thres ,color='b',alpha=0.2)
plt.fill_between(ang,thres,ksdenrand,where=ksdenrand>thres ,color='grey',alpha=0.2)

plt.axvline(100, color='grey', linestyle='--')
plt.text(105, 10**(-1), 'ED100',color='grey', fontsize=25)
plt.text(105, 10**(-2.5), r'$ >100^{\prime\prime}$',color='grey', fontsize=13)
#plt.legend(handles=[m2p,obsp,rp],loc=2)

plt.subplot(426)
plt.title('Size of the Largest Clique')
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
#plt.xlabel(r'Linking Length (arcsec)')
plt.yscale('log')
plt.ylabel(r'$p$ : K-S Test')
plt.axis([1,200,1e-3,1])
plt.plot(ang,kscqmone,'g')
plt.plot(ang,kscqmtwo,'b')
plt.plot(ang,kscqrand,color='grey')
plt.fill_between(ang,thres,kscqmone,where=kscqmone>thres ,color='g',alpha=0.2)
plt.fill_between(ang,thres,kscqmtwo,where=kscqmtwo>thres ,color='b',alpha=0.2)
plt.fill_between(ang,thres,kscqrand,where=kscqrand>thres ,color='grey',alpha=0.2)
plt.text(70, 0.05, r'$70^{\prime\prime}$',color='grey', fontsize=13)
plt.arrow(70, 0.05, 0, -0.022, head_width=5.0, head_length=0.02, fc='grey', ec='grey')

plt.subplot(427)
plt.title('Betweenness Centralization')
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.yscale('log')
plt.ylabel(r'$p$ : K-S Test')
plt.axis([1,200,1e-3,1])
plt.plot(ang,kscenbcmone,'g')
plt.plot(ang,kscenbcmtwo,'b')
plt.plot(ang,kscenbcrand,color='grey')
plt.fill_between(ang,thres,kscenbcmone,where=kscenbcmone>thres ,color='g',alpha=0.2)
plt.fill_between(ang,thres,kscenbcmtwo,where=kscenbcmtwo>thres ,color='b',alpha=0.2)
plt.fill_between(ang,thres,kscenbcrand,where=kscenbcrand>thres ,color='grey',alpha=0.2)



plt.subplot(428)
plt.title('Degree Centralization')
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.yscale('log')
plt.ylabel(r'$p$ : K-S Test')
plt.axis([1,200,1e-3,1])
plt.plot(ang,kscendgmone,'g')
plt.plot(ang,kscendgmtwo,'b')
plt.plot(ang,kscendgrand,color='grey')
plt.fill_between(ang,thres,kscendgmone,where=kscendgmone>thres ,color='g',alpha=0.2)
plt.fill_between(ang,thres,kscendgmtwo,where=kscendgmtwo>thres ,color='b',alpha=0.2)
plt.fill_between(ang,thres,kscendgrand,where=kscendgrand>thres ,color='grey',alpha=0.2)
plt.arrow(40, 0.05, 0, -0.022, head_width=5.0, head_length=0.02, fc='grey', ec='grey')
plt.text(40, 0.05, r'$40^{\prime\prime}$',color='grey', fontsize=13)

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
plt.savefig('kstest.pdf', dpi=600)
plt.show()
#%%

### plot network results on each separate pdf figure 

#plot the filter selections 
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman') 
plt.rcParams.update({'font.size': 14,'legend.fontsize': 14})
  
   
import matplotlib.patches as mp
obsp = mp.Patch(color='red', label='Observed LAEs')
m1p = mp.Patch(color='green', label='Model #1')
m2p = mp.Patch(color='blue', label='Model #2')
rp = mp.Patch(color='grey', label='Random')


fig = plt.figure(figsize=(5,3.5))

#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Diameter')
plt.axis([1,200,0,120])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,sdiabox1[30,:],color='green')
plt.fill_between(ang, sdiabox1[3,:],sdiabox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, sdiaboxobs[3,:],sdiaboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,sdiaboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, sdiaboxrand[3,:],sdiaboxrand[56,:], facecolor='grey', alpha=0.2)

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1diameter.pdf', dpi=600)
plt.show()

#%%

fig = plt.figure(figsize=(5,3.5))
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Giant Component Fraction')
plt.axis([1,200,0,1.1])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,sgcompbox1[30,:],color='green')
plt.fill_between(ang, sgcompbox1[3,:],sgcompbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, sgcompboxobs[3,:],sgcompboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,sgcompboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, sgcompboxrand[3,:],sgcompboxrand[56,:], facecolor='grey', alpha=0.2)

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1gcomp.pdf', dpi=600)
plt.show()

#%%

fig = plt.figure(figsize=(5,3.5))
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Average CC')
plt.axis([1,200,0.5,0.8])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=1)

#plt.plot(ang,savgccbox1[30,:],color='green')
plt.fill_between(ang, savgccbox1[3,:],savgccbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, savgccboxobs[3,:],savgccboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,savgccboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, savgccboxrand[3,:],savgccboxrand[56,:], facecolor='grey', alpha=0.2)
plt.axvline(170, color='grey', linestyle='--')
plt.text(170, 0.52, r'$ >170^{\prime\prime}$',color='grey', fontsize=13)
plt.text(105, 0.75, '\"CC170\"',color='green', fontsize=25)


plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1avgcc.pdf', dpi=600)
plt.show()

#%%

fig = plt.figure(figsize=(5,3.5))
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Transitivity')
plt.axis([1,200,0.5,0.8])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=1)
plt.axvline(70, color='grey', linestyle='--')
plt.text(75, 0.52, r'$70^{\prime\prime} = 1.42 h^{-1} $Mpc',color='grey', fontsize=20)
plt.text(75, 0.75, '\"TR70\"',color='green', fontsize=25)


#plt.plot(ang,strbox1[30,:],color='green')
plt.fill_between(ang, strbox1[3,:],strbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, strboxobs[3,:],strboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,strboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, strboxrand[3,:],strboxrand[56,:], facecolor='grey', alpha=0.2)

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1tr.pdf', dpi=600)
plt.show()

#%%

fig = plt.figure(figsize=(5,3.5))
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Edge Density')
plt.axis([1,200,0.0,0.01])
#plt.yscale('log')
#plt.xscale('log')
plt.axvline(100, color='grey', linestyle='--')
plt.text(105, 0.0015, '\"ED100\"',color='green', fontsize=25)
plt.text(105, 0.0003, r'$ >100^{\prime\prime}$',color='grey', fontsize=13)
plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,sdenbox1[30,:],color='green')
plt.fill_between(ang, sdenbox1[3,:],sdenbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, sdenboxobs[3,:],sdenboxobs[56,:], facecolor='red', alpha=0.3)
#plt.plot(ang,sdenboxobs[30,:],color='red',linewidth=0.5)
plt.fill_between(ang, sdenboxrand[3,:],sdenboxrand[56,:], facecolor='grey', alpha=0.2)

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1ed.pdf', dpi=600)
plt.show()

#%%

fig = plt.figure(figsize=(5,3.5))
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Size of the Largest Clique')
plt.axis([1,200,0,25])
#plt.yscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,smaxcliquebox1[30,:],color='green')
plt.fill_between(ang, smaxcliquebox1[3,:],smaxcliquebox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, smaxcliqueboxobs[3,:],smaxcliqueboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,smaxcliqueboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, smaxcliqueboxrand[3,:],smaxcliqueboxrand[56,:], facecolor='grey', alpha=0.2)
plt.text(75, 14, r'$70^{\prime\prime}$',color='grey', fontsize=20)
#plt.text(75, 2.3, r'Excess near $70^{\prime\prime}$',color='grey', fontsize=20)
plt.arrow(70, 16, 0, -2, head_width=5.0, head_length=1.0, fc='grey', ec='grey')

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1max.pdf', dpi=600)
plt.show()


#%%
fig = plt.figure(figsize=(5,3.5))
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Betweenness Centralization')
plt.axis([1,200,0.0,0.5])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,scenbcbox1[30,:],color='green')
plt.fill_between(ang, scenbcbox1[3,:],scenbcbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, scenbcboxobs[3,:],scenbcboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scenbcboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scenbcboxrand[3,:],scenbcboxrand[56,:], facecolor='grey', alpha=0.2)

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1bc.pdf', dpi=600)
plt.show()


#%%


fig = plt.figure(figsize=(5,3.5))
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Degree Centralization')
plt.axis([1,200,0.0,0.025])
#plt.yscale('log')
#plt.xscale('log')
#plt.legend(handles=[m1p,obsp,rp],loc=2)

#plt.plot(ang,scendgbox1[30,:],color='green')
plt.fill_between(ang, scendgbox1[3,:],scendgbox1[56,:], facecolor='green', alpha=0.7)
plt.fill_between(ang, scendgboxobs[3,:],scendgboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scendgboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scendgboxrand[3,:],scendgboxrand[56,:], facecolor='grey', alpha=0.2)
plt.text(45, 0.011, r'$40^{\prime\prime}$',color='grey', fontsize=20)
#plt.text(75, 2.3, r'Excess near $70^{\prime\prime}$',color='grey', fontsize=20)
plt.arrow(40, 0.013, 0, -0.0025, head_width=5.0, head_length=0.001, fc='grey', ec='grey')



plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.0)
plt.savefig('model1dc.pdf', dpi=600)
plt.show()

#%%


"""
## model 2 ....... 3x3 panel
fig = plt.figure(figsize=(17,13))


plt.subplot(331)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Diameter')
plt.axis([1,200,0,120])
#plt.yscale('log')
#plt.xscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=2)

#plt.plot(ang,sdiabox2[30,:],color='blue')
plt.fill_between(ang, sdiabox2[3,:],sdiabox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, sdiaboxobs[3,:],sdiaboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,sdiaboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, sdiaboxrand[3,:],sdiaboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(332)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Giant Component Fraction')
plt.axis([1,200,0,1.1])
#plt.yscale('log')
#plt.xscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=2)

#plt.plot(ang,sgcompbox2[30,:],color='blue')
plt.fill_between(ang, sgcompbox2[3,:],sgcompbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, sgcompboxobs[3,:],sgcompboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,sgcompboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, sgcompboxrand[3,:],sgcompboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(333)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Transitivity')
plt.axis([1,200,0.5,0.8])
#plt.yscale('log')
#plt.xscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=1)
plt.axvline(70, color='grey', linestyle='--')
plt.text(75, 0.52, r'$70^{\prime\prime}$',color='grey', fontsize=20)
plt.text(100, 0.52, '\"TR70\"',color='green', fontsize=25)


#plt.plot(ang,strbox2[30,:],color='blue')
plt.fill_between(ang, strbox2[3,:],strbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, strboxobs[3,:],strboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,strboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, strboxrand[3,:],strboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(334)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Average CC')
plt.axis([1,200,0.5,0.8])
#plt.yscale('log')
#plt.xscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=1)

#plt.plot(ang,savgccbox2[30,:],color='blue')
plt.fill_between(ang, savgccbox2[3,:],savgccbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, savgccboxobs[3,:],savgccboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,savgccboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, savgccboxrand[3,:],savgccboxrand[56,:], facecolor='grey', alpha=0.2)



plt.subplot(335)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Edge Density')
plt.axis([1,200,0.0,0.01])
#plt.yscale('log')
#plt.xscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=2)

#plt.plot(ang,sdenbox2[30,:],color='blue')
plt.fill_between(ang, sdenbox2[3,:],sdenbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, sdenboxobs[3,:],sdenboxobs[56,:], facecolor='red', alpha=0.3)
#plt.plot(ang,sdenboxobs[30,:],color='red',linewidth=0.5)
plt.fill_between(ang, sdenboxrand[3,:],sdenboxrand[56,:], facecolor='grey', alpha=0.2)



plt.subplot(336)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Size of the Largest Clique')
plt.axis([1,200,0,25])
#plt.yscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=2)

#plt.plot(ang,smaxcliquebox2[30,:],color='blue')
plt.fill_between(ang, smaxcliquebox2[3,:],smaxcliquebox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, smaxcliqueboxobs[3,:],smaxcliqueboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,smaxcliqueboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, smaxcliqueboxrand[3,:],smaxcliqueboxrand[56,:], facecolor='grey', alpha=0.2)
#plt.axvline(70, color='grey', linestyle='--')
plt.text(80, 14, r'$70^{\prime\prime}$',color='grey', fontsize=20)
#plt.text(75, 2.3, r'Excess near $70^{\prime\prime}$',color='grey', fontsize=20)
plt.arrow(75, 16, 0, -2, head_width=5.0, head_length=1.0, fc='grey', ec='grey')


plt.subplot(337)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Betweenness Centralization')
plt.axis([1,200,0.0,0.5])
#plt.yscale('log')
#plt.xscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=2)

#plt.plot(ang,scenbcbox2[30,:],color='blue')
plt.fill_between(ang, scenbcbox2[3,:],scenbcbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, scenbcboxobs[3,:],scenbcboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scenbcboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scenbcboxrand[3,:],scenbcboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(338)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Closeness Centralization')
plt.axis([1,200,0.00001,0.08])
plt.yscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=2)

#plt.plot(ang,scenclbox2[30,:],color='blue')
plt.fill_between(ang, scenclbox2[3,:],scenclbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, scenclboxobs[3,:],scenclboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scenclboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scenclboxrand[3,:],scenclboxrand[56,:], facecolor='grey', alpha=0.2)


plt.subplot(339)
#plt.scatter(r,r2xir,s=5.5,facecolor='none', edgecolor='red',marker="o")
plt.xlabel(r'Linking Length (arcsec)')
plt.ylabel(r'Degree Centralization')
plt.axis([1,200,0.0,0.025])
#plt.yscale('log')
#plt.xscale('log')
plt.legend(handles=[m2p,obsp,rp],loc=2)

#plt.plot(ang,scendgbox2[30,:],color='blue')
plt.fill_between(ang, scendgbox2[3,:],scendgbox2[56,:], facecolor='blue', alpha=0.7)
plt.fill_between(ang, scendgboxobs[3,:],scendgboxobs[56,:], facecolor='red', alpha=0.3)
plt.plot(ang,scendgboxobs[30,:],color='red',linewidth=3.0)
plt.fill_between(ang, scendgboxrand[3,:],scendgboxrand[56,:], facecolor='grey', alpha=0.2)
#plt.axvline(40, color='grey', linestyle='--')
#plt.text(45, 0.013, r'$40^{\prime\prime}$',color='grey', fontsize=20)
plt.text(45, 0.011, r'$40^{\prime\prime}$',color='grey', fontsize=20)
#plt.text(75, 2.3, r'Excess near $70^{\prime\prime}$',color='grey', fontsize=20)
plt.arrow(40, 0.013, 0, -0.0025, head_width=5.0, head_length=0.001, fc='grey', ec='grey')


plt.tight_layout(pad = -0.4, w_pad = -0.1, h_pad = -0.1)
plt.savefig('model2network.pdf', dpi=600)
plt.show()



print strboxobs[30,50:100]
print np.max(strboxobs[30,50:100])
print np.argmax(strboxobs[30,50:100])

print strboxobs[30,69]
"""





### clique section is obsolete, cuz the maximal clique is now in "netscalar.out"
"""
Generating edge files : R??arcsec??.edge from "dij.dat"

#####

shellscript = open("getClique.sh",'w')


cmdheader = "~/networkbin/netclique.out "
cmdtail = ""

catheader = "cat"
cmd =""
fedgename=""
curdir=""
curdijfile=""

#arcseclist = np.linspace(1,200,num=200,endpoint=True)
#arcseclist = np.linspace(1,5,num=5,endpoint=True)
arcseclist =[70,]


dirheader = "./geach12_model1/"
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    #print fedgename
    
    
    tmpstr = ''
    tmpfloat = 0.0
    for ridx, curarcsec in enumerate(arcseclist):
        curfname=str("R%04d.edge"%curarcsec)
        curnetname=str("R%04d.clique"%curarcsec)
 
        
        #make cmd 
        cmd = cmdheader+curdir+curfname+" "+curdir+curnetname
        print curdir+curfname
        print curdir+curnetname
        print cmd        
        shellscript.write(cmd+'\n')

dirheader = "./geach12_model2/"
for imodel in range(0,60):
    curdir = dirheader+catheader+str(imodel)+"/"
    #print fedgename
    
    
    tmpstr = ''
    tmpfloat = 0.0
    for ridx, curarcsec in enumerate(arcseclist):
        curfname=str("R%04d.edge"%curarcsec)
        curnetname=str("R%04d.clique"%curarcsec)
 
        
        #make cmd 
        cmd = cmdheader+curdir+curfname+" "+curdir+curnetname
        print curdir+curfname
        print curdir+curnetname
        print cmd        
        shellscript.write(cmd+'\n')
 
shellscript.close()



####### Parsing results
nummodel=60
arcseclist =[70,]

max1clique = np.zeros((nummodel,),dtype=np.int)
max1clique200 = np.zeros((nummodel,),dtype=np.int)

max2clique = np.zeros((nummodel,),dtype=np.int)


curfile=''
dirheader = "./geach12_model1/"
for imodel in range(0,nummodel):
    curdir = dirheader+catheader+str(imodel)+"/"
    #print fedgename
    
    
    for ridx, curarcsec in enumerate(arcseclist):
        curfile=curdir+str("R%04d.clique"%curarcsec)
        print curfile
 
        
        #read a networkfile list
        infile = open(curfile,"r")
        firstline = infile.readline()
        infile.close()
        print 'first line ', firstline
        print str(firstline).split(':')[1]

    max1clique[imodel] = np.int(str(firstline).split(':')[1])


print max1clique

curfile=''
dirheader = "./geach12_model2/"
for imodel in range(0,nummodel):
    curdir = dirheader+catheader+str(imodel)+"/"
    #print fedgename
    
    
    for ridx, curarcsec in enumerate(arcseclist):
        curfile=curdir+str("R%04d.clique"%curarcsec)
        print curfile
 
        
        #read a networkfile list
        infile = open(curfile,"r")
        firstline = infile.readline()
        infile.close()
        print 'first line ', firstline
        print str(firstline).split(':')[1]

    max2clique[imodel] = np.int(str(firstline).split(':')[1])


print max2clique



plt.axis([0,15,0,45])
plt.hist(max1clique,bins=range(0,15))
#plt.hist(max2clique,bins=range(0,15))
#plt.hist(maxclique200,bins=range(0,30))
plt.show()


"""

