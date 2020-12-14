import numpy as np
import networkx as nx
import tqdm
import percolate

f1=open('datos_municipios.dat')
f2=open('max2.dat', 'w')
f3=open('dmax.dat', 'w')

maxcl=np.zeros(226)
maxsuave=np.zeros(226)


nearest=[]

for line in f1:
    x1=line.split()
    tmpn=[]
    for x in x1:
        tmp=int(x)-1
        tmpn.append(tmp)
    nearest.append(tmpn)



nodes=np.array(range(226))

for i in tqdm.tqdm(range(100000)):
    c=0
    xUF=nx.utils.UnionFind()
    ncluster=np.zeros(226)
    occ=np.zeros(226)
    np.random.shuffle(nodes)
    for x in nodes:
        occ[x]+=1
        tmp=[]
        for y in nearest[x]:
            if occ[y]==1:
                tmp.append(y)
        if len(tmp)==0:
            ncluster[x]+=1
        else:
            ntmp=1
            for y in tmp:
                l2=xUF[y]
                xUF.union(x,y)
                ntmp+=ncluster[l2]
                ncluster[l2]=0
            l1=xUF[x]
            ncluster[l1]+=ntmp
            
    
        maxcl[c]=(i*maxcl[c]+max(ncluster))/(i+1)
#        s1=sum(ncluster)
#        s2=sum(ncluster**2)
#        stmp=s2/s1
#        stmp=0
#        ct=0
#        for i in range(226):
#            if ncluster[i]!=0:
#                stmp+=ncluster[i]
#                ct+=1
#        stmp=stmp/ct
        #smean[c]=(i*smean[c]+stmp)/(i+1)
        c+=1
for i in range(226):
    #print(i+1, maxcl[i], smean[i], file=f2)
    ptmp=(i+1)/226
    pbin=percolate.percolate._binomial_pmf(n=225,p=ptmp)
    smax=sum(pbin*maxcl)
    maxsuave[i]+=smax
    print(i+1, smax, file=f2)

for i in range(2,224):
    d=-maxsuave[i+2]+8*maxsuave[i+1]-8*maxsuave[i-1]+maxsuave[i-2]
    print(i+1, d, file=f3)