import numpy as np
import networkx as nx
import tqdm
import percolate
import sys

f1=open('bond-mc.dat')
f11=open('datos_municipios-mc.dat')

qa=int(sys.argv[1])

conexiones=[]
maxcl=np.zeros(226)
maxcl2=np.zeros(226)
maxsuave=np.zeros(226)
max2suave=np.zeros(226)


for line in f1:
    x1=line.split()
    i, j= int(x1[0]), int(x1[1])
    conexiones.append([i,j])
print(len(conexiones))

nearest=[]

for line in f11:
    x1=line.split()
    tmpn=[]
    for x in x1:
        tmp=int(x)-1
        tmpn.append(tmp)
    nearest.append(tmpn)


nodes=np.array(range(226))


for i in tqdm.tqdm(range(10000)):
    xUF=nx.utils.UnionFind()
    ncluster=np.zeros(226)
    occ=np.zeros(226)
    np.random.shuffle(nodes)
    subbonds=conexiones[:qa]
    for j in range(226):
        itmp=nodes[j]
        occ[itmp]+=1
        tmp=[]
        for y in nearest[itmp]:
            if occ[y]==1:
                if [itmp,y] in subbonds or [y,itmp] in subbonds:
                    tmp.append(y)
        if len(tmp)==0:
            ncluster[itmp]+=1
        else:
            ntmp=1
            for y in tmp:
                l2=xUF[y]
                xUF.union(itmp,y)
                ntmp+=ncluster[l2]
                ncluster[l2]=0
            l1=xUF[itmp]
            ncluster[l1]+=ntmp
            
    
        maxcl[j]=(i*maxcl[j]+max(ncluster))/(i+1)
        maxcl2[j]=(i*maxcl2[j]+max(ncluster)**2)/(i+1)

for i in range(226):
    
    ptmp=(i+1)/226
    pbin=percolate.percolate._binomial_pmf(n=225,p=ptmp)
    smax=sum(pbin*maxcl)
    smax2=sum(pbin*maxcl2)
    maxsuave[i]+=smax
    max2suave[i]+=smax2
suscp=(max2suave-maxsuave**2)/(maxsuave*226)
itmp=0
atmp=0
for i in range(226):
    
    if atmp<suscp[i]:

        atmp=suscp[i]
        itmp=i
print(qa, itmp+1)
