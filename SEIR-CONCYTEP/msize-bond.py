import numpy as np
import networkx as nx
import tqdm
import percolate

f1=open('bond.dat')

conexiones=np.zeros([575,2],int)
maxcl=np.zeros(575)
maxcl2=np.zeros(575)
maxsuave=np.zeros(575)
max2suave=np.zeros(575)

f2=open('max-bond.dat', 'w')
f3=open('dmax-bond.dat', 'w')

c=0
for line in f1:
    x1=line.split()
    i, j= int(x1[0]), int(x1[1])
    conexiones[c][0]+=i
    conexiones[c][1]+=j
    c+=1


bonds=np.array(range(575))

for i in tqdm.tqdm(range(100000)):
    xUF=nx.utils.UnionFind()
    added=np.zeros(226)
    np.random.shuffle(bonds)
    ncluster=np.ones(226)
    for j in range(575):
        itmp=bonds[j]
        N1=conexiones[itmp][0]
        N2=conexiones[itmp][1]
        added[N1]+=1
        added[N2]+=1
        l1=xUF[N1]
        l2=xUF[N2]
        xUF.union(N1, N2)
        if l1!=l2:
            cntmp=ncluster[l1]+ncluster[l2]
            ncluster[l1]=0
            ncluster[l2]=0
            l3=xUF[N1]
            ncluster[l3]+=cntmp
        maxcl[j]=(i*maxcl[j]+max(ncluster))/(i+1)
        maxcl2[j]=(i*maxcl2[j]+max(ncluster)**2)/(i+1)



for i in range(575):
    #print(i+1, maxcl[i], smean[i], file=f2)
    ptmp=(i+1)/575
    pbin=percolate.percolate._binomial_pmf(n=574,p=ptmp)
    smax=sum(pbin*maxcl)
    smax2=sum(pbin*maxcl2)
    maxsuave[i]+=smax
    max2suave[i]+=smax2

suscp=(max2suave-maxsuave**2)/(maxsuave*226)
itmp=0
atmp=0
for i in range(575):
    if atmp<suscp[i]:
        itmp=i+1
        atmp=suscp[i]
    print(i+1, maxsuave[i]/226, suscp[i], file=f2)

print(itmp)

for i in range(2,573):
    d=-maxsuave[i+2]+8*maxsuave[i+1]-8*maxsuave[i-1]+maxsuave[i-2]
    print(i+1, d, file=f3)
