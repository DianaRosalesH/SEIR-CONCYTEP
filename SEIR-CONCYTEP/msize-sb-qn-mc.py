import numpy as np
import networkx as nx
import tqdm
import percolate
import sys

f1=open('bond-mc.dat')

qn=int(sys.argv[1])

conexiones=np.zeros([662,2],int)
maxcl=np.zeros(662)
maxsuave=np.zeros(662)

f2=open('max-sb-qn-mc-'+str(qn)+'.dat', 'w')
f3=open('dmax-sb-qn-mc-'+str(qn)+'.dat', 'w')

c=0
for line in f1:
    x1=line.split()
    i, j= int(x1[0]), int(x1[1])
    conexiones[c][0]+=i
    conexiones[c][1]+=j
    c+=1


bonds=np.array(range(662))
nodes=np.array(range(226))
subnodes=nodes[:qn]

for i in tqdm.tqdm(range(100000)):
    xUF=nx.utils.UnionFind()
    np.random.shuffle(bonds)
    ncluster=np.ones(226)
    for j in range(662):
        itmp=bonds[j]
        N1=conexiones[itmp][0]
        N2=conexiones[itmp][1]
        if N1 in subnodes and N2 in subnodes:
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

for i in range(662):
    #print(i+1, maxcl[i], smean[i], file=f2)
    ptmp=(i+1)/662
    pbin=percolate.percolate._binomial_pmf(n=661,p=ptmp)
    smax=sum(pbin*maxcl)
    maxsuave[i]+=smax
    print(i+1, smax, file=f2)

for i in range(2,660):
    d=-maxsuave[i+2]+8*maxsuave[i+1]-8*maxsuave[i-1]+maxsuave[i-2]
    print(i+1, d, file=f3)