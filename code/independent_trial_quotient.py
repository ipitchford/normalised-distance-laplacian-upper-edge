"""Independent integer evaluation of manuscript trial quotient for alpha=1."""
from math import sqrt,comb

def trial(n,r,tree=False):
    m=(n-r+1)//2
    masses=(m,n-r+1-m)
    # Distances from the attachment vertex within the block are 0 and 1.
    E=H=0
    for size in masses:
        within=comb(size,2)+(comb(size-1,2) if tree else 0)
        E+=4*r*r*within
        H+=2*r*r*within
    l,m=masses
    W=l*m*r+(l-1)*m+(m-1)*l
    H+=2*r*r*W
    for i in range(1,r):
        y=r-2*i
        left=l*i+l-1
        right=m*(r-i)+m-1
        E+=left*(r+y)**2+right*(y-r)**2
        H+=(left+right)*(r*r+y*y)
        for j in range(1,i):
            v=r-2*j
            E+=(i-j)*(y+v)**2
            H+=(i-j)*(y*y+v*v)
    return E/H*sqrt(n)

if __name__=="__main__":
    for n in (256,4096,65536):
        for tree in (False,True):
            r=round(sqrt(n*(6 if tree else 3)))
            print(n, 'tree' if tree else 'clique', r, f'{trial(n,r,tree):.8f}')
