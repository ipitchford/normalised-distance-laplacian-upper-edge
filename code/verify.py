"""Finite independent checks for the normalised distance Laplacian paper.

Run with Python 3, NumPy, SciPy and SymPy. Large-n quotients use
only integer sums and fractions. They certify upper bounds on the deficit at
alpha=1; no finite calculation establishes an asymptotic lower bound.
"""
from fractions import Fraction
from math import ceil, comb, log, log1p, sqrt
from collections import deque
import numpy as np
import scipy.linalg as sla
from scipy.integrate import quad
import sympy as sy



def require(condition):
    """Fail loudly even under python -O (no reliance on assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED")

def exact_path_polynomial(n):
    distance = sy.Matrix([[abs(i-j) for j in range(n)] for i in range(n)])
    transmission = sy.diag(*[sum(distance.row(i)) for i in range(n)])
    transition = transmission.inv() * distance
    x = sy.Symbol('x')
    return sy.factor(transition.charpoly(x).as_expr())


def pair_sums(n, r, internal, weight=lambda d: d):
    """Trial E,H for profile r-2k. Integer weights give exact results."""
    left = (n - r - 1) // 2
    right = n - r - 1 - left
    require(min(left, right) >= 1)
    within_weight = weight(internal) * (comb(left, 2) + comb(right, 2))
    e = 4*r*r*within_weight
    h = 2*r*r*(within_weight + left*right*weight(r+2))
    for k in range(r+1):
        value = r - 2*k
        to_left = left*weight(k+1)
        to_right = right*weight(r-k+1)
        e += to_left*(value+r)**2 + to_right*(value-r)**2
        h += (to_left+to_right)*(r*r+value*value)

    # Group geodesic pairs by their separation g. Left indices run 0..m-1.
    for g in range(1, r+1):
        m = r+1-g
        sum_i = m*(m-1)//2
        sum_ii = m*(m-1)*(2*m-1)//6

        def sum_affine_squares(intercept, slope):
            return m*intercept**2 + 2*intercept*slope*sum_i + slope*slope*sum_ii

        e += weight(g)*sum_affine_squares(2*r-2*g, -4)
        h += weight(g)*(sum_affine_squares(r, -2) + sum_affine_squares(r-2*g, -2))
    return e, h


def construction(n, r, internal):
    """Join two cliques or stars to a path; return graph and integer profile."""
    left = (n-r-1)//2
    graph = [set() for _ in range(n)]

    def add_edge(u,v):
        graph[u].add(v)
        graph[v].add(u)

    for k in range(r):
        add_edge(k,k+1)
    profile = [r-2*k for k in range(r+1)]
    for root, bulk, value in [(0, range(r+1, r+1+left), r),
                              (r, range(r+1+left, n), -r)]:
        for u in bulk:
            add_edge(root,u)
        if internal == 1:
            for u in bulk:
                for v in range(u+1,bulk.stop):
                    add_edge(u,v)
        profile.extend([value]*len(bulk))
    return graph, profile


def brute_check(n, r, internal):
    graph, profile = construction(n,r,internal)
    distances = []
    for source in range(n):
        row = {source: 0}
        queue = deque([source])
        while queue:
            u = queue.popleft()
            for v in graph[u]:
                if v not in row:
                    row[v] = row[u]+1
                    queue.append(v)
        distances.append(row)
    e = sum(distances[i][j]*(profile[i]+profile[j])**2
            for i in range(n) for j in range(i+1,n))
    h = sum(sum(distances[i].values())*profile[i]**2 for i in range(n))
    require((e,h) == pair_sums(n,r,internal))
    d = np.array([[distances[i][j] for j in range(n)] for i in range(n)])
    t = d.sum(axis=1)
    least = sla.eigh(d,np.diag(t),subset_by_index=(0,0),eigvals_only=True)[0]
    require(1+least <= e/h+1e-12)


def power_check():
    """Numerical checks for general alpha, using a full distance matrix."""
    for alpha in (0.5, 1., 2., 3.):
        def integrand(z):
            a,b = z**alpha,(1-z)**alpha
            return 4*a*b/(a+b)

        coefficient = quad(integrand, 0., 1., epsabs=1e-12)[0]
        for kind in ("clique", "bipartite", "tree"):
            b = {"clique":1., "bipartite":(1+2**alpha)/2, "tree":2**alpha}[kind]
            s = (2*alpha*b/coefficient)**(1/(alpha+1))
            n = 211
            r = min(round(s*n**(1/(alpha+1))),n-6)
            left = (n-r-1)//2
            sides = [list(range(r+1,r+1+left)), list(range(r+1+left,n))]
            D = np.zeros((n,n))
            path = np.arange(r+1)
            D[:r+1,:r+1] = np.abs(path[:,None]-path[None,:])
            off = []
            for side_no, side in enumerate(sides):
                idx = np.arange(len(side))
                if kind == 'tree':
                    block = np.full((len(side),len(side)),2.)
                    delta = np.ones(len(side))
                elif kind == 'clique':
                    block = np.ones((len(side),len(side)))
                    delta = np.ones(len(side))
                else:
                    # Root is in part 0. Number ordinary vertices of each
                    # colour so that including the root balances the parts.
                    colours = (idx+1)%2
                    block = 1.+(colours[:,None]==colours[None,:])
                    delta = 1.+(colours==0)
                off.append(delta)
                np.fill_diagonal(block,0.)
                D[np.ix_(side,side)] = block
                to_path = delta[:,None] + (path if side_no==0 else r-path)
                D[np.ix_(side,path)] = to_path
                D[np.ix_(path,side)] = to_path.T
            cross = off[0][:,None]+r+off[1][None,:]
            D[np.ix_(*sides)] = cross
            D[np.ix_(*sides[::-1])] = cross.T
            W = D**alpha
            t = W.sum(axis=1)
            z = path/r
            profile = np.r_[( (1-z)**alpha-z**alpha )/((1-z)**alpha+z**alpha),
                            np.ones(len(sides[0])), -np.ones(len(sides[1]))]
            E = 0.5*np.sum(W*(profile[:,None]+profile[None,:])**2)
            ratio = E/np.dot(t,profile**2)
            adjoint = W/np.sqrt(t[:,None]*t[None,:])
            ev = sla.eigvalsh(adjoint)
            require(np.isclose(ev[-1],1.) and 0 <= 1+ev[0] <= ratio+1e-12)
            require(np.isclose((ev*ev).sum(), (adjoint*adjoint).sum()))
            print(f"{kind:9s} {alpha:3.1f} {r:3d} {1+ev[0]:.8f} {ratio:.8f}")


def kernel_checks():
    print("Slowly varying weight w(d)=log(1+d). The printed number is a scaled upper bound")
    print("from clique sizes, minimum across-clique weight and maximum clique-to-path weight.")
    print("log10(n), log(n) * upper bound, r=floor(n/ceil(log(n)^2))")
    m = log(2)
    for exponent in (4,16,64,256,1024,4096):
        n = 10**exponent
        r = n//ceil(log(n)**2)
        left = (n-r+1)//2
        right = n-r+1-left
        within = (comb(left,2)+comb(right,2))/(n*n)
        block_path_count = (left+right)*(r-1)/(n*n)
        cross_count = 2*left*right/(n*n)
        upper = (4*m*within+block_path_count*log(r+3))/(cross_count*log(r+1))
        print(exponent, f"{log(n)*upper:.9f}")
    print("predicted limit", f"{2*m:.9f}")
    print("Positive index w(d)=d log(1+d), exact pair counts with floating weights")
    print("n h r w(h)*trial, r=round(sqrt(3 log(2))*h)")
    weight = lambda d: d*log1p(d)
    for n in (1000,100000,10000000):
        lo,hi = 1,n
        while lo<hi:
            mid = (lo+hi+1)//2
            if mid*weight(mid)<=n:
                lo=mid
            else:
                hi=mid-1
        h = lo
        r = round(sqrt(3*m)*h)
        e,denom = pair_sums(n,r,1,weight=weight)
        print(n,h,r,f"{weight(h)*e/denom:.9f}")
    print("predicted limit", f"{4*sqrt(m/3):.9f}")
    rng = np.random.default_rng(5728391)
    for n in range(2,32):
        for m,M in ((1., 1.),(1.,4.),(0.4,12.)):
            f = lambda a: m*(a-1)/(m*(a-1)+M*(n-a))
            threshold = f(n//2)+f(n-n//2)
            for balanced in (False,True):
                if balanced:
                    group = np.arange(n)<(n//2)
                    W = np.where(group[:,None]==group[None,:],m,M)
                else:
                    upper = np.triu(rng.uniform(m,M,(n,n)),k=1)
                    W = upper+upper.T
                np.fill_diagonal(W,0.)
                t = W.sum(axis=1)
                eps = 1+sla.eigvalsh(W/np.sqrt(t[:,None]*t[None,:]),
                                    subset_by_index=(0,0))[0]
                require(eps >= threshold-1e-12)
                if balanced:
                    require(abs(eps-threshold)<1e-12)
    print("PASS finite bounded-weight comparison on 180 matrices, including all 90 balanced equality cases.")


def main():
    print("Exact characteristic polynomial det(xI-P) on paths (alpha=1)")
    x = sy.Symbol('x')
    for n in (2,3,4,5):
        poly = exact_path_polynomial(n)
        require(poly.subs(x,1)==0)
        print(n,poly)
    for n in range(7,20):
        for r in range(1,n-2):
            for internal in (1,2):
                brute_check(n,r,internal)
    print(f"Direct all-pairs integer energy, denominator, and numerical spectral inequality checked for {sum(2*(n-3) for n in range(7,20))} parameter triples.")
    print("Large-n exact trial upper bounds; final column is their rounded scaled value.")
    print("class n r E/H sqrt(n)*E/H")
    for internal,label in ((1,"cliques"),(2,"stars")):
        for n in (100,10000,1000000,100000000):
            r = round(sqrt(3*internal*n))
            trial = Fraction(*pair_sums(n,r,internal))
            print(label,n,r,trial,f"{sqrt(n)*float(trial):.9f}")
        print("predicted limit", f"{4*sqrt(internal/3):.9f}")
    print("Full numerical distance-power checks at n=211: class, power, r, true deficit, trial bound")
    power_check()
    kernel_checks()
    print("PASS")


if __name__ == '__main__':
    main()
