import sys
inf = float('inf')
def bellman_ford(W,s):
    D,P = {v:inf for v in W},{}; D[s] = 0
    for _ in W:
        changed = False
        for u in W:
            for v in W[u]:
                if relax(W,u,v,D,P):
                    changed = True
        if not changed: break
    return D,P

def relax(W,u,v,D,P):
    d = D[u]+W[u][v]
    if d < D[v]:
        D[v],P[v] = d, u
        return True

def main():
    with open(sys.argv[1]) as file:
        n,m = map(int, file.readline().split())
        W = {v:{} for v in range(1,n+1)}
        for line in file:
            u,v,e = map(int, line.split()); W[u][v] = e
        D, P = bellman_ford(W,1)
        with open('output.txt','w') as outfile:
            for x in D.values():
                if x is inf: outfile.write('x ')
                else: outfile.write(str(x)+' ')

if __name__ == '__main__':
    main()
