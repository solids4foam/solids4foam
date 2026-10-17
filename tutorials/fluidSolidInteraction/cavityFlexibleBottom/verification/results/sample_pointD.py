import re,sys,numpy as np
def readvec(path):
    t=open(path).read()
    m=re.search(r'internalField\s+nonuniform\s+List<vector>\s*\n?(\d+)\s*\n\(',t)
    n=int(m.group(1)); s=m.end()
    body=t[s:]
    vals=re.findall(r'\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)',body[:n*80+n*10])[:n]
    a=np.array(vals,float); assert len(a)==n; return a
def samp(case,t='200'):
    P=readvec(case+'/constant/solid/polyMesh/points'.replace('internalField','')) if False else None
def readpts(path):
    t=open(path).read()
    m=re.search(r'\n(\d+)\s*\n\(',t); n=int(m.group(1))
    vals=re.findall(r'\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)',t[m.end():])[:n]
    return np.array(vals,float)
def at(case,x,y,t='200'):
    P=readpts(case+'/constant/solid/polyMesh/points')
    pD=readvec(case+'/%s/solid/pointD'%t)
    idx=np.where((abs(P[:,0]-x)<1e-9)&(abs(P[:,1]-y)<1e-9))[0]
    assert len(idx)==2,(len(idx))
    return pD[idx].mean(axis=0), pD[idx]
if __name__=='__main__':
    case=sys.argv[1]
    for y in (-1.0,-1.2,-1.4):
        m,a=at(case,4.0,y); print(y,m[0],m[1],a[:,1].tolist())
