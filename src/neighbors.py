import numpy as np

def _norm(v):
    v=np.asarray(v,float); n=np.linalg.norm(v); return v/max(n,1e-12)

def rod_distance(a,b):
    p1=np.asarray(a['centroid']); p2=np.asarray(b['centroid'])
    u=_norm(a['axis']); v=_norm(b['axis'])
    h1=min(max(a.get('length',12)/2,4),12); h2=min(max(b.get('length',12)/2,4),12)
    x1=p1+np.linspace(-h1,h1,9)[:,None]*u; x2=p2+np.linspace(-h2,h2,9)[:,None]*v
    return min(float(np.min(np.linalg.norm(x2-p,axis=1))) for p in x1)

def graph(candidates, center_cutoff=12.0, rod_cutoff=7.0):
    g={i:set() for i in range(len(candidates))}; edges=[]
    for i in range(len(candidates)):
        for j in range(i+1,len(candidates)):
            cd=float(np.linalg.norm(np.asarray(candidates[i]['centroid'])-np.asarray(candidates[j]['centroid'])))
            rd=rod_distance(candidates[i],candidates[j])
            if cd<=center_cutoff or rd<=rod_cutoff:
                g[i].add(j); g[j].add(i); edges.append((i,j,cd,rd))
    return g,edges

def groups(g):
    rem=set(g); out=[]
    while rem:
        s=rem.pop(); stack=[s]; grp=[s]
        while stack:
            i=stack.pop()
            for j in g[i]:
                if j in rem: rem.remove(j); stack.append(j); grp.append(j)
        out.append(sorted(grp))
    return out

def atom_clash(a,b,cutoff=2.2):
    aa=np.asarray(a['coords']); bb=np.asarray(b['coords']); c2=cutoff*cutoff
    return any(np.any(np.sum((bb-p)**2,axis=1)<c2) for p in aa)

def pose_score(p):
    return float(p.get('v7_score',p.get('v6_score',p.get('score',0)))) + 1.0 + float(p.get('support_fraction',0)) - 2.0*float(p.get('foreign_fraction',0))

def optimize(group, pose_map, cutoff=2.2, beam_width=250):
    # Beam search works for both small and large groups and avoids exponential blow-up.
    states=[(0.0,{},[])]
    for site in group:
        new=[]
        for score,assign,selected in states:
            a=dict(assign); a[site]=None; new.append((score,a,list(selected)))
            for p in pose_map.get(site,[])[:5]:
                if any(atom_clash(p,q,cutoff) for q in selected): continue
                a=dict(assign); a[site]=p; s=list(selected); s.append(p)
                new.append((score+pose_score(p),a,s))
        new.sort(key=lambda x:x[0],reverse=True); states=new[:beam_width]
    if not states: return {},0.0,float('-inf')
    second=states[1][0] if len(states)>1 else float('-inf')
    return states[0][1],states[0][0],second
