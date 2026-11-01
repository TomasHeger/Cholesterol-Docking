import os, csv, math, random
from collections import deque
import numpy as np
from .diagnostics import Diagnostics
from . import neighbors
from .confidence import classify

DEFAULTS=dict(
 detection_level=.141, support_level=.05, orientations=1500, maxsites=60,
 max_tilt=45., mask_radius=4., protein_map_resolution=6., blob_smoothing=3.,
 blob_threshold_fraction=.10, protein_mask_radius=2.5, seed=12345,
 protein_exclusion=1.5, max_protein_distance=8., min_island_volume=8., max_island_volume=600.,
 island_group_distance=4.5, min_candidate_volume=25., max_candidate_volume=500.,
 mask_half_length=11., translation_along_axis=4., translation_across_axis=2.,
 pre_refine_keep=100, refine_cycles=6, translation_step=.8, rotation_step=6.,
 clash_distance=2., contact_distance=4.5, poses_per_site=5,
 strict_support=.80, strict_inside=.80, strict_foreign=.30,
 rescue_support=.70, rescue_inside=.70, rescue_foreign=.35,
 membrane_site_distance=8., minimum_outward_oh_cosine=.05,
 cholesterol_clash_distance=2.2,
)

def norm(v):
    v=np.asarray(v,float); return v/max(np.linalg.norm(v),1e-12)

def select_volume(models):
    x=[m for m in models if hasattr(m,'data') and hasattr(m,'interpolated_values')]
    if len(x)!=1: raise ValueError('Expected exactly one Volume; found %d'%len(x))
    return x[0]

def select_structure(models):
    x=[]
    for m in models:
        try:
            if hasattr(m,'atoms') and hasattr(m,'residues') and len(m.atoms)>0: x.append(m)
        except Exception: pass
    if len(x)!=1: raise ValueError('Expected exactly one AtomicStructure; found %d'%len(x))
    return x[0]

def heavy_atoms(model):
    out=[]
    for a in model.atoms:
        try:
            if a.element.name.upper()=='H': continue
        except Exception: pass
        out.append(a)
    return out

def atom_xyz(atoms):
    out=[]
    for a in atoms:
        try: p=a.scene_coord
        except Exception: p=a.coord
        out.append(np.asarray(p,float))
    return np.asarray(out)

class SpatialHash:
    def __init__(self,points,cell=4.):
        self.points=np.asarray(points,float); self.cell=float(cell); self.grid={}
        for i,p in enumerate(self.points): self.grid.setdefault(self.key(p),[]).append(i)
    def key(self,p): return tuple(np.floor(np.asarray(p)/self.cell).astype(int))
    def nearby(self,p,r):
        k=self.key(p); n=int(math.ceil(r/self.cell)); out=[]
        for a in range(-n,n+1):
            for b in range(-n,n+1):
                for c in range(-n,n+1): out.extend(self.grid.get((k[0]+a,k[1]+b,k[2]+c),[]))
        return out
    def mindist(self,p,r=20.):
        idx=self.nearby(p,r)
        if not idx:return float('inf')
        return float(np.min(np.linalg.norm(self.points[idx]-np.asarray(p),axis=1)))

def ijk_xyz(vol,ijk):
    p=np.asarray(vol.data.ijk_to_xyz(tuple(ijk)),float)
    try:p=np.asarray(vol.scene_position*tuple(p),float)
    except Exception:pass
    return p

def map_values(vol,scene_xyz):
    xyz=np.asarray(scene_xyz,float)
    try:
        inv=vol.scene_position.inverse(); xyz=np.asarray([inv*tuple(p) for p in xyz],np.float32)
    except Exception: xyz=np.asarray(xyz,np.float32)
    return np.asarray(vol.interpolated_values(xyz),float)

def pca(points,weights=None):
    x=np.asarray(points,float); w=np.ones(len(x)) if weights is None else np.asarray(weights,float)
    c=np.average(x,axis=0,weights=w); d=x-c; cov=(d*w[:,None]).T@d/max(w.sum(),1e-12)
    vals,vecs=np.linalg.eigh(cov); o=np.argsort(vals)[::-1]
    return c,vals[o],vecs[:,o]

def rotation_axis(axis,deg):
    x,y,z=norm(axis); a=math.radians(deg); c=math.cos(a); s=math.sin(a); C=1-c
    return np.array([[c+x*x*C,x*y*C-z*s,x*z*C+y*s],[y*x*C+z*s,c+y*y*C,y*z*C-x*s],[z*x*C-y*s,z*y*C+x*s,c+z*z*C]])

def rotation_vectors(a,b):
    a=norm(a); b=norm(b); v=np.cross(a,b); c=float(np.dot(a,b))
    if c>.999999:return np.eye(3)
    if c<-.999999:
        t=np.array([1.,0,0]) if abs(a[0])<.9 else np.array([0.,1,0]); return rotation_axis(np.cross(a,t),180)
    s=np.linalg.norm(v); vx=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    return np.eye(3)+vx+vx@vx*((1-c)/(s*s))

def transform(coords,R,t): return np.asarray(coords)@R.T+np.asarray(t)

def obtain_clr(session):
    before=set(session.models.list()); from chimerax.core.commands import run; run(session,'open ccd:CLR')
    new=list(set(session.models.list())-before); a=[m for m in new if hasattr(m,'atoms') and len(m.atoms)>0]
    if not a: raise RuntimeError('Could not fetch CCD:CLR')
    a[0].name='CLR docking template'
    try:a[0].display=False
    except Exception:pass
    return a[0]

def detect_candidates(session,vol,protein_xyz,p):
    mat=np.asarray(vol.data.matrix(),np.float32); nk,nj,ni=mat.shape; ph=SpatialHash(protein_xyz,4.)
    vox=set(); level=p['detection_level']
    for k in range(nk):
        for j in range(nj):
            for i in np.where(mat[k,j,:]>=level)[0]:
                xyz=ijk_xyz(vol,(i,j,k)); d=ph.mindist(xyz,p['max_protein_distance']+2)
                if p['protein_exclusion']<=d<=p['max_protein_distance']: vox.add((i,j,k))
    neigh=[(a,b,c) for a in(-1,0,1) for b in(-1,0,1) for c in(-1,0,1) if (a,b,c)!=(0,0,0)]
    rem=set(vox); comps=[]
    while rem:
        s=rem.pop(); q=deque([s]); comp=[s]
        while q:
            i,j,k=q.popleft()
            for a,b,c in neigh:
                x=(i+a,j+b,k+c)
                if x in rem: rem.remove(x); q.append(x); comp.append(x)
        comps.append(comp)
    vv=abs(np.prod(vol.data.step)); islands=[]
    for comp in comps:
        volume=len(comp)*vv
        if not(p['min_island_volume']<=volume<=p['max_island_volume']):continue
        xyz=np.asarray([ijk_xyz(vol,x) for x in comp]); weights=np.asarray([max(float(mat[k,j,i])-level,1e-6) for i,j,k in comp])
        islands.append(dict(xyz=xyz,centroid=np.average(xyz,axis=0,weights=weights),volume=volume,peak=max(float(mat[k,j,i]) for i,j,k in comp)))
    # conservative grouping of fragmented collinear islands
    unused=set(range(len(islands))); candidates=[]
    while unused:
        seed=unused.pop(); group=[seed]; changed=True
        while changed:
            changed=False; pts=np.concatenate([islands[i]['xyz'] for i in group]); center,vals,vecs=pca(pts); axis=vecs[:,0]; keep=set()
            for j in unused:
                d=islands[j]['centroid']-center; axial=abs(np.dot(d,axis)); radial=np.linalg.norm(d-np.dot(d,axis)*axis); direct=np.linalg.norm(d)
                if direct<=p['island_group_distance'] or (axial<=6 and radial<=2.5):
                    test=np.concatenate([islands[q]['xyz'] for q in group+[j]]); _,tv,_=pca(test)
                    if tv[0]>0 and tv[1]/tv[0]<=.70: group.append(j); changed=True
                    else: keep.add(j)
                else: keep.add(j)
            unused=keep
        volume=sum(islands[i]['volume'] for i in group)
        if not(p['min_candidate_volume']<=volume<=p['max_candidate_volume']):continue
        pts=np.concatenate([islands[i]['xyz'] for i in group]); center,vals,vecs=pca(pts); axis=norm(vecs[:,0]); proj=(pts-center)@axis
        radial=np.linalg.norm(pts-center-np.outer(proj,axis),axis=1)
        candidates.append(dict(centroid=center,axis=axis,eigenvalues=vals,elongation=math.sqrt(vals[0]/max(vals[1],1e-12)),length=float(proj.max()-proj.min()),radius=float(np.percentile(radial,90)),volume=volume,points=pts,n_islands=len(group),peak=max(islands[i]['peak'] for i in group)))
    candidates.sort(key=lambda c:c['volume']*c['peak'],reverse=True)
    return candidates[:p['maxsites']]

def protein_helix_normal(candidates):
    # Rod tensor is robust and was the dominant reliable signal in development.
    T=np.zeros((3,3))
    for c in candidates:
        if c['elongation']>=2 and c['length']>=10 and c['volume']>=60:
            a=norm(c['axis']); T+=math.sqrt(c['volume'])*c['elongation']*np.outer(a,a)
    if np.trace(T)<=0:
        for c in candidates: a=norm(c['axis']); T+=np.outer(a,a)
    vals,vecs=np.linalg.eigh(T); return norm(vecs[:,np.argmax(vals)])

def current_volumes(session): return [m for m in session.models.list() if hasattr(m,'data') and hasattr(m,'interpolated_values')]

def membrane_blob(session,expmap,protein,candidates,protein_xyz,p):
    from chimerax.core.commands import run
    before=set(current_volumes(session)); run(session,'molmap #%s %.2f gridSpacing %.4f'%(protein.id_string,p['protein_map_resolution'],min(expmap.data.step)))
    pm=[m for m in current_volumes(session) if m not in before][-1]; pm.name='V7 calculated protein density'
    before=set(current_volumes(session)); run(session,'volume subtract #%s #%s minRms true onGrid #%s valueType float32'%(expmap.id_string,pm.id_string,expmap.id_string))
    residual=[m for m in current_volumes(session) if m not in before][-1]; residual.name='V7 raw non-protein residual'
    before=set(current_volumes(session)); run(session,'volume gaussian #%s sDev %.2f'%(residual.id_string,p['blob_smoothing']))
    smooth=[m for m in current_volumes(session) if m not in before][-1]; smooth.name='V7 smoothed non-protein residual'
    mat=np.asarray(smooth.data.matrix(),np.float32); threshold=max(float(mat.max())*p['blob_threshold_fraction'],1e-8); step=2; ph=SpatialHash(protein_xyz,4.)
    cent=np.asarray([c['centroid'] for c in candidates]); bmin=cent.min(0)-18; bmax=cent.max(0)+18
    sampled=mat[::step,::step,::step]; binary=np.zeros(sampled.shape,bool)
    for sk,k in enumerate(range(0,mat.shape[0],step)):
        for sj,j in enumerate(range(0,mat.shape[1],step)):
            for si in np.where(sampled[sk,sj,:]>=threshold)[0]:
                i=int(si*step); xyz=ijk_xyz(smooth,(i,j,k))
                if np.any(xyz<bmin) or np.any(xyz>bmax):continue
                if ph.mindist(xyz,p['protein_mask_radius']+1)<p['protein_mask_radius']:continue
                binary[sk,sj,si]=True
    try:
        from scipy.ndimage import label,generate_binary_structure
        labels,n=label(binary,structure=generate_binary_structure(3,2))
    except Exception as e: raise RuntimeError('V7 membrane segmentation requires scipy.ndimage: %s'%e)
    vv=abs(np.prod(smooth.data.step))*step**3; comps=[]
    for lid in range(1,n+1):
        idx=np.argwhere(labels==lid); volume=len(idx)*vv
        if volume<1500:continue
        xyz=np.asarray([ijk_xyz(smooth,(int(si*step),int(sj*step),int(sk*step))) for sk,sj,si in idx]); weights=np.asarray([max(float(mat[int(sk*step),int(sj*step),int(si*step)]),1e-8) for sk,sj,si in idx])
        center,vals,vecs=pca(xyz,weights); test=xyz[::max(1,int(math.ceil(len(xyz)/15000)))]
        hits=sum(float(np.min(np.linalg.norm(test-c['centroid'],axis=1)))<=14 for c in candidates)
        flat=math.sqrt(max(vals[1],1e-8)/max(vals[2],1e-8)); comps.append(dict(indices=idx,xyz=xyz,weights=weights,center=center,values=vals,vectors=vecs,volume=volume,candidate_hits=hits,flatness=flat))
    if not comps:raise RuntimeError('No membrane residual component passed filtering')
    mv=max(c['volume'] for c in comps); mh=max(c['candidate_hits'] for c in comps)
    for c in comps:c['selection_score']=5*c['candidate_hits']/max(len(candidates),1)+2*c['candidate_hits']/max(mh,1)+.5*c['volume']/mv+.5*min(c['flatness']/3,1)
    blob=max(comps,key=lambda c:c['selection_score']); structural=protein_helix_normal(candidates); bn=norm(blob['vectors'][:,2]);
    if np.dot(bn,structural)<0:bn=-bn
    conf=np.clip((blob['flatness']-1)/1.5,0,1); normal=norm(.75*conf*bn+(1-.75*conf)*structural)
    proj=blob['xyz']@normal; w=blob['weights']
    def wp(q):
        o=np.argsort(proj); x=proj[o]; ww=w[o]; cum=np.cumsum(ww); return float(x[min(np.searchsorted(cum,q/100*cum[-1]),len(x)-1)])
    lo,hi=wp(7.5),wp(92.5); center=(lo+hi)/2
    # visible selected map
    try:
        from chimerax.map_data import ArrayGridData
        from chimerax.map import volume_from_grid_data
        out=np.zeros_like(mat)
        for sk,sj,si in blob['indices']:
            k,j,i=int(sk*step),int(sj*step),int(si*step); out[k:min(k+step,out.shape[0]),j:min(j+step,out.shape[1]),i:min(i+step,out.shape[2])]=mat[k,j,i]
        grid=ArrayGridData(out,origin=smooth.data.origin,step=smooth.data.step,cell_angles=smooth.data.cell_angles,rotation=smooth.data.rotation,name='V7 SELECTED MEMBRANE BLOB')
        bm=volume_from_grid_data(grid,session); bm.position=smooth.position
    except Exception as e: session.logger.warning('Could not create selected blob map: %s'%e)
    return dict(blob=blob,normal=normal,lower=lo,center=center,upper=hi,protein_map=pm,residual=residual,smooth=smooth)

def clr_template(session):
    m=obtain_clr(session); atoms=heavy_atoms(m); xyz=atom_xyz(atoms); center,vals,vecs=pca(xyz); axis=norm(vecs[:,0]); centered=xyz-center
    oxy=[i for i,a in enumerate(atoms) if getattr(a.element,'name','').upper()=='O']
    if len(oxy)!=1:raise RuntimeError('Expected exactly one CLR oxygen')
    return m,atoms,centered,axis,oxy[0]

def random_perp(axis,rng):
    for _ in range(10):
        v=np.array([rng.gauss(0,1),rng.gauss(0,1),rng.gauss(0,1)]); v-=np.dot(v,axis)*axis
        if np.linalg.norm(v)>1e-6:return norm(v)
    return norm(np.cross(axis,[1,0,0] if abs(axis[0])<.8 else [0,1,0]))

def guided_R(lig_axis,cand_axis,max_tilt,rng):
    target=np.array(cand_axis,float)*(1 if rng.random()>.5 else -1); perp=random_perp(target,rng); theta=math.degrees(math.acos(1-rng.random()*(1-math.cos(math.radians(max_tilt)))))
    target=rotation_axis(perp,theta)@target; return rotation_axis(target,rng.random()*360)@rotation_vectors(lig_axis,target)

def mask_fraction(coords,c):
    d=np.asarray(coords)-c['centroid']; ax=d@c['axis']; rad=np.linalg.norm(d-np.outer(ax,c['axis']),axis=1); return float(np.mean((np.abs(ax)<=11)&(rad<=4.0)))

def evaluate(coords,R,lig_axis,vol,c,ph,p):
    vals=map_values(vol,coords)
    if len(vals)==0 or not np.all(np.isfinite(vals)):return None
    support=float(np.mean(vals>=p['support_level'])); inside=mask_fraction(coords,c); clashes=0; contacts=0
    for q in coords:
        idx=ph.nearby(q,p['contact_distance'])
        if idx:
            md=float(np.min(np.linalg.norm(ph.points[idx]-q,axis=1))); clashes+=md<p['clash_distance']; contacts+=md<p['contact_distance']
    density=float(np.mean(np.clip(vals/max(p['support_level'],1e-8),-2,4)))
    fitted=norm(R@lig_axis); angle=math.degrees(math.acos(np.clip(abs(np.dot(fitted,c['axis'])),-1,1)))
    score=density+2*support-2*(1-support)-3*(1-inside)-5*clashes-.3*(angle/max(p['max_tilt'],1))**2
    return dict(score=score,mean_density=float(np.mean(vals)),support_fraction=support,inside_mask_fraction=inside,clashes=int(clashes),contacts=int(contacts),axis_angle=angle)

def refine(centered,lig_axis,R,t,scorefun,p):
    curR=R.copy(); curt=np.array(t,float); cur=scorefun(transform(centered,curR,curt),curR)
    if cur is None:return None
    axes=np.eye(3); ts=p['translation_step']; rs=p['rotation_step']
    for _ in range(p['refine_cycles']):
        improved=True; loops=0
        while improved and loops<20:
            improved=False; loops+=1
            for a in axes:
                for s in(-1,1):
                    nt=curt+s*ts*a; r=scorefun(transform(centered,curR,nt),curR)
                    if r and r['score']>cur['score']:curt=nt;cur=r;improved=True
            for a in axes:
                for s in(-1,1):
                    nr=rotation_axis(a,s*rs)@curR; r=scorefun(transform(centered,nr,curt),nr)
                    if r and r['score']>cur['score']:curR=nr;cur=r;improved=True
        ts*=.55;rs*=.55
    return dict(R=curR,t=curt,coords=transform(centered,curR,curt),**cur)

def foreign_fraction(pose,own,candidates):
    n=0
    for q in pose['coords']:
        invaded=False
        for j,c in enumerate(candidates):
            if j==own:continue
            d=q-c['centroid']; ax=np.dot(d,c['axis']); rad=np.linalg.norm(d-ax*c['axis'])
            if abs(ax)<=10 and rad<=4.5:invaded=True;break
        n+=invaded
    return n/len(pose['coords'])

def oh_cos(coords,oi,aqu): return float(np.dot(norm(coords[oi]-np.mean(coords,axis=0)),norm(aqu)))

def flip_R(pose,lig_axis):
    fa=norm(pose['R']@lig_axis); test=np.array([1.,0,0]) if abs(fa[0])<.85 else np.array([0.,1,0]); return rotation_axis(np.cross(fa,test),180)@pose['R']

def dock_site(session,site,c,allc,vol,ph,centered,lig_axis,oi,aqu,p,rng,diag,rescue=False):
    n=p['orientations']*(2 if rescue else 1); diag.generated(site,n,rescue)
    scorefun=lambda coords,R:evaluate(coords,R,lig_axis,vol,c,ph,p)
    initial=[]
    for _ in range(n):
        R=guided_R(lig_axis,c['axis'],p['max_tilt'],rng); along=rng.uniform(-4,4); rp=random_perp(c['axis'],rng); t=c['centroid']+along*c['axis']+rng.random()**.5*(2 if not rescue else 3)*rp
        r=scorefun(transform(centered,R,t),R)
        if r and r['support_fraction']>=.2:initial.append(dict(R=R,t=t,coords=transform(centered,R,t),**r))
    initial.sort(key=lambda x:x['score'],reverse=True); initial=initial[:p['pre_refine_keep']]; diag.refined(site,len(initial),rescue)
    out=[]
    for base in initial:
        choices=[]
        for label,R in [('ORIGINAL',base['R']),('FLIPPED',flip_R(base,lig_axis))]:
            rr=refine(centered,lig_axis,R,base['t'],scorefun,p)
            if rr:
                rr['orientation_source']=label; rr['oh_cosine']=oh_cos(rr['coords'],oi,aqu); choices.append(rr)
        if not choices:continue
        outward=[x for x in choices if x['oh_cosine']>.05]
        chosen=max(outward,key=lambda x:x['score']+1.5*x['oh_cosine']) if outward else max(choices,key=lambda x:x['score'])
        chosen['foreign_fraction']=foreign_fraction(chosen,site-1,allc); chosen['v6_score']=chosen['score']+1.5*chosen['oh_cosine']-1.5*chosen['foreign_fraction']; chosen['rescue_pass']=rescue
        # detailed hard filters
        reason=None
        if chosen['clashes']>0:reason='PROTEIN_CLASH'
        elif chosen['support_fraction']<(p['rescue_support'] if rescue else p['strict_support']):reason='LOW_SUPPORT'
        elif chosen['inside_mask_fraction']<(p['rescue_inside'] if rescue else p['strict_inside']):reason='OUTSIDE_MASK'
        elif chosen['foreign_fraction']>(p['rescue_foreign'] if rescue else p['strict_foreign']):reason='NEIGHBOUR_INTRUSION'
        elif chosen['oh_cosine']<=.05:reason='INWARD_OH'
        if reason:diag.reject(site,reason,chosen)
        else:out.append(chosen)
    out.sort(key=lambda x:x['v6_score'],reverse=True)
    # cluster roughly by RMSD/translation
    selected=[]
    for x in out:
        if any(np.linalg.norm(x['t']-y['t'])<2 and np.sqrt(np.mean(np.sum((x['coords']-y['coords'])**2,axis=1)))<2 for y in selected):continue
        selected.append(x)
        if len(selected)>=5:break
    return selected

def write_pdb(path,items,atoms):
    with open(path,'w') as f:
        serial=1
        for resid,item in enumerate(items,1):
            for atom,xyz in zip(atoms,item['pose']['coords']):
                elem=getattr(atom.element,'name','C'); name=atom.name[:4]; x,y,z=xyz
                f.write('HETATM%5d %-4s CLR Z%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%(serial,name,resid,x,y,z,1.,20.,elem));serial+=1
        f.write('END\n')

def run_docking(session,map_models,protein_models,out,**kwargs):
    p=dict(DEFAULTS); p.update({k:v for k,v in kwargs.items() if v is not None}); rng=random.Random(int(p['seed'])); diag=Diagnostics()
    vol=select_volume(map_models); protein=select_structure(protein_models); os.makedirs(out,exist_ok=True)
    patoms=heavy_atoms(protein); pxyz=atom_xyz(patoms); ph=SpatialHash(pxyz,4.5); clr,latoms,centered,laxis,oi=clr_template(session)
    candidates=detect_candidates(session,vol,pxyz,p); session.logger.info('V7 high-resolution candidate rods: %d'%len(candidates))
    mem=membrane_blob(session,vol,protein,candidates,pxyz,p); normal=mem['normal']; center=mem['center']; session.logger.info('V7 membrane interfaces %.2f / %.2f; thickness %.2f A'%(mem['lower'],mem['upper'],mem['upper']-mem['lower']))
    site_sets=[]
    for i,c in enumerate(candidates,1):
        proj=float(np.dot(c['centroid'],normal)); leaflet='UPPER' if proj>=center else 'LOWER'; aqu=normal if leaflet=='UPPER' else -normal
        # blob association
        pts=mem['blob']['xyz'][::max(1,int(math.ceil(len(mem['blob']['xyz'])/20000)))]; bd=float(np.min(np.linalg.norm(pts-c['centroid'],axis=1))); assoc=bd<=p['membrane_site_distance']
        poses=dock_site(session,i,c,candidates,vol,ph,centered,laxis,oi,aqu,p,rng,diag,False)
        if not poses:
            poses=dock_site(session,i,c,candidates,vol,ph,centered,laxis,oi,aqu,p,rng,diag,True)
        for q in poses:q.update(site=i,leaflet=leaflet,membrane_associated=assoc,blob_distance=bd)
        site_sets.append(dict(site=i,valid=poses,candidate=c,leaflet=leaflet,associated=assoc))
    # joint neighboring optimization
    g,edges=neighbors.graph(candidates); grps=neighbors.groups(g); posemap={s['site']-1:s['valid'] for s in site_sets}; accepted=[]; group_rows=[]
    for gid,grp in enumerate(grps,1):
        assignment,best,second=neighbors.optimize(grp,posemap,p['cholesterol_clash_distance']); gap=best-second if np.isfinite(second) else float('inf')
        for si in grp:
            pose=assignment.get(si); group_rows.append((gid,si+1,'OCCUPIED' if pose else 'EMPTY',best,second,gap))
            if pose:
                conf,why=classify(pose,gap); accepted.append(dict(site=si+1,pose=pose,confidence=conf,confidence_reason=why,group=gid,gap=gap)); diag.final(si+1,'ACCEPT')
            else:
                # distinguish no pose from competition
                if posemap.get(si): diag.final(si+1,'JOINT_OPTIMIZATION_EMPTY'); diag.reject(si+1,'JOINT_OPTIMIZATION_EMPTY',posemap[si][0])
                else: diag.final(si+1,'NO_VALID_POSE')
    accepted.sort(key=lambda x:x['site'])
    # outputs
    write_pdb(os.path.join(out,'all_accepted_cholesterols.pdb'),accepted,latoms)
    for level in ('HIGH','MEDIUM','LOW'):
        subset=[x for x in accepted if x['confidence']==level]
        if subset:write_pdb(os.path.join(out,'%s_confidence_cholesterols.pdb'%level.lower()),subset,latoms)
    with open(os.path.join(out,'final_selection.csv'),'w',newline='') as f:
        w=csv.writer(f);w.writerow(['site','status','confidence','confidence_reason','neighbor_group','group_score_gap','v6_score','density_score','leaflet','orientation','OH_cosine','support_fraction','inside_mask_fraction','foreign_fraction','protein_clashes','rescue_pass','membrane_associated','blob_distance_A'])
        amap={x['site']:x for x in accepted}
        for i in range(1,len(candidates)+1):
            if i in amap:
                x=amap[i];q=x['pose'];w.writerow([i,'ACCEPT',x['confidence'],x['confidence_reason'],x['group'],x['gap'],q.get('v6_score',''),q.get('score',''),q.get('leaflet',''),q.get('orientation_source',''),q.get('oh_cosine',''),q.get('support_fraction',''),q.get('inside_mask_fraction',''),q.get('foreign_fraction',''),q.get('clashes',''),q.get('rescue_pass',''),q.get('membrane_associated',''),q.get('blob_distance','')])
            else:w.writerow([i,'UNFILLED','','','','','','','','','','','','','','','',''])
    with open(os.path.join(out,'neighbor_groups.csv'),'w',newline='') as f:
        w=csv.writer(f);w.writerow(['group','site','selection','group_best_score','group_second_score','group_score_gap']);w.writerows(group_rows)
    with open(os.path.join(out,'site_diagnostics.csv'),'w',newline='') as f:
        w=csv.writer(f);w.writerow(['site','generated','refined','rescue_generated','rescue_refined','protein_clash','low_support','outside_mask','inward_OH','neighbor_intrusion','joint_optimization_empty','final_reason','best_rejected_score','best_rejected_support','best_rejected_OH','best_rejected_clashes'])
        for i in range(1,len(candidates)+1):
            d=diag.sites[i];br=d['best_rejected'][1] if d['best_rejected'] else {};w.writerow([i,d['generated'],d['refined'],d['rescue_generated'],d['rescue_refined'],d['reasons']['PROTEIN_CLASH'],d['reasons']['LOW_SUPPORT'],d['reasons']['OUTSIDE_MASK'],d['reasons']['INWARD_OH'],d['reasons']['NEIGHBOUR_INTRUSION'],d['reasons']['JOINT_OPTIMIZATION_EMPTY'],d['final_reason'],br.get('v6_score',br.get('score','')),br.get('support_fraction',''),br.get('oh_cosine',''),br.get('clashes','')])
    with open(os.path.join(out,'membrane_model.csv'),'w',newline='') as f:
        w=csv.writer(f);w.writerow(['normal_x','normal_y','normal_z','lower_interface','center','upper_interface','thickness_A','blob_volume_A3','blob_candidate_hits']);w.writerow([*normal,mem['lower'],center,mem['upper'],mem['upper']-mem['lower'],mem['blob']['volume'],mem['blob']['candidate_hits']])
    # open final
    from chimerax.core.commands import run
    final=os.path.join(out,'all_accepted_cholesterols.pdb')
    if accepted:run(session,'open "%s"'%final.replace('\\','/'))
    session.logger.info('V7 complete: %d accepted; HIGH %d, MEDIUM %d, LOW %d'%(len(accepted),sum(x['confidence']=='HIGH' for x in accepted),sum(x['confidence']=='MEDIUM' for x in accepted),sum(x['confidence']=='LOW' for x in accepted)))
    return dict(accepted=accepted,candidates=candidates,diagnostics=diag,output=out,membrane=mem,groups=grps)
