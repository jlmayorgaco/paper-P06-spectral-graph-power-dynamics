from geometry import *

def main():
    lock={'source':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'protocol':hashlib.sha256((OUT/'PROTOCOL.md').read_bytes()).hexdigest()}
    target=OUT/'RADIUS_SOURCE_LOCK.json'
    if target.exists(): assert json.loads(target.read_text())==lock
    else:save(target.name,lock)
    corners=[-.05+30.7j,.02+30.7j,.02+31.2j,-.05+31.2j,-.05+30.7j]
    points=[a+(b-a)*t/8 for a,b in zip(corners,corners[1:]) for t in range(8)]
    rows=[];summary=[]
    for name in ['joint','corrected']:
        d=json.loads((OLD/f'CERT_{name}.json').read_text());p=np.array(d['p']);m=Model(p)
        gs=[return_matrix(m,s) for s in points]
        for r in [0,.001,.005,.01,.025,.05]:
            lo=np.maximum(LO,p[10:]*(1-r));hi=np.minimum(HI,p[10:]*(1+r))
            witnesses=[]
            for j,(s,G) in enumerate(zip(points,gs)):
                Q,lab=inequalities(G,s,lo,hi);res=exclude(Q)
                witnesses.append({'point':j,'s':[s.real,s.imag],'labels':lab,**res})
                rows.append(dict(design=name,radius=r,point=j,slack=res['slack'],seconds=res['seconds'],status=res['status']))
            summary.append(dict(design=name,radius=r,min_slack=min(x['slack'] for x in witnesses),
                                positive_points=sum(x['slack']>1e-9 for x in witnesses),count=len(points)))
            save(f'RADIUS_{name}_{r}.json',{'lo':lo.tolist(),'hi':hi.tolist(),'rho':p[:10].tolist(),'witnesses':witnesses})
            pd.DataFrame(rows).to_csv(OUT/'TABLE_03_RADIUS_POINTS.csv',index=False)
            pd.DataFrame(summary).to_csv(OUT/'TABLE_04_RADIUS_SUMMARY.csv',index=False)
            print(summary[-1],flush=True)
if __name__=='__main__':main()
