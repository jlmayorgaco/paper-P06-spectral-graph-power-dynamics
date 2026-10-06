from certify_roots import *

def rectangle_count(bm,box,n):
    x0,x1,y0,y1=[arb(float(x)) for x in box]
    corners=[acb(x0,y0),acb(x1,y0),acb(x1,y1),acb(x0,y1),acb(x0,y0)]
    total=acb(0);rows=[]
    for a,b in zip(corners,corners[1:]):
        for k in range(n):
            aa=a+(b-a)*k/n;bb=a+(b-a)*(k+1)/n;ds=bb-aa;mid=(aa+bb)/2
            s=mid+acb(arb(0,(abs(ds.real)/2).upper()),arb(0,(abs(ds.imag)/2).upper()))
            K,Ks=bm.matrix(s);R=acb_mat([[K[i,j].mid() for j in range(10)] for i in range(10)]).inv()
            W=(R*K).inv()*(R*Ks);val=sum((W[i,i] for i in range(10)),acb(0))*ds
            total+=val;rows.append(bounds(val))
    count=total/(2*arb.pi()*acb(0,1))
    ok=bool(count.contains(acb(1)) and count.real>arb('.5') and count.real<arb('1.5'))
    return count,ok,rows

def main():
    box=(-.105,.01,30.40,31.05);center=(box[0]+box[1])/2+1j*(box[2]+box[3])/2
    rows=[]
    for name in ['anchor','single30','single37','joint']:
        result=json.loads((OUT/('CERT_'+name+'.json')).read_text());p=np.array(result['p'])
        bm=BallModel(p);bm.prepare(center,.35,order=14);start=time.perf_counter()
        for n in [64,128,256,512,1024]:
            try:
                count,ok,panels=rectangle_count(bm,box,n)
                print(name,n,count,ok,flush=True)
            except (ValueError,ZeroDivisionError) as ex:
                print(name,n,'interval inverse inconclusive',flush=True);continue
            if ok:break
        else:raise RuntimeError('common contour unresolved')
        row={'name':name,'box':box,'count':bounds(count),'one_root':ok,'panels':4*n,'seconds':time.perf_counter()-start}
        rows.append(row);save('COMMON_CONTOUR_CERTIFICATES.json',rows);save('COMMON_PANELS_'+name+'.json',panels)
    print('COMMON_COUNT_GATE_PASSED',flush=True)
if __name__=='__main__':main()
