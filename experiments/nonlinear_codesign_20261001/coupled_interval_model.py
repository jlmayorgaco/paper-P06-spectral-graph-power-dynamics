"""Full SG/GFL AC DAE in an exact rotating chart, with interval first derivatives.

Only common angular position is removed. Common physical frequency, every device
state and 39 phase-frequency/RoCoF instruments remain. The initial implementation
uses interval enclosures over boxes, which may be very conservative.
"""
from pathlib import Path
import math
import tomllib
import numpy as np
import mpmath as mp

mp.iv.dps = 60
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/nonlinear_codesign_20261001/coupled_network_contract"
EPS = np.finfo(float).eps
down = lambda x: np.nextafter(x, -np.inf)
up = lambda x: np.nextafter(x, np.inf)


def add(a, b):
    return down(a[0] + b[0]), up(a[1] + b[1])


def neg(a):
    return -a[1], -a[0]


def mul(a, b):
    terms = np.stack(np.broadcast_arrays(a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]))
    return down(terms.min(axis=0)), up(terms.max(axis=0))


def inv(a):
    if np.any((a[0] <= 0) & (a[1] >= 0)):
        raise ValueError("interval division across zero")
    return down(1/a[1]), up(1/a[0])


def magnitude(a):
    return np.maximum(np.abs(a[0]), np.abs(a[1]))


def positive_product(a, b):
    """Upper bound for products of nonnegative finite binary64 matrices.

    gamma uses 4*n*eps, exceeding the standard 2*n-operation dot-product bound;
    nextafter covers the final scale. Underflow is covered by a tiny additive term.
    """
    assert np.all(a >= 0) and np.all(b >= 0)
    n = a.shape[-1]
    gam = up((4*n*EPS)/(1-4*n*EPS))
    result = up((a @ b) * up(1+gam) + n*np.nextafter(0., 1.))
    assert np.isfinite(result).all()
    return result


def mm(a, b):
    am = a[0]/2 + a[1]/2
    bm = b[0]/2 + b[1]/2
    ar = up(np.maximum(am-a[0], a[1]-am))
    br = up(np.maximum(bm-b[0], b[1]-bm))
    mid = am @ bm
    rad = up(positive_product(np.abs(am), br) + positive_product(ar, np.abs(bm)))
    rad = up(rad + positive_product(ar, br))
    gam = up((4*am.shape[-1]*EPS)/(1-4*am.shape[-1]*EPS))
    rad = up(rad + up(gam*positive_product(np.abs(am), np.abs(bm))) + np.nextafter(0., 1.))
    return down(mid-rad), up(mid+rad)


def inverse_enclosure(g):
    n = len(g[0]); ident = np.eye(n)
    k = np.linalg.inv(g[0]/2 + g[1]/2)
    residual = add((ident, ident), neg(mm((k, k), g)))
    e = magnitude(residual)
    rows = up(np.sum(e, axis=1) * (1+4*n*EPS))
    if rows.max() >= 1:
        raise ValueError(f"inverse Neumann check failed: {rows.max()}")
    rhs = positive_product(e, np.abs(k))
    radius = np.maximum(np.linalg.solve(ident-e, rhs), 0.)
    radius = up(radius*1.00001 + 1e-22)
    for _ in range(8):
        bound = positive_product(e, up(np.abs(k)+radius))
        if np.all(bound <= radius):
            return (down(k-radius), up(k+radius)), float(rows.max())
        radius = up(np.maximum(radius, bound)*1.00001 + 1e-22)
    raise ValueError("inverse enclosure fixed-point inequality failed")


def ivfun(name, a):
    x = mp.iv.mpf([float(a[0]), float(a[1])])
    r = mp.iv.atan2(x, mp.iv.mpf(1)) if name == "atan" else getattr(mp.iv, name)(x)
    return down(float(r.a)), up(float(r.b))


class Jet:
    n = 0

    def __init__(self, lo, hi=None, dl=None, dh=None):
        self.v = (float(lo), float(lo if hi is None else hi))
        self.d = (np.zeros(self.n), np.zeros(self.n)) if dl is None else (dl, dh)

    @classmethod
    def from_parts(cls, v, d):
        return cls(v[0], v[1], d[0], d[1])

    @staticmethod
    def cast(x):
        return x if isinstance(x, Jet) else Jet(x)

    def __add__(self, other):
        b = self.cast(other)
        return self.from_parts(add(self.v, b.v), add(self.d, b.d))
    __radd__ = __add__

    def __neg__(self):
        return self.from_parts(neg(self.v), neg(self.d))

    def __sub__(self, other):
        return self + (-self.cast(other))

    def __rsub__(self, other):
        return self.cast(other) + (-self)

    def __mul__(self, other):
        b = self.cast(other)
        return self.from_parts(mul(self.v, b.v), add(mul(self.v, b.d), mul(b.v, self.d)))
    __rmul__ = __mul__

    def reciprocal(self):
        val = inv(self.v)
        der = neg(mul(val, val))
        return self.from_parts(val, mul(der, self.d))

    def __truediv__(self, other):
        return self * self.cast(other).reciprocal()

    def __rtruediv__(self, other):
        return self.cast(other) * self.reciprocal()

    def __pow__(self, n):
        if n != 2:
            raise ValueError("only square used by this model")
        val = mul(self.v, self.v)
        if self.v[0] <= 0 <= self.v[1]:
            val = (0., val[1])
        return self.from_parts(val, mul(mul((2., 2.), self.v), self.d))

    def fun(self, name):
        val = ivfun(name, self.v)
        if name == "sin":
            der = ivfun("cos", self.v)
        elif name == "cos":
            der = neg(ivfun("sin", self.v))
        elif name == "sqrt":
            der = mul((.5, .5), inv(val))
        elif name == "atan":
            der = inv(add((1., 1.), (self**2).v))
        else:
            raise ValueError(name)
        return self.from_parts(val, mul(der, self.d))


def sf(name, x):
    if hasattr(x, "_mpi_"):
        return mp.iv.atan2(x,mp.iv.mpf(1)) if name=="atan" else getattr(mp.iv,name)(x)
    return x.fun(name) if isinstance(x, Jet) else getattr(math, name)(x)


def lower(x):
    if hasattr(x, "_mpi_"): return down(float(x.a))
    return x.v[0] if isinstance(x, Jet) else float(x)


def upper(x):
    if hasattr(x, "_mpi_"): return up(float(x.b))
    return x.v[1] if isinstance(x, Jet) else float(x)


class Model:
    def __init__(self, path=OUT / "model.toml", center=None):
        self.raw = tomllib.loads(path.read_text(encoding="utf-8"))
        d = self.raw
        self.nx = d["state_count"]; self.nw = d["input_count"]
        self.x0 = np.array(d["x0"]); self.Y = np.array(d["Y"])
        self.sg = [np.array(x)-1 for x in d["sg_indices"]]
        self.gf = [np.array(x)-1 for x in d["gfl_indices"]]
        self.keep = np.array(d["keep_indices"])-1
        self.ig = d["gauge_index"]-1
        self.center = np.zeros(self.nx) if center is None else np.asarray(center,dtype=float)
        self.last_inverse_residual = None

    def voltage(self, coeff, source):
        if hasattr(source[0], "_mpi_"):
            g = mp.iv.matrix(self.Y.tolist())
            for i,k in enumerate(coeff):
                j=2*(i+29);g[j,j+1]-=k;g[j+1,j]+=k
            vv = mp.iv.lu_solve(g,mp.iv.matrix([-q for q in source]))
            return list(vv)
        interval = isinstance(source[0], Jet)
        if not interval:
            g = self.Y.copy()
            for i, k in enumerate(coeff):
                j = 2*(i+29); g[j,j+1] -= k; g[j+1,j] += k
            return list(-np.linalg.solve(g, np.array(source)))
        gl = self.Y.copy(); gh = self.Y.copy()
        for i, k in enumerate(coeff):
            j = 2*(i+29)
            gl[j,j+1], gh[j,j+1] = add((gl[j,j+1], gh[j,j+1]), neg(k.v))
            gl[j+1,j], gh[j+1,j] = add((gl[j+1,j], gh[j+1,j]), k.v)
        kin, check = inverse_enclosure((gl, gh)); self.last_inverse_residual = check
        hv = (np.array([x.v[0] for x in source])[:,None], np.array([x.v[1] for x in source])[:,None])
        # Centered residual form avoids multiplying uncertain inverse by a large trim source.
        vc = -np.linalg.solve(gl/2+gh/2, (hv[0]+hv[1]).ravel()/2)
        residual = add(hv, mm((gl, gh), (vc[:,None], vc[:,None])))
        vin = add((vc[:,None], vc[:,None]), neg(mm(kin, residual)))
        vals = [Jet(vin[0][j,0], vin[1][j,0]) for j in range(78)]
        dsource = [Jet.from_parts((0.,0.), x.d) for x in source]
        for i, k in enumerate(coeff):
            j = 2*(i+29); kd = Jet.from_parts((0.,0.), k.d)
            dsource[j] = dsource[j] - vals[j+1]*kd
            dsource[j+1] = dsource[j+1] + vals[j]*kd
        di = (np.array([x.d[0] for x in dsource]), np.array([x.d[1] for x in dsource]))
        dv = neg(mm(kin, di))
        return [Jet.from_parts(vals[j].v, (dv[0][j], dv[1][j])) for j in range(78)]

    def evaluate(self, z, w):
        interval = isinstance(z[0], Jet)
        C = Jet if interval else (mp.iv.mpf if hasattr(z[0], "_mpi_") else float)
        z = [a+C(b) for a,b in zip(z,self.center)]
        d = self.raw
        x = [C(v) for v in self.x0]
        for k, idx in enumerate(self.keep):
            x[idx] = x[idx] + z[k]
        source = [C(v) if not isinstance(v, Jet) else v for v in w[:78]]
        coeff = []; sgaux = []
        for i, ix in enumerate(self.sg):
            p = {k: (v if isinstance(v, bool) else C(v)) for k,v in d["sg_parameters"][i].items()}
            pq,pd,ed,eq,omega,delta = [x[k] for k in ix[-6:]]
            assert lower(omega) > 0
            gd1 = (p["xdpp"]-p["xls"])/(p["xdp"]-p["xls"])
            gq1 = (p["xqpp"]-p["xls"])/(p["xqp"]-p["xls"])
            cd = gd1*eq+(1-gd1)*pd; cq = -gq1*ed+(1-gq1)*pq
            sn = sf("sin", delta); cs = sf("cos", delta)
            scale = (1-C(d["rho"][i])) * p["rating_mva"]/p["system_base_mva"]
            coeff.append(scale/(omega*p["xdpp"]))
            j = 2*(i+29)
            source[j] = source[j] + scale/p["xdpp"]*(sn*cd+cs*cq)
            source[j+1] = source[j+1] + scale/p["xdpp"]*(-cs*cd+sn*cq)
            gf = self.gf[i]
            source[j] = source[j] + C(d["rho"][i])*x[gf[6]]
            source[j+1] = source[j+1] + C(d["rho"][i])*x[gf[5]]
            sgaux.append((p,pq,pd,ed,eq,omega,cd,cq,sn,cs,gd1,gq1))
        v = self.voltage(coeff, source)
        nu = C(d["sg_parameters"][-1]["omega_base"])*(x[self.ig-1]-C(d["sg_parameters"][-1]["omega_frame"]))
        dx = [C(0.) for _ in x]
        for i, ix in enumerate(self.sg):
            p,pq,pd,ed,eq,omega,cd,cq,sn,cs,gd1,gq1 = sgaux[i]
            ur,ui = v[2*(i+29):2*(i+30)]
            vd = sn*ur-cs*ui; vq = cs*ur+sn*ui
            id_ = cd/p["xdpp"]-vq/(omega*p["xdpp"])
            iq = cq/p["xqpp"]+vd/(omega*p["xqpp"])
            gd2 = (p["xdp"]-p["xdpp"])/(p["xdp"]-p["xls"])**2
            gq2 = (p["xqp"]-p["xqpp"])/(p["xqp"]-p["xls"])**2
            if p["controlled"]:
                xg1,xg2,vfb,vfout,vr,vm = [x[k] for k in ix[:6]]
                assert lower(xg1) > upper(p["gov_vmin"]) and upper(xg1) < lower(p["gov_vmax"])
                assert lower(vr) > upper(p["avr_vr_min"]) and upper(vr) < lower(p["avr_vr_max"])
                assert lower(vfout) > 0
                a = sf("sqrt", p["avr_se1"]*p["avr_e1"]/(p["avr_se2"]*p["avr_e2"]))
                thresh = p["avr_e2"]-(p["avr_e1"]-p["avr_e2"])/(a-1)
                bc = p["avr_se2"]*p["avr_e2"]*(a-1)**2/(p["avr_e1"]-p["avr_e2"])**2
                if upper(vfout) < lower(thresh):
                    ceil = C(0.)
                elif lower(vfout) > upper(thresh):
                    ceil = bc*(vfout-thresh)**2
                else:
                    raise ValueError("AVR saturation branch crossed by interval domain")
                vfoutdot = (vr-ceil-p["avr_ke"]*vfout)/p["avr_te"]
                vfbdot = (p["avr_kf"]*vfoutdot-vfb)/p["avr_tf"]
                vrdot = (p["avr_ka"]*(p["avr_vref"]-vm-vfb)-vr)/p["avr_ta"]
                vmdot = (sf("sqrt", ur**2+ui**2)-vm)/p["avr_tr"]
                xg1dot = ((p["gov_p_ref"]-(omega-p["gov_omega_ref"]))/p["gov_r"]-xg1)/p["gov_t1"]
                xg2dot = (xg1+p["gov_t2"]*xg1dot-xg2)/p["gov_t3"]
                taum = (xg2-p["gov_dt"]*(omega-p["gov_omega_ref"]))/omega
                vf = vfout
                for j,val in zip(ix[:6], (xg1dot,xg2dot,vfbdot,vfoutdot,vrdot,vmdot)):
                    dx[j] = val
            else:
                taum = p["tau_m_set"]; vf = p["vf_set"]
            taue = (cd*vd+cq*vq)/(omega*p["xdpp"])
            omdot = (taum-taue-p["damping"]*(omega-1))/(2*p["inertia"])
            eqdot = (-eq-(p["xd"]-p["xdp"])*(id_-gd2*pd-(1-gd1)*id_+gd2*eq)+vf)/p["tdp"]
            eddot = (-ed+(p["xq"]-p["xqp"])*(iq-gq2*pq-(1-gq1)*iq-gq2*ed))/p["tqp"]
            pddot = (-pd+eq-(p["xdp"]-p["xls"])*id_)/p["tdpp"]
            pqdot = (-pq-ed-(p["xqp"]-p["xls"])*iq)/p["tqpp"]
            deldot = p["omega_base"]*(omega-p["omega_frame"])-nu
            for j,val in zip(ix[-6:], (pqdot,pddot,eddot,eqdot,omdot,deldot)):
                dx[j] = val

            ixg = self.gf[i]
            gammaq,gammad,theta,om,xi,ifi,ifr,vdi,vdc = [x[k] for k in ixg]
            assert lower(vdc) > 0
            gp = {k: C(v) for k,v in d["gfl_parameters"][i].items()}
            c = sf("cos",theta); s = sf("sin",theta)
            idg = c*ifr+s*ifi; iqg = -s*ifr+c*ifi
            edg = (vdc-gp["Vdc"])*gp["dc_kp"]+vdi-idg; eqg = gp["iset_q"]-iqg
            vid = gp["cc_kp"]*edg+gp["cc_ki"]*gammad
            viq = gp["cc_kp"]*eqg+gp["cc_ki"]*gammaq
            e = -s*ur+c*ui
            fg = [eqg,edg,om-nu,(xi+C(d["Kp"][i])*e-om)/gp["pll_tau"],C(d["Ki"][i])*e,
                  gp["omega_base"]/gp["Xf"]*(s*vid+c*viq-ui-gp["Rf"]*ifi-gp["omega_frame"]*gp["Xf"]*ifr)-nu*ifr,
                  gp["omega_base"]/gp["Xf"]*(c*vid-s*viq-ur-gp["Rf"]*ifr+gp["omega_frame"]*gp["Xf"]*ifi)+nu*ifi,
                  (vdc-gp["Vdc"])*gp["dc_ki"],(gp["Pdc"]+w[78+i]-vid*idg-viq*iqg)/(gp["Cdc"]*vdc)]
            for j,val in zip(ixg,fg): dx[j] = val
        out = [dx[j] for j in self.keep]+[C(0.) for _ in range(78)]
        freq=[];rocof=[]
        for b in range(39):
            c = sf("cos",C(d["phase0"][b])); s = sf("sin",C(d["phase0"][b]))
            ur,ui = v[2*b:2*b+2]; den = c*ur+s*ui
            assert lower(den) > 0
            ph = sf("atan",(c*ui-s*ur)/den)
            eta = z[len(self.keep)+b]; chi = z[len(self.keep)+39+b]
            ff = (ph-eta)/(C(2*math.pi)*C(d["tau_f"]))
            rr = (ff-chi)/C(d["tau_r"])
            out[len(self.keep)+b] = (ph-eta)/C(d["tau_f"])-nu
            out[len(self.keep)+39+b] = rr
            freq.append(ff);rocof.append(rr)
        return out,v,freq,rocof

    def intervals(self, xr, wr):
        Jet.n = self.nx+self.nw
        zz=[];ww=[]
        for j,r in enumerate(np.r_[xr,wr]):
            grad = np.zeros(Jet.n);grad[j]=1.
            q = Jet(-r,r,grad,grad.copy())
            (zz if j<self.nx else ww).append(q)
        return self.evaluate(zz,ww)
