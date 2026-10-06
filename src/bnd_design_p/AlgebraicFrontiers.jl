module AlgebraicFrontiers

using LinearAlgebra, DataFrames
using ..ExpP
const N = ExpP.PDExactDesignN

export pll_pair_identity, single_sg_identity, two_sg_identity,
       pi_boundary_cases, retention_roots, retention_roots_all_buses

_rel(a, b) = norm(a - b) / max(norm(b), eps(Float64))

function _state_rows(m, bus)
    r = only(eachrow(m.state_map[(m.state_map.bus .== bus) .&
                                  (m.state_map.kind .== "GFL"), :]))
    first = Int(r.first)
    first + 3, first + 4
end

function _det_ratio(P0, P1)
    det(P0 \ P1)
end

function _det3(A)
    A[1,1] * (A[2,2]*A[3,3]-A[2,3]*A[3,2]) -
    A[1,2] * (A[2,1]*A[3,3]-A[2,3]*A[3,1]) +
    A[1,3] * (A[2,1]*A[3,2]-A[2,2]*A[3,1])
end

function symbolic_rank_one_check()
    R = Rational{BigInt}
    A = R[1 2 0; 0 1 1; 1 0 1]
    u = R[1,-1,2]; v = R[0,1,1]
    q = R[2,0,-1]
    a = R(3); b = R(-2)
    d0 = _det3(A)
    dp = _det3(A + a*u*transpose(q))
    di = _det3(A + b*v*transpose(q))
    dpi = _det3(A + (a*u+b*v)*transpose(q))
    same_cross = dpi - dp - di + d0

    q2 = R[1,1,0]
    d1 = _det3(A + a*u*transpose(q))
    d2 = _det3(A + b*v*transpose(q2))
    d12 = _det3(A + a*u*transpose(q) + b*v*transpose(q2))
    cross_devices = d12-d1-d2+d0
    (;same_device_cross=Float64(same_cross),
      different_device_cross=Float64(cross_devices),
      exact= iszero(same_cross) && !iszero(cross_devices))
end

"""Verify the pairwise PLL pencil update on all ten physical GFL rows."""
function pll_pair_identity(ctx, rho, kp0, ki0;
        frequencies=ComplexF64[-0.05 + 0.001im, 4.0 + 9.0im])
    m0 = N.descriptor(ctx, rho, kp0, ki0)
    rows = NamedTuple[]
    factors = Dict{Int,NamedTuple}()
    Δp = 0.37N.K0P
    Δi = 0.29N.K0I
    for bus in 30:39
        k = bus-29
        ip, ii = _state_rows(m0, bus)
        kp = copy(kp0); ki = copy(ki0)
        kp[k] += Δp; ki[k] += Δi
        kp_only = copy(kp0); kp_only[k] += Δp
        ki_only = copy(ki0); ki_only[k] += Δi
        mp = N.descriptor(ctx, rho, kp_only, ki0)
        mi = N.descriptor(ctx, rho, kp0, ki_only)
        mb = N.descriptor(ctx, rho, kp, ki)
        dAp = mp.Ared - m0.Ared
        dAi = mi.Ared - m0.Ared
        dAb = mb.Ared - m0.Ared
        tau = N.trim_gfl(ctx, bus).pars.pll_tau
        qP = vec(dAp[ip, :]) .* (tau / Δp)
        qI = vec(dAi[ii, :]) ./ Δi
        qerr = _rel(qP, qI)
        offP = norm(dAp[setdiff(1:size(dAp,1), [ip]), :]) / max(norm(dAp), eps())
        offI = norm(dAi[setdiff(1:size(dAi,1), [ii]), :]) / max(norm(dAi), eps())
        sv = svdvals(dAb)
        rank_resid = length(sv) >= 2 ? sv[2] / max(sv[1], eps()) : 0.0
        affine_resid = _rel(dAb, dAp+dAi)
        factors[bus] = (;q=qI, ip, ii, tau)
        u = zeros(Float64, size(m0.Ared,1))
        u[ip] = -Δp/tau
        u[ii] = -Δi
        for s in frequencies
            P0 = s*I(size(m0.Ared,1)) - m0.Ared
            Pp = s*I(size(mp.Ared,1)) - mp.Ared
            Pi = s*I(size(mi.Ared,1)) - mi.Ared
            Pb = s*I(size(mb.Ared,1)) - mb.Ared
            rp = _det_ratio(P0, Pp)
            ri = _det_ratio(P0, Pi)
            rb = _det_ratio(P0, Pb)
            lemma = 1 + transpose(qI) * (P0 \ u)
            push!(rows, (;bus, s_real=real(s), s_imag=imag(s),
                rank2_over_rank1=rank_resid, common_detector_rel_error=qerr,
                offrow_Kp=offP, offrow_Ki=offI, matrix_affine_residual=affine_resid,
                determinant_same_device_cross=abs(rb-rp-ri+1),
                determinant_lemma_abs_error=abs(rb-lemma[1]),
                reference_pencil_condition=cond(P0)))
        end
    end
    symbolic = symbolic_rank_one_check()
    # A two-device update is genuinely rank two and its determinant contains a
    # cross term for a coupled network; report the largest deterministic witness.
    cross_rows = NamedTuple[]
    for s in frequencies, i in 30:39, j in (i+1):39
        fi=factors[i];fj=factors[j]
        ui=zeros(Float64,size(m0.Ared,1));uj=similar(ui);fill!(uj,0)
        ui[fi.ip]=-Δp/fi.tau;uj[fj.ip]=-Δp/fj.tau
        U=hcat(ui,uj);V=hcat(fi.q,fj.q)
        P0=s*I(size(m0.Ared,1)) - m0.Ared
        ratio12=det(Matrix{ComplexF64}(I,2,2)+transpose(V)*(P0\U))
        ratio1=1+transpose(fi.q)*(P0\ui)
        ratio2=1+transpose(fj.q)*(P0\uj)
        push!(cross_rows,(;bus_i=i,bus_j=j,s_real=real(s),s_imag=imag(s),
            cross=abs(ratio12-ratio1[1]*ratio2[1])))
    end
    sort!(cross_rows; by=x->x.cross, rev=true)
    (;table=DataFrame(rows), factors, symbolic,
      max_between_device_cross=first(cross_rows), frequencies)
end

function _selector(n, bus)
    Π=zeros(Float64,n,2)
    Π[2bus-1,1]=1;Π[2bus,2]=1
    Π
end

function _gfl_port(ctx, bus, kp, ki)
    op=N.trim_gfl(ctx,bus)
    N.gfl_jacobians(op,kp[bus-29],ki[bus-29])
end

"""Check the exact rank-two single-retained-SG port determinant identity."""
function single_sg_identity(ctx, kp, ki;
        frequencies=ComplexF64[3+6im, -0.05+0.2im], epsilon=0.37)
    allgfl=ones(10);Trows=NamedTuple[]
    for s in frequencies
        TF=ExpP.nodal_operator(ctx,allgfl,kp,ki,s)
        for bus in 30:39
            Π=_selector(size(TF,1),bus)
            YF=N.port_transfer(_gfl_port(ctx,bus,kp,ki),s)
            YSG=N.port_transfer(ctx.net.sg[bus].J,s)
            # The installed component port is current into the device. In
            # network-injection convention Yinj=-Yraw, so the ε update is
            # Yraw_SG-Yraw_GFL.
            D=YSG-YF
            N2=D*transpose(Π)*(TF\Π)
            rho=ones(10);rho[bus-29]=1-epsilon
            Tmix=ExpP.nodal_operator(ctx,rho,kp,ki,s)
            Tpred=TF+epsilon*Π*D*transpose(Π)
            relop=_rel(Tmix,Tpred)
            fullratio=det(TF\Tmix)
            smallratio=det(Matrix{ComplexF64}(I,2,2)+epsilon*N2)
            push!(Trows,(;bus,s_real=real(s),s_imag=imag(s),epsilon,
                operator_relative_error=relop,determinant_identity_abs_error=abs(fullratio-smallratio),
                determinant_ratio=fullratio,quadratic_coeff_b=tr(N2),
                quadratic_coeff_c=det(N2),TF_condition=cond(TF)))
        end
    end
    DataFrame(Trows)
end

"""Check the 4×4 update, including cross-path terms for two retained SG ports."""
function two_sg_identity(ctx,kp,ki;
        frequencies=ComplexF64[3+6im,-0.05+0.2im],buses=(31,38),epsilons=(0.23,0.37))
    rows=NamedTuple[];rho=ones(10)
    rho[buses[1]-29]=1-epsilons[1];rho[buses[2]-29]=1-epsilons[2]
    for s in frequencies
        TF=ExpP.nodal_operator(ctx,ones(10),kp,ki,s)
        pis=[_selector(size(TF,1),b) for b in buses]
        Π=hcat(pis...)
        Ds=map(buses) do bus
            N.port_transfer(ctx.net.sg[bus].J,s)-N.port_transfer(_gfl_port(ctx,bus,kp,ki),s)
        end
        De=zeros(ComplexF64,4,4)
        De[1:2,1:2].=epsilons[1].*Ds[1]
        De[3:4,3:4].=epsilons[2].*Ds[2]
        N4=De*transpose(Π)*(TF\Π)
        Tmix=ExpP.nodal_operator(ctx,rho,kp,ki,s)
        Tpred=TF+Π*De*transpose(Π)
        det4=det(Matrix{ComplexF64}(I,4,4)+N4)
        ind1=det(Matrix{ComplexF64}(I,2,2)+epsilons[1].*Ds[1]*transpose(pis[1])*(TF\pis[1]))
        ind2=det(Matrix{ComplexF64}(I,2,2)+epsilons[2].*Ds[2]*transpose(pis[2])*(TF\pis[2]))
        push!(rows,(;s_real=real(s),s_imag=imag(s),bus_i=buses[1],bus_j=buses[2],
            operator_relative_error=_rel(Tmix,Tpred),
            determinant_identity_abs_error=abs(det(TF\Tmix)-det4),
            cross_term=abs(det4-ind1*ind2),determinant_4x4=det4,TF_condition=cond(TF)))
    end
    DataFrame(rows)
end

function _line_box_point(a, lo, hi)
    # Closest Euclidean point to the origin on aᵀx=-1 intersected with a box.
    pts=Vector{Vector{Float64}}()
    for j in 1:2
        k=3-j
        abs(a[j])>1e-14 || continue
        for edge in (lo[k],hi[k])
            x=zeros(2);x[k]=edge;x[j]=(-1-a[k]*edge)/a[j]
            lo[j]-1e-10<=x[j]<=hi[j]+1e-10 && push!(pts,x)
        end
    end
    if norm(a)>0
        x=-a/(dot(a,a))
        all(lo.-1e-10 .<= x .<= hi.+1e-10) && push!(pts,x)
    end
    isempty(pts) ? nothing : pts[argmin(norm.(pts))]
end

"""Compute selected conditioned PI boundary points; ω=0 remains a real line."""
function pi_boundary_cases(ctx,rho,kp_base,ki_base,bus;
        frequencies=[0.0,0.1,0.5,2.0],sigma=0.05)
    i=bus-29
    m0=N.descriptor(ctx,rho,kp_base,ki_base)
    ip,ii=_state_rows(m0,bus)
    Δp=(kp_base[i]+N.K0P<=ctx.kpmax[i]) ? N.K0P : -N.K0P
    Δi=(ki_base[i]+N.K0I<=ctx.kimax[i]) ? N.K0I : -N.K0I
    kp=copy(kp_base);kp[i]+=Δp
    ki=copy(ki_base);ki[i]+=Δi
    mp=N.descriptor(ctx,rho,kp,ki_base)
    mi=N.descriptor(ctx,rho,kp_base,ki)
    rowp=vec((mp.Ared-m0.Ared)[ip,:])./Δp
    rowi=vec((mi.Ared-m0.Ared)[ii,:])./Δi
    tau=N.trim_gfl(ctx,bus).pars.pll_tau
    qP=tau.*rowp;qI=rowi
    qerr=_rel(qP,qI)
    up=zeros(Float64,size(m0.Ared,1));up[ip]=1/tau
    ui=zeros(Float64,size(m0.Ared,1));ui[ii]=1
    lo=[ctx.kpmin[i]-kp_base[i],ctx.kimin[i]-ki_base[i]] ./ [N.K0P,N.K0I]
    hi=[ctx.kpmax[i]-kp_base[i],ctx.kimax[i]-ki_base[i]] ./ [N.K0P,N.K0I]
    rows=NamedTuple[]
    for omega in frequencies
        s=complex(-sigma,omega)
        P=s*I(size(m0.Ared,1))-m0.Ared
        cp=-sum(qI .* (P\up))*N.K0P
        ci=-sum(qI .* (P\ui))*N.K0I
        M=[real(cp) real(ci);imag(cp) imag(ci)]
        if omega==0
            a=[real(cp),real(ci)]
            x=_line_box_point(a,lo,hi)
            status=x===nothing ? "REAL_LINE_NO_BOX_INTERSECTION" : "REAL_LINE"
            condM=Inf;resid=x===nothing ? NaN : abs(1+cp*x[1]+ci*x[2])
        else
            condM=cond(M)
            if !isfinite(condM) || condM>1e12
                x=nothing;status="SINGULAR_OR_ILL_CONDITIONED_2X2";resid=NaN
            else
                x=M\[-1.0,0.0]
                inbox=all(lo.-1e-9 .<= x .<= hi.+1e-9)
                resid=abs(1+cp*x[1]+ci*x[2])
                status=inbox ? "POINT_IN_GAIN_BOX" : "POINT_OUTSIDE_GAIN_BOX"
            end
        end
        λerr=NaN;alpha=NaN;right_count=-1;boundary_pass=false
        cross_lo_alpha=NaN;cross_hi_alpha=NaN;cross_lo_count=-1;cross_hi_count=-1
        cross_lo_in_box=false;cross_hi_in_box=false
        in_box = x!==nothing && all(lo.-1e-9 .<= x .<= hi.+1e-9)
        if in_box && all(isfinite,x)
            kt=copy(kp_base);it=copy(ki_base)
            kt[i]+=x[1]*N.K0P;it[i]+=x[2]*N.K0I
            st=N.spectrum(ctx,rho,kt,it)
            λerr=minimum(abs.(st.lambda .- s));alpha=st.alpha
            right_count=count(real.(st.lambda).>=-sigma)
            boundary_pass=λerr<1e-5
            normal=[real(cp),real(ci)]
            norm(normal)>0 && (normal./=norm(normal))
            step=1e-5
            for side in (-1,1)
                xs=x+side*step*normal
                ks=copy(kp_base);is=copy(ki_base)
                ks[i]+=xs[1]*N.K0P;is[i]+=xs[2]*N.K0I
                within=all(lo.-1e-10 .<= xs .<= hi.+1e-10)
                # For geometric falsification only, permit exactly this one
                # off-box point by widening that coordinate just enough. It is
                # explicitly not a design-feasible gain vector.
                cside=ctx
                if !within
                    kpmin=copy(ctx.kpmin);kpmax=copy(ctx.kpmax)
                    kimin=copy(ctx.kimin);kimax=copy(ctx.kimax)
                    kpmin[i]=min(kpmin[i],ks[i]);kpmax[i]=max(kpmax[i],ks[i])
                    kimin[i]=min(kimin[i],is[i]);kimax[i]=max(kimax[i],is[i])
                    cside=N.DesignContext(ctx.root,ctx.net,ctx.original,ctx.power,
                        kpmin,kpmax,kimin,kimax)
                end
                ss=N.spectrum(cside,rho,ks,is)
                if side<0
                    cross_lo_alpha=ss.alpha;cross_lo_count=count(real.(ss.lambda).>=-sigma)
                    cross_lo_in_box=within
                else
                    cross_hi_alpha=ss.alpha;cross_hi_count=count(real.(ss.lambda).>=-sigma)
                    cross_hi_in_box=within
                end
            end
        end
        push!(rows,(;bus,omega,status,condition_2x2=condM,detector_row_rel_error=qerr,
            delta_Kp_normalized=x===nothing ? NaN : x[1],
            delta_Ki_normalized=x===nothing ? NaN : x[2],
            determinant_residual=resid,boundary_pole_error=λerr,
            alpha=alpha,physical_poles_at_or_right_of_margin=right_count,
            boundary_pole_verified=boundary_pass,cross_minus_alpha=cross_lo_alpha,
            cross_plus_alpha=cross_hi_alpha,cross_minus_pole_count=cross_lo_count,
            cross_plus_pole_count=cross_hi_count,cross_minus_in_gain_box=cross_lo_in_box,
            cross_plus_in_gain_box=cross_hi_in_box))
    end
    DataFrame(rows)
end

"""All physical retention roots from eig(N), verified with the complete spectrum."""
function retention_roots(ctx,kp,ki,bus;sigma=0.05,delta_num=1e-9)
    s=-sigma
    TF=ExpP.nodal_operator(ctx,ones(10),kp,ki,s)
    Π=_selector(size(TF,1),bus)
    D=N.port_transfer(ctx.net.sg[bus].J,s)-N.port_transfer(_gfl_port(ctx,bus,kp,ki),s)
    N2=D*transpose(Π)*(TF\Π)
    evals=eigvals(Matrix(N2))
    roots=ComplexF64[-1/x for x in evals if abs(x)>eps(Float64)]
    records=NamedTuple[]
    for (ri,root) in enumerate(roots)
        abs(imag(root))<1e-7 || continue
        epsv=real(root)
        (-1e-9<=epsv<=1+1e-9) || continue
        epsv=clamp(epsv,0,1)
        rho=ones(10);rho[bus-29]=1-epsv
        sp=N.spectrum(ctx,rho,kp,ki)
        near=findall(abs.(sp.lambda .- s).<1e-5)
        target_pass=sp.alpha<=-sigma+1e-7 && !isempty(near)
        b=tr(N2);c=det(N2)
        polynomial=abs(1+b*epsv+c*epsv^2)
        push!(records,(;bus,root_index=ri,epsilon=epsv,rho=1-epsv,
            retained_SG_MW=ctx.power[bus-29]*epsv,
            converted_GFL_MW=sum(ctx.power)-ctx.power[bus-29]*epsv,
            quadratic_b=real(b),quadratic_c=real(c),polynomial_residual=polynomial,
            real_eigenvalue=minimum(real.(sp.lambda)),alpha=sp.alpha,
            boundary_pole_count=length(near),all_physical_poles_pass=target_pass,
            valid_architecture=0<epsv<=1,TF_condition=cond(TF)))
    end
    # Separately include the prescribed numerical guard at sigma+delta_num.
    sg=-sigma-delta_num
    TFg=ExpP.nodal_operator(ctx,ones(10),kp,ki,sg)
    Dg=N.port_transfer(ctx.net.sg[bus].J,sg)-N.port_transfer(_gfl_port(ctx,bus,kp,ki),sg)
    Ng=Dg*transpose(Π)*(TFg\Π)
    guarded=ComplexF64[-1/x for x in eigvals(Matrix(Ng)) if abs(x)>eps(Float64)]
    grecords=NamedTuple[]
    for (ri,r) in enumerate(guarded)
        abs(imag(r))<1e-7 || continue
        e=real(r);0<e<=1 || continue
        rho=ones(10);rho[bus-29]=1-e
        sp=N.spectrum(ctx,rho,kp,ki)
        push!(grecords,(;bus,root_index=ri,epsilon=e,
            retained_SG_MW=ctx.power[bus-29]*e,alpha=sp.alpha,
            target=-sigma-delta_num,
            guarded_boundary_pass=sp.alpha<=-sigma-delta_num+1e-7,
            boundary_pole_count=count(abs.(sp.lambda.-sg).<1e-5)))
    end
    (;roots=DataFrame(records),guarded=DataFrame(grecords),N2,b=tr(N2),c=det(N2))
end

function retention_roots_all_buses(ctx,kp,ki)
    roots=DataFrame[];guarded=DataFrame[]
    for bus in 30:39
        r=retention_roots(ctx,kp,ki,bus)
        !isempty(r.roots) && push!(roots,r.roots)
        !isempty(r.guarded) && push!(guarded,r.guarded)
    end
    (;roots=isempty(roots) ? DataFrame() : vcat(roots...),
      guarded=isempty(guarded) ? DataFrame() : vcat(guarded...))
end

end
