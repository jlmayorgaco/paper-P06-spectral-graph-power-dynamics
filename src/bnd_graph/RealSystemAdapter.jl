module RealSystemAdapter

using LinearAlgebra
using CSV
using DataFrames

export load_expA_primary, case_from_reduced, build_ieee39_cases,
       descriptor_reduction, named_coordinate_partition, angle_operator,
       write_matrix_csv, matrix_from_csv, model_versions

function matrix_from_csv(path)
    df = CSV.read(path, DataFrame)
    cols = filter(!=("row_index"), String.(names(df)))
    return Matrix{Float64}(Matrix(df[:, cols]))
end

function write_matrix_csv(path, A::AbstractMatrix)
    names = ["x$i" for i in axes(A, 2)]
    df = DataFrame(A, names)
    insertcols!(df, 1, :row_index => collect(1:size(A,1)))
    CSV.write(path, df)
    return path
end

function named_coordinate_partition(A::AbstractMatrix, names::Vector{String};
                                    bus_filter=nothing)
    length(names) == size(A,1) || throw(ArgumentError("state-name count must match Ared"))
    parsed = NamedTuple[]
    for (i,name) in enumerate(names)
        m = match(r"VIndex\((\d+),\s*:(.*)\)", name)
        m === nothing && throw(ArgumentError("cannot parse PowerDynamics state name: $name"))
        bus = parse(Int,m.captures[1])
        var = m.captures[2]
        if occursin("machine₊δ",var) || occursin("pll₊θ",var)
            kind = "q"
            device = occursin("pll₊θ",var) ? "gfl" : "machine"
        elseif occursin("machine₊ω",var) || occursin("pll₊Δω_rad_s",var)
            kind = "v"
            device = occursin("pll₊Δω_rad_s",var) ? "gfl" : "machine"
        else
            continue
        end
        push!(parsed,(index=i,bus=bus,device=device,kind=kind,name=name))
    end
    if bus_filter !== nothing
        parsed = filter(r -> r.bus in bus_filter, parsed)
    end
    qs = sort(filter(r -> r.kind=="q",parsed); by=r -> (r.bus,r.device))
    vs = sort(filter(r -> r.kind=="v",parsed); by=r -> (r.bus,r.device))
    [(r.bus,r.device) for r in qs] == [(r.bus,r.device) for r in vs] ||
        throw(ArgumentError("angle and speed state maps do not pair"))
    q = getproperty.(qs,:index)
    v = getproperty.(vs,:index)
    C = Matrix(A[q,v])
    offdiag = C - Diagonal(diag(C))
    norm(offdiag) <= 1e-9*max(norm(C),1.0) ||
        throw(ArgumentError("measured qdot=Ckin*v map is not diagonal"))
    all(diag(C) .> 0) || throw(ArgumentError("kinematic diagonal is not positive"))
    residual = norm(A[q,:] - sparse_row_map(A,q,v,C)) /
        max(norm(A[q,:]),eps(Float64))
    residual <= 1e-12 || throw(ArgumentError("qdot equations include additional linear terms"))
    retained = vcat(q,v)
    c = setdiff(collect(1:size(A,1)),retained)
    return (q=q,v=v,condensed=c,retained=retained,Ckin=C,
            pairs=qs,kinematic_residual=residual)
end

function sparse_row_map(A,q,v,C)
    B = zeros(Float64,length(q),size(A,2))
    B[:,v] = C
    return B
end

function descriptor_reduction(sys)
    n=size(sys.A,1)
    E=sys.M isa UniformScaling ? Matrix{Float64}(sys.M,n,n) : Matrix{Float64}(sys.M)
    A=Matrix{Float64}(sys.A)
    d=findall(==(1.0),diag(E))
    a=findall(==(0.0),diag(E))
    length(d)+length(a)==n || throw(ArgumentError("descriptor mask must be binary diagonal"))
    Fx=A[d,d]; Fy=A[d,a]; Gx=A[a,d]; Gy=A[a,a]
    if isempty(a)
        Ared=Fx
        gycond=1.0
        gysmin=Inf
    else
        sv=svdvals(Gy)
        gysmin=minimum(sv)
        gycond=maximum(sv)/max(gysmin,eps(Float64))
        gysmin>eps(Float64)*max(opnorm(Gy),1.0) ||
            throw(ArgumentError("algebraic descriptor block is singular"))
        Ared=Fx-Fy*(Gy\Gx)
    end
    return (E=E,A=A,Ared=Ared,d=d,a=a,Fx=Fx,Fy=Fy,Gx=Gx,Gy=Gy,
            Gy_condition=gycond,Gy_sigma_min=gysmin)
end

function angle_operator(A::AbstractMatrix, part)
    q,v,c=part.q,part.v,part.condensed
    C=part.Ckin
    M=C\Matrix{Float64}(I,size(C,1),size(C,2))
    Avv=A[v,v]; Avq=A[v,q]; Avc=A[v,c]
    Acq=A[c,q]; Acv=A[c,v]; Acc=A[c,c]
    D0=-Avv*M
    L0=-Avq
    return (M=Matrix(M),D0=Matrix(D0),L0=Matrix(L0),Avc=Matrix(Avc),
            Acc=Matrix(Acc),Acq=Matrix(Acq),Acv=Matrix(Acv))
end

function load_expA_primary(root)
    out=joinpath(root,"reports","experiment_A")
    A=matrix_from_csv(joinpath(out,"matrices","bus33_A_reduced.csv"))
    df=CSV.read(joinpath(out,"tables","TABLE_A03_state_partition.csv"),DataFrame)
    dynamic=filter(r -> r.differential_or_algebraic=="differential",df)
    n=size(A,1)
    names=fill("",n)
    for r in eachrow(dynamic)
        i=Int(r.reduced_state_index)
        names[i]=String(r.state_name)
    end
    all(!isempty,names) || throw(ArgumentError("ExpA state map does not cover Ared"))
    part=named_coordinate_partition(A,names)
    mats=angle_operator(A,part)
    return (Ared=A,state_names=names,partition=part,matrices=mats,
            source=joinpath(out,"matrices","bus33_A_reduced.csv"))
end

function _description(sym)
    text=string(sym)
    m=match(r"^VIndex\((\d+),\s*:(.*)\)$",text)
    m===nothing && return (bus=missing,var=text)
    return (bus=parse(Int,m.captures[1]),var=m.captures[2])
end

function _case_from_sys(label,bus,sys)
    blk=descriptor_reduction(sys)
    allnames=String.(string.(sys.sym))
    names=String[]
    for ix in blk.d
        push!(names,allnames[ix])
    end
    part=named_coordinate_partition(blk.Ared,names)
    mats=angle_operator(blk.Ared,part)
    return (label=label,bus=bus,sys=sys,blocks=blk,Ared=blk.Ared,
            state_names=names,partition=part,matrices=mats,
            equilibrium_state=nothing,equilibrium_residual=NaN,
            source="PowerDynamics.jl 5.0.0 IEEE-39 official example, exact A pipeline")
end

"""Rebuild baseline and frozen nominal single-GFL cases without changing ExpA."""
function build_ieee39_cases(PD39, PowerDynamics, NetworkDynamics; on_case=(id,case)->nothing)
    cases=Dict{String,Any}()
    base_nw=PD39.baseline_network()
    base_pf=PowerDynamics.solve_powerflow(base_nw;verbose=false,sparse=false)
    base_s=PowerDynamics.initialize_from_pf!(base_nw;pfs=base_pf,verbose=false,sparsepf=false)
    base_sys=NetworkDynamics.linearize_network(base_s)
    cases["C0"]=_case_from_sys("C0_all_SG",0,base_sys)
    cases["C0"]=merge(cases["C0"],(equilibrium_state=base_s,))
    on_case("C0",cases["C0"])
    template=PD39.simple_gfldc_template()
    for bus in (30,33,35,37)
        nw=PD39.replace_bus(base_nw,bus;template=template)
        pf=PowerDynamics.solve_powerflow(nw;verbose=false,sparse=false)
        s0=PowerDynamics.initialize_from_pf!(nw;pfs=pf,verbose=false,sparsepf=false)
        sys=NetworkDynamics.linearize_network(s0)
        c=_case_from_sys("C$(bus)_bus$(bus)_SimpleGFLDC",bus,sys)
        cases["C$(bus)"]=merge(c,(equilibrium_state=s0,))
        on_case("C$(bus)",cases["C$(bus)"])
    end
    return cases
end

function case_from_reduced(label,bus,Ared,state_names)
    part=named_coordinate_partition(Ared,state_names)
    mats=angle_operator(Ared,part)
    return (label=label,bus=bus,sys=nothing,blocks=nothing,Ared=Ared,
            state_names=state_names,partition=part,matrices=mats,
            equilibrium_state=nothing,equilibrium_residual=NaN,source="Experiment-A frozen export")
end

function model_versions(PowerDynamics,NetworkDynamics)
    return (powerdynamics=string(Base.pkgversion(PowerDynamics)),
            networkdynamics=string(Base.pkgversion(NetworkDynamics)))
end

end
