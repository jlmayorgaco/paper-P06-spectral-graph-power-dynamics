using LinearAlgebra
using TOML
BLAS.set_num_threads(1)

include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","IdealDamping.jl"))
using .IdealDamping

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_F1")
mkpath(joinpath(OUT,"tables")); mkpath(joinpath(OUT,"matrices"))

function deterministic_spd(n)
    A=[sin(0.37i+0.19j)+cos(0.23i-0.41j) for i in 1:n,j in 1:n]
    return A'*A+2I
end
function damp_matrix(Mhalf,U,nu,factor)
    Mhalf*U*Diagonal(factor .* (2sqrt.(nu)))*U'*Mhalf
end
function state_matrix(M,K,D)
    n=size(M,1)
    [zeros(n,n) Matrix{Float64}(I,n,n); -(M\K) -(M\D)]
end
function write_table(path,headers,rows)
    open(path,"w") do io
        println(io,join(headers,','))
        for row in rows
            vals=row isa NamedTuple ? collect(values(row)) : row
            println(io,join(string.(vals),','))
        end
    end
end
function write_matrix(path,A)
    headers=["row_index";["x$j" for j in 1:size(A,2)]]
    rows=[[i;A[i,:]] for i in 1:size(A,1)]
    write_table(path,headers,rows)
end

n=5
M=deterministic_spd(n)
ME=eigen(Symmetric(M)); Mhalf=ME.vectors*Diagonal(sqrt.(ME.values))*ME.vectors'
Q=Matrix(qr([sin(0.11i*j)+cos(0.29i+0.17j) for i in 1:n,j in 1:n]).Q)
nu=[0.6,1.5,4.0,9.0,16.0]
K=Mhalf*Q*Diagonal(nu)*Q'*Mhalf
design=modal_critical_damping(M,K)
Dcrit=design.D
factors=(underdamped=0.7,modal_critical=1.0,overdamped=1.8)
rows=NamedTuple[]
max_pole_error=0.0
max_modal_residual=0.0
for (scenario,factor) in pairs(factors)
    D=damp_matrix(design.Mhalf,design.U,design.nu,factor)
    Dtilde=design.Minvhalf*D*design.Minvhalf
    global max_modal_residual=max(max_modal_residual,norm(Dtilde-factor*design.U*Diagonal(2sqrt.(design.nu))*design.U',Inf))
    numeric=eigvals(state_matrix(M,K,D))
    expected=ComplexF64[]
    for mode in eachindex(nu)
        d=factor*2sqrt(nu[mode])
        roots=modal_poles(nu[mode],d)
        decay=isolated_decay_rate(nu[mode],d)
        append!(expected,roots)
        err=minimum(abs.(numeric .- roots[1]))
        err2=minimum(abs.(numeric .- roots[2]))
        global max_pole_error=max(max_pole_error,err,err2)
        push!(rows,(scenario=String(scenario),mode=mode,nu=nu[mode],damping=d,
            lambda_1_real=real(roots[1]),lambda_1_imag=imag(roots[1]),
            lambda_2_real=real(roots[2]),lambda_2_imag=imag(roots[2]),
            analytic_decay_rate=decay,numerical_spectral_abscissa=maximum(real.(numeric)),
            root_matching_error=max(err,err2)))
    end
    length(expected)==length(numeric) || error("state dimension and modal root count differ")
end
write_table(joinpath(OUT,"tables","TABLE_F01_modal_critical_damping.csv"),
    collect(string.(propertynames(first(rows)))),rows)
write_matrix(joinpath(OUT,"matrices","M_synthetic.csv"),M)
write_matrix(joinpath(OUT,"matrices","K_synthetic.csv"),K)
write_matrix(joinpath(OUT,"matrices","D_modal_critical.csv"),Dcrit)
pass=max_pole_error<=1e-6 && max_modal_residual<=1e-10
result=Dict("experiment"=>"BND_EXP_F1","status"=> (pass ? "PASS" : "FAIL"),
    "synthetic_dimension"=>n,"nu"=>nu,"damping_factors"=>Dict(String(k)=>v for (k,v) in pairs(factors)),
    "max_analytic_vs_numeric_pole_error"=>max_pole_error,
    "max_modal_damping_transform_residual"=>max_modal_residual,
    "D_modal_formula"=>"2 M^(1/2) (M^(-1/2) K M^(-1/2))^(1/2) M^(1/2)",
    "claim_scope"=>"Each isolated commuting mode is optimally damped; no global optimum over arbitrary noncommuting D is claimed.",
    "nonproportional_reference"=>"not implemented; no verified Weyl-Horn construction was available in the local audited theory")
open(joinpath(OUT,"RESULTS_EXP_F1.toml"),"w") do io TOML.print(io,result) end
open(joinpath(OUT,"REPORT_EXP_F1.md"),"w") do io
    println(io,"# Experiment F1 — Ideal second-order graph damping\n")
    println(io,"**F1_STATUS: ",pass ? "PASS" : "FAIL","**\n")
    println(io,"For `M>0`, `K>=0`, mass-normalization gives `K̃=M^{-1/2}KM^{-1/2}`. In the commuting modal basis, each scalar equation is `s²+dₖs+νₖ=0`. For underdamping, the decay rate is `dₖ/2`; after critical damping, it is `(dₖ-√(dₖ²-4νₖ))/2`, which decreases as `dₖ` increases. Therefore the isolated-mode maximum occurs at `dₖ*=2√νₖ`, yielding `D*=2M^{1/2}K̃^{1/2}M^{1/2}`. This is not a claim of global optimality over arbitrary noncommuting damping matrices.\n")
    println(io,"A deterministic synthetic 5-mode SPD case compared underdamped, critical, and overdamped modal damping. Maximum analytic-to-numeric pole error: `$(max_pole_error)`; modal transformation residual: `$(max_modal_residual)`.\n")
    println(io,"No nonproportional optimum or Weyl–Horn matrix is proposed: this experiment limits itself to the proved modal result and explicitly marks the broader reference unimplemented. No GFL claim is made.\n")
end
println("EXPERIMENT_F1_STATUS: ",pass ? "PASS" : "FAIL")
println("MAX_ANALYTIC_NUMERIC_POLE_ERROR: ",max_pole_error)
println("MODAL_TRANSFORM_RESIDUAL: ",max_modal_residual)
