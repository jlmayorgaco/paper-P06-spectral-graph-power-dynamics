module BNDGraph

include("GraphBasis.jl")
include("DynamicSelfEnergy.jl")
include("GraphModalOperator.jl")
include("IntermodalSelfEnergy.jl")
include("PolePredictor.jl")
include("CouplingMetrics.jl")
include("SyntheticSystems.jl")
include("ExpAAdapter.jl")
include("Reporting.jl")
include("GeneralizedDynamicSelfEnergy.jl")
include("GraphBackbone.jl")
include("IntermodalPathways.jl")
include("RealPoleAnalysis.jl")
include("RealSystemAdapter.jl")

using .GraphBasis
using .DynamicSelfEnergy
using .GraphModalOperator
using .IntermodalSelfEnergy
using .PolePredictor
using .CouplingMetrics
using .SyntheticSystems
using .ExpAAdapter
using .Reporting
using .GeneralizedDynamicSelfEnergy
using .GraphBackbone
using .IntermodalPathways
using .RealPoleAnalysis
using .RealSystemAdapter

export GraphBasis, DynamicSelfEnergy, GraphModalOperator, IntermodalSelfEnergy,
       PolePredictor, CouplingMetrics, SyntheticSystems, ExpAAdapter, Reporting
export GeneralizedDynamicSelfEnergy, GraphBackbone, IntermodalPathways,
       RealPoleAnalysis, RealSystemAdapter

end
