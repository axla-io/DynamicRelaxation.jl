using Graphs
using MetaGraphsNext
using StaticArrays
using LinearAlgebra
using NBodySimulator
using StaticGraphs

using Plots, GraphRecipes

include("src/include_lib.jl")

# Define a simple graph and plot
n_elem = 16
n_pt = n_elem + 1
graph = StaticGraph(path_graph(n_pt))

#plot(graph, curves = false)

E = 210 * 1e9               # [Pa]
Iy = Iz = 4.7619 * 1e-7     # [m^4]
A = 4.7619 * 1e-4           # [m^2]
G = 78 * 1e9                # [Pa]
It = 2 * Iy                 # [m^4]

ep = ElementProperties{Float64}(E, A, Iy, Iz, G, It)

M = SMatrix{3, 3, Float64}(I)
M_inv = SMatrix{3, 3, Float64}(I)
J = SMatrix{3, 3, Float64}(I)
J_inv  = SMatrix{3, 3, Float64}(I)

np1 = NodeProperties{Float64}(M, M_inv, J, J_inv , false, @SVector(zeros(Bool, 7)));

np_fix = NodeProperties{Float64}(M, M_inv, J, J_inv , true, @SVector(zeros(Bool, 7)));




nodes_description = vcat(1 => np_fix, [i => np1 for i in 2:n_pt]...) #fill with properties

edges_description = [(src(edge), dst(edge)) => ep for edge in edges(graph)] #fill with properties

structure = MetaGraph(graph, nodes_description, edges_description, "simple structure")
typeof(graph)