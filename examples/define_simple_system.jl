using Graphs
using MetaGraphsNext
using StaticArrays

using Plots, GraphRecipes

include("src/include_lib.jl")

# Define a simple graph and plot
n_elem = 16
n_pt = n_elem + 1
g = path_graph(n_pt)

plot(g, curves = false)

nodes_description = [] #fill with properties
edges_description = [] #fill with properties

colors2 = MetaGraph(graph, nodes_description, edges_description, "simple structure")