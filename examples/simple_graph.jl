using Plots, GraphRecipes

include("../src/include_lib.jl")

# Define a simple graph and plot
n_elem = 17
n_pt = n_elem + 1
graph = StaticGraph(path_graph(n_pt))

plot(path_graph(n_pt), curves=false)




