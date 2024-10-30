function neighbors(graph, vtx)
	neighbors = findnz(graph[vtx, :])[1]
	return neighbors
end

function edges(graph)
    is, js, ids = findnz(graph)

    # Filter list
    pos_ids = ids .>= 0.0
    is_pos = is[pos_ids]
    js_pos = js[pos_ids]
    ids_pos = ids[pos_ids]

    # Order list
    order_dict = Dict(v => i for (i, v) in enumerate(ids_pos))

    ids_sort = sortperm(ids_pos, by=x -> order_dict[x])
    is_sort = is_pos[ids_sort]
    js_sort = js_pos[ids_sort]
	return zip(is_sort, js_sort)
end

function create_graph(edges, pts)
    n_pts = length(pts)
	graph = spzeros(Int, n_pts, n_pts)
	for (i, edge) in enumerate(edges)
		graph[edge[1], edge[2]] = i
		graph[edge[2], edge[1]] = -i
	end
	return graph
end

function edge_index(edge, graph)
	return abs(graph[edge[1], edge[2]])
end
