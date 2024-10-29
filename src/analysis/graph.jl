function neighbors(graph, vtx)
	neighbors = findnz(graph[vtx, :])[2]
	return neighbors
end

function edges(graph)
	return findnz(graph)
end

function create_graph(edges, pts)
    n_pts = length(pts)
	graph = spzeros(Int, n_pts, n_pts)
	for edge in edges
		graph[edge[1], edge[2]] = 1
		graph[edge[2], edge[1]] = 1
	end
	return graph
end

function edge_index(adj_matrix, i, j)
	row_indices, col_indices, _ = findnz(adj_matrix)
	for idx in 1:length(row_indices)
		if row_indices[idx] == i && col_indices[idx] == j
			return idx  # Position of the edge in the sparse storage
		end
	end
	return nothing  # Edge not found
end
