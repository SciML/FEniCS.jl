"""
Creates a type which incorporates the mesh attributes alongside the necessary values
to solve finite element problems. It works by taking in a FEniCS mesh, working out
the necessary values, and then creating an ordering of the internal/boundary nodes
and applying the DirichletBC (via a function) to the respective nodes.
"""
mutable struct feMesh
    n_nodes::Int
    n_elements::Int
    nodes::Array{Float64, 2}
    internal_nodes::Vector
    boundary_nodes::Vector
    elements::Array{Int, 2}
    n_internal_nodes::Int
    n_boundary_nodes::Int
    node_vals::Vector

    function feMesh(mesh::Mesh,boundaryCondition::Base.Function)
        n_nodes = num_vertices(mesh)
        n_elements = num_cells(mesh)
        nodes = coordinates(mesh)
        bmesh = BoundaryMesh(mesh,"exterior")
        external_nodes = coordinates(bmesh)
        n_bc_nodes = size(external_nodes)[1]
        n_internal_nodes = n_nodes - n_bc_nodes
        node_vals = zeros(n_nodes)
        #this applied the boundary values to node_vals, and also calculates the numbering for internal nodes
        boundary_nodes, internal_nodes = find_node_number(n_nodes, n_bc_nodes, external_nodes, nodes,boundaryCondition,node_vals)
        elements_temp = cells(mesh)
        #the following commands convert the element numbering to Int64 and handle the 0/1 difference between FEniCS(Python) and Julia
        elements_temp2 = convert.(Int,elements_temp)
        elements = elements_temp2 .+1
        new(n_nodes, n_elements, nodes, internal_nodes, boundary_nodes, elements, n_internal_nodes, n_bc_nodes, node_vals)
    end
end

export feMesh

"finds node numbering for a specific FeMesh. Not currently exported as it is only
used to create the FeMesh type."
function find_node_number(n_nodes, n_bc_nodes, external_nodes, nodes, boundaryCondition::Core.Function, node_vals)
    j=1
    k=1
    count_int=0
    boundary_nodes = zeros(Int64, n_nodes)
    internal_nodes = zeros(Int64, n_nodes)

    node_ext = [external_nodes[1,:]]
    ##maybe it's more efficient to create and then pass values instead of push?
    for i=2:n_bc_nodes
        push!(node_ext,external_nodes[i,:])
    end
    node_overall = [nodes[1,:]]
    for i =2:n_nodes
        push!(node_overall,nodes[i,:])
    end

    node_number = findall((in)(node_ext),node_overall)
    for i = 1:n_nodes
        if i in node_number
            @inbounds node_vals[i] = boundaryCondition(nodes[i,1],nodes[i,2])
        else
            count_int+=1
            @inbounds internal_nodes[i]=count_int
        end
        i+=1
    end
    return boundary_nodes, internal_nodes
end
