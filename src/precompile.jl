@setup_workload begin
    @compile_workload begin
        if !PyCall.ispynull(fenics)
            mesh = UnitSquareMesh(2, 2)
            space = FunctionSpace(mesh, "P", 1)
            expression = Expression("x[0] + x[1]", degree = 1)
            interpolate(expression, space)
        end
    end
end
