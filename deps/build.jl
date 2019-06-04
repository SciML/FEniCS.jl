using PyCall

	try
		pyimport_conda("fenics", "fenics", "conda-forge")
		pyimport_conda("mshr", "mshr", "conda-forge")
	catch ee
		typeof(ee) <: PyCall.PyError || rethrow(ee)
		@warn("""
						Python Dependancies not installed
						Please either:
						 - Rebuild PyCall using the path to FEniCS using
						 - `ENV["PYTHON"]="/path/to/FEniCS"; Pkg.build("PyCall"); Pkg.build("FEniCS")`
						 - Or install the Dependancies using the Docker image provided
				 """
		)
	end
