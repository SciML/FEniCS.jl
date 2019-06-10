using PyCall

if Sys.iswindows()
	@warn("""
					FEniCS is not known to run on Windows!
					Please try installing from Windows Subsystem for Linux.
		""")
end
	try
		pyimport_conda("fenics", "fenics=2019.1.0", "conda-forge")
		pyimport_conda("mshr", "mshr=2019.1.0", "conda-forge")
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
