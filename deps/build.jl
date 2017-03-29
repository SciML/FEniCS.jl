using PyCall
using Conda

if PyCall.conda
	Conda.add_channel("conda-forge")
	Conda.add("fenics")
else
	try
		pyimport("fenics")
	catch ee
		typeof(ee) <: PyCall.PyError || rethrow(ee)
		warn("""
						Python Dependancies not installed
						Please either:
						 - Rebuild PyCall to use Conda, by running in the julia REPL:
						 - `ENV["PYTHON"]=""; Pkg.build("PyCall"); Pkg.build("FEniCS")`
						 - Or install the depencences separately
				 """
		)
	end

end
