#!/bin/bash

JULIAVER=$1                     # the first and only argument to the script is the version

PKGNAME="FEniCS"


cd /mnt && if [[ -a .git/shallow ]]; then git fetch --unshallow; fi

# run tests
#julia -e "Pkg.clone(\"https://github.com/JuliaPy/PyCall.jl\"); ENV[\"PYTHON\"] = \"/usr/bin/python3\"; Pkg.build(\"PyCall\")"
julia -e "import Pkg; Pkg.clone(\"/mnt/\", \"$PKGNAME\"); Pkg.build(\"$PKGNAME\"); Pkg.test(\"$PKGNAME\"; coverage=true)"
TEST_EXIT=$?                    # return with this

# save coverage results back to host
PKGDIR=`julia -e "import Pkg; print(Pkg.dir(\"$PKGNAME\"))"`
rsync --no-o --no-g --no-perms -mav --include="*/" --include="*.cov" --exclude="*" $PKGDIR/ /mnt/
exit $TEST_EXIT
