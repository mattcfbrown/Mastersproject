FROM julia:1.9-buster

LABEL image.author.name "Matthew Brown"
LABEL image.author.email "brownmc@student.unimelb.edu.au"

RUN apt-get update && apt-get install -y procps

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.update(); Pkg.add("Adapt"); Pkg.precompile();'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.update(); Pkg.add("DelimitedFiles"); Pkg.precompile();'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.update(); Pkg.add("NetworkInference"); Pkg.precompile();'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.update(); Pkg.add(url="https://github.com/ananth-pallaseni/EmpiricalBayes.jl"); Pkg.precompile();'