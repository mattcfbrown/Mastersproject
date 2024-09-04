FROM julia:1.9-buster

LABEL image.author.name "Matthew Brown"
LABEL image.author.email "brownmc@student.unimelb.edu.au"

RUN apt-get update && apt-get install -y procps

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.add("Adapt");'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.add("NetworkInference");'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.develop(url="https://github.com/ananth-pallaseni/EmpiricalBayes.jl");'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.add("CSV");'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.add("DataFrames");'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.add("DelimitedFiles");'

RUN JULIA_DEPOT_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/share/julia \
    julia -e 'using Pkg; Pkg.add("Distributions");'



