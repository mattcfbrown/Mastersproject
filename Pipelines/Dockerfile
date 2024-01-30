FROM "rocker/r-base"

LABEL image.author.name "Matthew Brown"
LABEL image.author.email "brownmc@student.unimelb.edu.au"

RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install("GENIE3")'

ENV PATH=$PATH:/usr/games/