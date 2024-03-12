FROM "arm64v8/r-base"

LABEL image.author.name "Matthew Brown"
LABEL image.author.email "brownmc@student.unimelb.edu.au"

RUN Rscript -e 'install.packages("fdrtool")'
RUN Rscript -e 'install.packages("MASS")'

ENV PATH=$PATH:/usr/games/