FROM "arm64v8/r-base"

LABEL image.author.name "Matthew Brown"
LABEL image.author.email "brownmc@student.unimelb.edu.au"

RUN Rscript -e 'install.packages("Hmisc")'

ENV PATH=$PATH:/usr/games/