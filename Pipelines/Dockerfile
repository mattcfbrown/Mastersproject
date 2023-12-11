FROM "rocker/r-base"

LABEL image.author.name "Matthew Brown"
LABEL image.author.email "brownmc@student.unimelb.edu.au"

RUN R -e "install.packages('igraph',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('nlnet',dependencies=TRUE, repos='http://cran.rstudio.com/')"

ENV PATH=$PATH:/usr/games/