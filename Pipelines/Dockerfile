FROM "r-base"

LABEL image.author.name "Matthew Brown"
LABEL image.author.email "brownmc@student.unimelb.edu.au"

RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install("slingshot")'
RUN Rscript -e 'install.packages("mclust")'

ENV PATH=$PATH:/usr/games/
ENV DOCKER_DEFAULT_PLATFORM=linux/amd64