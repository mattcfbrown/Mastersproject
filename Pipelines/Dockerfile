FROM "amd64/r-base"

LABEL image.author.name "Matthew Brown"
LABEL image.author.email "brownmc@student.unimelb.edu.au"

RUN apt-get update
RUN apt-get install -y -t unstable libcurl4-openssl-dev
RUN apt-get install -y -t unstable libxml2-dev
RUN apt-get install -y -t unstable libglpk-dev

RUN Rscript -e 'install.packages("RCurl")'
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'install.packages("mclust")'
RUN Rscript -e 'install.packages("igraph")'
RUN Rscript -e 'BiocManager::install("slingshot", lib="/usr/local/lib/R/site-library", INSTALL_opts = "--no-staged-install")'
RUN Rscript -e 'library("RCurl")'
RUN Rscript -e 'library("igraph")'
RUN Rscript -e 'library("slingshot")'