################################################################################
# Dockerfile to build solids4foam with foam-extend-4.0 and Ubuntu-18.04-LTS
# philip.cardiff@gmail.com, July 2019
################################################################################

# Set the base image to foam-extend-4.0 and Ubuntu-18.04-LTS
FROM philippic/foam-extend-4.0-ubuntu18.04

# File Author / Maintainer
MAINTAINER Philip Cardiff <philip.cardiff@gmail.com>

# Change default shell to bash
SHELL ["/bin/bash", "-c"]

# Change to root and update the repository sources list and then change back to
# the app user
USER root
RUN apt update && apt install -y sudo
USER app

# Set working directory to the app home directory
WORKDIR /home/app

####################### solids4Foam INSTALLATION ###############################

# Create FOAM_RUN directory
RUN echo "Creating FOAM_RUN directory" && \
    source /home/app/foam/foam-extend-4.0/etc/bashrc && \
    mkdir -p $FOAM_RUN

# Add solids4Foam and change permissions
RUN echo "Adding and compiling solids4Foam-4.0"
COPY . /home/app/foam/app-4.0/solids4foam
USER root
RUN chown -R app: /home/app/foam/app-4.0/solids4foam
USER app

# Compile solids4Foam
RUN echo "Compiling solids4Foam" && \
    source /home/app/foam/foam-extend-4.0/etc/bashrc && \
    cd $FOAM_RUN/../solids4foam && \
    ./Allwmake