FROM kyanagis/restretto:1.0

SHELL ["/bin/bash", "-c"]

# To suppress interactive installation of boost
ARG DEBIAN_FRONTEND=noninteractive

# post processing
RUN echo "alias ls='ls --color=auto'" >> /etc/bash.bashrc

