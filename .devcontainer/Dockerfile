FROM kyanagis/restretto:1.0

# To suppress interactive installation of boost
ARG DEBIAN_FRONTEND=noninteractive

# post processing
RUN echo "alias ls='ls --color=auto'" >> /etc/bash.bashrc

# ====== settings for git / ssh =======
RUN echo "mkdir -p /root/.ssh" >> /etc/bash.bashrc
RUN echo "cp /root/.ssh_host/* /root/.ssh/" >> /etc/bash.bashrc
RUN echo "chmod 600 /root/.ssh/*" >> /etc/bash.bashrc