# Build from intel-hpckit image
FROM intel/hpckit:2024.0.1-devel-ubuntu22.04

# Intel did not provide new key (23.05.2024), delete it from source to
# avoid errors.
RUN rm /etc/apt/sources.list.d/intel-graphics.list

# Use apt-get to install all necessary things
RUN apt-get update
RUN apt-get install -y vim
RUN apt-get install -y less
RUN apt-get install -y libsm6 libxrender1 libfontconfig1
RUN apt-get install -y unzip

# Make it possible to also have the hdf5 output
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y autoconf
RUN apt-get install -y libtool

# Add a new user without root access, this is beneficial when
# mounting folders to the docker as it will avoid requiring root
# access to open from outside the container.
RUN useradd -g users winnet && passwd -d winnet

# Create the home folder and change permissions
RUN mkdir /home/winnet

# Modify Matplotlib default params
COPY matplotlibrc /home/winnet/.config/matplotlib/
RUN chown -R winnet:users /home/winnet/.config/matplotlib/

# Change the access rights
RUN chown winnet:users /home/winnet/

# Take care of FFMPEG and all necessary to run the movie script
RUN apt-get install -y ffmpeg

# Switch the user
USER winnet:users

# Clone hdf5 library and compile it with the intel compiler
RUN cd /home/winnet && git clone https://github.com/HDFGroup/hdf5.git
RUN cd /home/winnet/hdf5 && ./autogen.sh
RUN cd /home/winnet/hdf5 && ./configure --enable-fortran CC=icx FC=ifx CXX=icpx
RUN cd /home/winnet/hdf5 && make -j8
RUN cd /home/winnet/hdf5 && make install
ENV PATH="/home/winnet/hdf5/hdf5/bin:$PATH"

# Install miniconda
RUN mkdir -p /home/winnet/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/winnet/miniconda3/miniconda.sh
RUN bash /home/winnet/miniconda3/miniconda.sh -b -u -p /home/winnet/miniconda3
RUN rm -rf /home/winnet/miniconda3/miniconda.sh
RUN /home/winnet/miniconda3/bin/conda init bash
RUN /home/winnet/miniconda3/bin/conda init zsh
RUN echo 'export PATH="/home/winnet/miniconda3/bin:$PATH"' >> /home/winnet/.bashrc
ENV PATH="/home/winnet/miniconda3/bin:$PATH"

# Clone WinNet
RUN cd /home/winnet && git clone https://github.com/nuc-astro/WinNet.git
# Setup the makefile to correctly use ifx in order to get rid of annoying warnings. Also delete all unnecessary stuff
RUN cd /home/winnet/WinNet && mv Makefile.example_linux Makefile && sed 's/FC\ =\ ifort/FC\ =\ h5fc/g' -i Makefile && rm Makefile.example_mac && mv makerun.py.example makerun.py

# Setup conda environment for WinNet
RUN conda install -y pip
# Now install all needed packages
RUN cd /home/winnet/WinNet/install && pip install -r requirements.txt

# Install everything for the movie
RUN pip install mpi4py
ENV PATH="/home/winnet/WinNet/bin/movie_script:$PATH"


# Use the WinNet directory as default one
RUN echo 'cd /home/winnet/WinNet/' >> /home/winnet/.bashrc
RUN echo 'export OMP_NUM_THREADS=1' >> /home/winnet/.bashrc
RUN echo "export LS_OPTIONS='--color=auto'" >> /home/winnet/.bashrc
RUN echo 'eval "$(dircolors -b)"' >> /home/winnet/.bashrc
RUN echo "alias ls='ls \$LS_OPTIONS'" >> /home/winnet/.bashrc
