# REST2 Enhanced Sampling Docker Environment
# Based on GROMACS and PLUMED

FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    python3 \
    python3-pip \
    python3-dev \
    libfftw3-dev \
    liblapack-dev \
    libblas-dev \
    libgsl-dev \
    libboost-all-dev \
    libeigen3-dev \
    libhdf5-dev \
    libnetcdf-dev \
    libx11-dev \
    libxext-dev \
    libxrender-dev \
    libxrandr-dev \
    libxss-dev \
    libxtst-dev \
    libxi-dev \
    libxfixes-dev \
    libxcb1-dev \
    libxcb-render0-dev \
    libxcb-shape0-dev \
    libxcb-xfixes0-dev \
    libxcb-randr0-dev \
    libxcb-xtest0-dev \
    libxcb-xinerama0-dev \
    libxcb-xkb-dev \
    libxkbcommon-dev \
    libxkbcommon-x11-dev \
    libxrandr-dev \
    libxss-dev \
    libxtst-dev \
    libxi-dev \
    libxfixes-dev \
    libxcb1-dev \
    libxcb-render0-dev \
    libxcb-shape0-dev \
    libxcb-xfixes0-dev \
    libxcb-randr0-dev \
    libxcb-xtest0-dev \
    libxcb-xinerama0-dev \
    libxcb-xkb-dev \
    libxkbcommon-dev \
    libxkbcommon-x11-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt

# Install GROMACS (example for 2023.2)
WORKDIR /opt
RUN wget https://ftp.gromacs.org/pub/gromacs/gromacs-2023.2.tar.gz && \
    tar -xzf gromacs-2023.2.tar.gz && \
    cd gromacs-2023.2 && \
    mkdir build && cd build && \
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && \
    rm -rf gromacs-2023.2*

# Install PLUMED (example for 2.8.3)
RUN wget https://github.com/plumed/plumed2/releases/download/v2.8.3/plumed-2.8.3.tgz && \
    tar -xzf plumed-2.8.3.tgz && \
    cd plumed-2.8.3 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf plumed-2.8.3*

# Set up environment
ENV PATH="/usr/local/gromacs/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/gromacs/lib:${LD_LIBRARY_PATH}"
ENV PKG_CONFIG_PATH="/usr/local/gromacs/lib/pkgconfig:${PKG_CONFIG_PATH}"

# Create working directory
WORKDIR /workspace

# Copy project files
COPY . /workspace/

# Make main script executable
RUN chmod +x /workspace/main.py

# Set default command
CMD ["python3", "main.py", "--help"] 