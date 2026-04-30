#FROM gitlab-registry.cern.ch/delphi/deployment/delphi/al9_aarch64
FROM gitlab-registry.cern.ch/delphi/deployment/delphi/al9_64_xrd

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install GDB and common build/runtime deps useful for ROOT.
RUN if command -v microdnf >/dev/null 2>&1; then \
      microdnf install -y \
        gdb \
        gcc \
        gcc-c++ \
        make \
        cmake \
        git \
        wget \
        tar \
        bzip2 \
        xz \
        python3 \
        python3-pip \
        python3-devel \
        pkgconf-pkg-config \
        libX11 \
        libXpm \
        libXft \
        libXext \
        mesa-libGLU \
        mesa-libGL \
        libXrender \
        libXcursor \
        libXrandr \
        libXi \
        libXinerama \
        libuuid \
        libuuid-devel \
        openssl-devel \
        pcre-devel \
        zlib-devel \
        lz4 \
        zstd \
        glibc-langpack-en \
      && microdnf clean all; \
    else \
      dnf install -y \
        gdb \
        gcc \
        gcc-c++ \
        make \
        cmake \
        git \
        wget \
        tar \
        bzip2 \
        xz \
        python3 \
        python3-pip \
        python3-devel \
        pkgconf-pkg-config \
        libX11 \
        libXpm \
        libXft \
        libXext \
        mesa-libGLU \
        mesa-libGL \
        libXrender \
        libXcursor \
        libXrandr \
        libXi \
        libXinerama \
        libuuid \
        libuuid-devel \
        openssl-devel \
        pcre-devel \
        zlib-devel \
        lz4 \
        zstd \
        glibc-langpack-en \
      && dnf clean all; \
    fi

# Build and install CERN ROOT 6.34.00 from source (traditional CMake build).
ARG ROOT_VERSION=6.34.00
ARG ROOT_SRC_URL=https://root.cern/download/root_v6.34.00.source.tar.gz
RUN mkdir -p /tmp/root_src /tmp/root_build \
    && wget -qO /tmp/root_src.tar.gz "${ROOT_SRC_URL}" \
    && tar -xf /tmp/root_src.tar.gz -C /tmp/root_src --strip-components=1 \
    && rm -f /tmp/root_src.tar.gz \
    && cmake -S /tmp/root_src -B /tmp/root_build \
         -DCMAKE_INSTALL_PREFIX=/opt/root \
         -DCMAKE_BUILD_TYPE=Release \
         -Dgnuinstall=ON \
         -Dmathmore=ON \
         -Dpyroot=ON \
         -Dpython3=ON \
    && cmake --build /tmp/root_build -j"$(nproc)" \
    && cmake --install /tmp/root_build \
    && rm -rf /tmp/root_src /tmp/root_build

ENV ROOTSYS=/opt/root
ENV PATH=/opt/root/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/root/lib

# Install Boost C++ development libraries in a later layer.
RUN if command -v microdnf >/dev/null 2>&1; then \
      microdnf install -y boost-devel \
      && microdnf clean all; \
    else \
      dnf install -y boost-devel \
      && dnf clean all; \
    fi

# Build and install FastJet 3.3.3 (provides fastjet-config).
ARG FASTJET_VERSION=3.3.3
ARG FASTJET_SRC_URL=https://fastjet.fr/repo/fastjet-${FASTJET_VERSION}.tar.gz
RUN mkdir -p /tmp/fastjet_src /tmp/fastjet_build \
    && wget -qO /tmp/fastjet_src.tar.gz "${FASTJET_SRC_URL}" \
    && tar -xf /tmp/fastjet_src.tar.gz -C /tmp/fastjet_src --strip-components=1 \
    && rm -f /tmp/fastjet_src.tar.gz \
    && cd /tmp/fastjet_src \
    && ./configure --prefix=/opt/fastjet \
    && make -j"$(nproc)" \
    && make install \
    && rm -rf /tmp/fastjet_src /tmp/fastjet_build

ENV PATH=/opt/fastjet/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/fastjet/lib:$LD_LIBRARY_PATH

# Default to bash for interactive work.
CMD ["/bin/bash"]
