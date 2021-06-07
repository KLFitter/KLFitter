FROM rootproject/root:6.24.00-ubuntu20.04

LABEL maintainer.name="KLFitter developer team"
LABEL maintainer.email="none"

# Set the BAT version here and propagate this throughout the script.
ARG BAT_VERSION=1.0.0
ARG BAT_SHA=184ff58f16fa35e73a4c48c9b218c2d8daac0a387e3400e8b2f995beb18f3b66

# Set the googletest version here and propagate this throughout the script.
ARG GTEST_VERSION=1.10.0
ARG GTEST_SHA=9dc9157a9a1551ec7a7e43daea9a694a0bb5fb8bec81235d8a1e6ef64c716dcb

USER root

WORKDIR /tmp

# Download, unpack and install BAT, then clean up thoroughly.
RUN wget https://github.com/bat/bat/archive/refs/tags/v${BAT_VERSION}.tar.gz \
    && echo "${BAT_SHA} v${BAT_VERSION}.tar.gz" | sha256sum --check --status \
    && tar xf v${BAT_VERSION}.tar.gz \
    && cd bat-${BAT_VERSION} \
    && ./autogen.sh \
    && ./configure --prefix=/opt/bat \
    && make -j 4 \
    && make install \
    && cd /tmp \
    && rm -rf bat-${BAT_VERSION} \
    && rm v${BAT_VERSION}.tar.gz

ENV PATH=${PATH}${PATH:+:}/opt/bat/bin
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:}/opt/bat/lib
ENV PKG_CONFIG_PATH=${PKG_CONFIG_PATH}${PKG_CONFIG_PATH:+:}/opt/bat/lib/pkgconfig

RUN wget https://github.com/google/googletest/archive/refs/tags/release-${GTEST_VERSION}.tar.gz \
    && echo "${GTEST_SHA} release-${GTEST_VERSION}.tar.gz" | sha256sum --check --status \
    && tar xf release-${GTEST_VERSION}.tar.gz \
    && mkdir build \
    && cd build \
    && cmake ../googletest-release-${GTEST_VERSION} -DBUILD_GMOCK=Off -DGTEST_HAS_PTHREAD=0 -DCMAKE_INSTALL_PREFIX=/opt/gtest \
    && make -j 4 \
    && make install \
    && cd /tmp \
    && rm -rf ./*

ENV GTEST_ROOT=/opt/gtest
ENV PATH=${PATH}${PATH:+:}/opt/gtest/bin
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:}/opt/gtest/lib
ENV PKG_CONFIG_PATH=${PKG_CONFIG_PATH}${PKG_CONFIG_PATH:+:}/opt/gtest/lib/pkgconfig

# Create a non-privileged user 'docker' which is used for running this image.
RUN groupadd -r docker -g 901 && useradd -u 901 -r -g docker docker

ENV HOME=/user/docker
USER docker
WORKDIR ${HOME}

ENTRYPOINT /bin/bash
CMD /bin/bash
