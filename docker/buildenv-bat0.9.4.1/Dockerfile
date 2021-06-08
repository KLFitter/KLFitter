# Copyright (c) 2009--2021, the KLFitter developer team
#
# This file is part of KLFitter.
#
# KLFitter is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# KLFitter is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
#
FROM rootproject/root:6.24.00-ubuntu20.04

LABEL maintainer.name="KLFitter developer team"
LABEL maintainer.email="none"

# Set the BAT version here and propagate this throughout the script.
ARG BAT_VERSION=0.9.4.1
ARG BAT_SHA=0f363bfbebc6902f451580a141dfdcba0edeaf7022e50fa88df014750f06cb74

# Set the googletest version here and propagate this throughout the script.
ARG GTEST_VERSION=1.10.0
ARG GTEST_SHA=9dc9157a9a1551ec7a7e43daea9a694a0bb5fb8bec81235d8a1e6ef64c716dcb

USER root

WORKDIR /tmp

# We need to apply some minor patches to BAT in order to avoid a stream of
# deprecation warnings and an actual compilation error due to a missing include.
RUN apt-get update \
    && apt-get install -y patch

COPY BAT-0.9.4.1-cxx17.patch BAT-0.9.4.1-cxx17.patch
COPY BAT-0.9.4.1-root624.patch BAT-0.9.4.1-root624.patch

# Download, unpack and install BAT, then clean up thoroughly.
RUN wget https://github.com/bat/bat/archive/refs/tags/v${BAT_VERSION}.tar.gz \
    && echo "${BAT_SHA} v${BAT_VERSION}.tar.gz" | sha256sum --check --status \
    && tar xf v${BAT_VERSION}.tar.gz \
    && cd bat-${BAT_VERSION} \
    && patch -ruN -p1 < /tmp/BAT-0.9.4.1-cxx17.patch \
    && patch -ruN -p1 < /tmp/BAT-0.9.4.1-root624.patch \
    && ./autogen.sh \
    && ./configure --prefix=/opt/bat \
    && make -j 4 \
    && make install \
    && cd /tmp \
    && rm -rf ./*

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
