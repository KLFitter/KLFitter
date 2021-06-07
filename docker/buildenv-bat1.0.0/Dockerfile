FROM rootproject/root:6.24.00-ubuntu20.04

LABEL maintainer.name="KLFitter developer team"
LABEL maintainer.email="none"

# Set the BAT version here and propagate this throughout the script.
ARG BAT_VERSION=1.0.0
ARG BAT_SHA=184ff58f16fa35e73a4c48c9b218c2d8daac0a387e3400e8b2f995beb18f3b66

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

# Create a non-privileged user 'docker' which is used for running this image.
RUN groupadd -r docker -g 901 && useradd -u 901 -r -g docker docker

ENV HOME=/user/docker
USER docker
WORKDIR ${HOME}

ENTRYPOINT /bin/bash
CMD /bin/bash
