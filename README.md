# KLFitter – The Kinematic Likelihood Fitter

[![Build Status](https://img.shields.io/travis/KLFitter/KLFitter/dev-2.0)](https://travis-ci.org/KLFitter/KLFitter)
[![Development Branch](https://img.shields.io/badge/dev%20branch-2.0-red)](https://github.com/KLFitter/KLFitter/branches)
[![GitHub release](https://img.shields.io/github/release/KLFitter/KLFitter.svg)](https://github.com/KLFitter/KLFitter/releases)
[![Github commits (since latest release)](https://img.shields.io/github/commits-since/KLFitter/KLFitter/latest.svg)](https://github.com/KLFitter/KLFitter/commits/)


KLFitter is a library for kinematic fitting using a likelihood approach. It is
primarily developed for the case of top quark reconstruction, but it can be
easily modified to fit other processes. Detailed documentation:

- [List of authors and contributors](doc/Authors.md)
- [What is KLFitter and how does it work?](doc/WhatIsKLF.md)
- [Out-of-the-box example](doc/Example.md)
- [Frequently Asked Questions](doc/FAQ.md)

An auto-generated class reference guide for the KLFitter library is also
uploaded to https://KLFitter.github.io. The uploaded version always corresponds
to the latest version of the master branch.


### Licensing terms

KLFitter is licensed under the GNU Lesser General Public License v3.0. For more
information about the licensing terms and conditions, please refer to
[COPYING](COPYING) and [COPYING.lesser](COPYING.LESSER) for the GNU General
Public License and the additional terms of the GNU Lesser General Public
License, respectively.


### Dependencies

KLFitter depends on the ROOT and BAT libraries. For information about ROOT,
please consult the [ROOT web page](https://root.cern.ch/). Information about
BAT, installation instructions and the latest BAT releases can be found on the
[library's web page](http://www.mppmu.mpg.de/bat/). The following versions have
been tested with KLFitter and are working:

- ROOT v6.24.00
- BAT v0.9.4.1 (developments for compatibility with v1.0.0 ongoing)


### Obtaining KLFitter

The KLFitter source code can be obtained from this repository
(https://github.com/KLFitter/KLFitter). A list of all releases of KLFitter can
be found under
[KLFitter/releases](https://github.com/KLFitter/KLFitter/releases). To download
the source code, you can use the git clone command:

```sh
git clone https://github.com/KLFitter/KLFitter.git
```

and, for checking out a specific release, use the git checkout procedure:

```sh
cd KLFitter
git checkout $VERSION_TAG
```

Before KLFitter was made public, it was maintained in non-public SVN
repositories. For completeness we also keep a list of all SVN tags of KLFitter
and their associated commits in [this file](doc/SVN-tags.md).


### Installation via cmake

Before proceeding please make sure that you have a valid installation of ROOT on
your system and that cmake is correctly configured (version 3.1 or higher is
required). An installed version of BAT is recommended, too. However, as cmake
does fully automated configuration, it is possible to build KLFitter and the BAT
library simultaneously. In that case, BAT will be downloaded during the cmake
build process. For doing so, create a build directory and call the cmake command
with the `KLF_BUILD_BAT` flag set:

```sh
mkdir build && cd build
cmake -DKLF_BUILD_BAT=TRUE /path/to/KLFitter/source
make -j
```

This will build both BAT and KLFitter in the _build_ directory. Please be aware
that with this setup KLFitter can only be run inside the build tree, but an
attempt to `make install` KLFitter into permanent installation tree will not
copy over the on-the-fly BAT build!

NB: For executing KLFitter inside the build tree, you might have to -- depending
on your system -- update the linker paths in order to find the KLFitter shared
library, e.g. by doing:

``` shsh
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:}/path/to/build/tree/lib
```

For a permanent installation, BAT and KLFitter should be installed separately,
and the `KLF_BUILD_BAT` flag should not be used. With a pre-installed version of
BAT to link against, you can omit the `KLF_BUILD_BAT` flag and build against BAT
directly. cmake will attempt to find the BAT installation through its package
configuration files (if that fails, you may edit the file `CMakeLists.txt` and
set the variable `BAT_ROOT` to the root directory of the BAT installation, e.g.
`/usr/local`, which will be used as a hint where to look for BAT). With a
pre-built BAT, the KLFitter installation procedure simply becomes:

```sh
mkdir build && cd build
cmake /path/to/KLFitter/source
make -j
```

And then you may install it via `make install` which, on a UNIX-based system,
installs to `/usr/local` unless `CMAKE_BUILD_PREFIX` is set to a different value
when calling the cmake command.

For problems with the cmake configuration, please also refer to the [Frequently
Asked Questions](doc/FAQ.md) and don't hesitate to get in touch!


### Installation in a docker image

> TODO: Eventually we would like to provide a docker image with a pre-installed
> version of KLFitter installed that is deployed automatically with the CI
> runners, but this isn't set up yet.

KLFitter may also be installed inside a docker image that already provides all
necessary dependencies to build and run KLFitter. If you wish to do so, first
download the required docker image:

``` sh
docker pull klfitter/klfitter-buildenv:ubuntu20.04-root6.24.00-bat0.9.4.1
```

Once the image is downloaded you can start an interactive docker container based
on it with the following command:

``` sh
docker run -i -d --name klf-docker -v /path/to/KLFitter/source:/user/docker/src -w /user/docker/build klfitter/klfitter-buildenv:ubuntu20.04-root6.24.00-bat0.9.4.1
```

where you should replace `/path/to/KLFitter/source` with the actual path to the
KLFitter source code on your machine. This starts the container, mounts the
KLFitter source code to the path `/user/docker/src` inside the docker image and
sets the docker working directory to be `/user/docker/build`. 

You can now repeat the above cmake configuration steps by preceeding every
command with `docker exec klf-docker` and by specifying the mount point of the
KLFitter source code, e.g. by executing the following commands:

``` sh
docker exec build cmake /user/docker/src
docker exec build make -j
```

Once that is done, the KLFitter shared library and the example executable should
be available inside the docker container under `/user/docker/build/lib` and
`/user/docker/build/bin`, respectively.
