# KLFitter – The Kinematic Likelihood Fitter

[![CI](https://github.com/KLFitter/KLFitter/actions/workflows/ci.yml/badge.svg)](https://github.com/KLFitter/KLFitter/actions/workflows/ci.yml)
[![GitHub release](https://img.shields.io/github/release/KLFitter/KLFitter.svg)](https://github.com/KLFitter/KLFitter/releases)
[![Github commits (since latest release)](https://img.shields.io/github/commits-since/KLFitter/KLFitter/latest.svg)](https://github.com/KLFitter/KLFitter/commits/)


KLFitter is a library for kinematic fitting using a likelihood approach. It is
primarily developed for the case of top quark reconstruction, but it can be
easily modified to fit other processes. Detailed documentation:

- [List of authors and contributors](doc/Authors.md)
- [What is KLFitter and how does it work?](doc/WhatIsKLF.md)
- [Out-of-the-box example](doc/Example.md)
- [Frequently Asked Questions](doc/FAQ.md)

An auto-generated class reference guide for the KLFitter library is also uploaded to https://KLFitter.github.io/KLFitter/. The uploaded version always corresponds to the latest version of the master branch.


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

- ROOT v5.34.10 or later
- BAT v0.9.3 and v0.9.4.1


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
your system and that cmake is correctly configured (version 3.10 or higher is
required).

The standard way to build KLFitter is against an existing BAT installation.
CMake will automatically search for BAT in standard system paths (e.g.,
`/usr/local`). If BAT is installed in a non-standard location, you should
provide its path either by setting the `BAT_ROOT` CMake variable or the
`$BAT_ROOT` (or legacy `$BATINSTALLDIR`) environment variable.

The build procedure is:

```sh
mkdir build && cd build
cmake .. # Add -DBAT_ROOT=/path/to/bat if needed
make -j
```

Alternatively, it is possible to build KLFitter and the BAT library
simultaneously by letting cmake download and compile BAT automatically. For
doing so, use the `-DBUILTIN_BAT=TRUE` flag:

```sh
mkdir build && cd build
cmake -DBUILTIN_BAT=TRUE ..
make -j
```

This will build both BAT and KLFitter in the sub-folder of the KLFitter
directory. For problems with the cmake configuration, please also refer to the
[Frequently Asked Questions](doc/FAQ.md).
