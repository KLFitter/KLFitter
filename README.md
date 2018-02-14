# KLFitter – The Kinematic Likelihood Fitter

[![pipeline status](https://gitlab.cern.ch/KLFitter/KLFitter/badges/master/pipeline.svg)](https://gitlab.cern.ch/KLFitter/KLFitter/commits/master)
[![coverage report](https://gitlab.cern.ch/KLFitter/KLFitter/badges/master/coverage.svg)](https://gitlab.cern.ch/KLFitter/KLFitter/commits/master)

KLFitter is a library for kinematic fitting using a likelihood approach. It is
primarily developed for the case of top quark reconstruction, but it can be
easily modified to fit other processes. KLFitter is experiment-independent,
that is, different experiments can be parameterized.

 - [List of authors and contributors](doc/Authors.md) to KLFitter
 - [What is KLFitter and how does it work?](doc/WhatIsKLF.md)
 - For implementation and examples and instructions on how to use KLFitter in
   a software framework, please check out the
   [KLFitterExtras repository](https://gitlab.cern.ch/KLFitter/KLFitterExtras).
 - KLFitter is licensed under the GNU Lesser General Public License v3.0. More
   information about the licensing terms and conditions can be found under
   [COPYING](COPYING) for the GNU General Public License and under
   [COPYING.lesser](COPYING.LESSER) for the additional terms of the GNU Lesser
   General Public License.


## Dependencies

KLFitter depends on the ROOT and BAT libraries. For information about ROOT,
please consult the [https://root.cern.ch/](Root webpage). BAT releases and
information about BAT and its installation can be found on the
[http://www.mppmu.mpg.de/bat/](library's webpage). The following versions have
been tested with KLFitter and are working:
 - ROOT v5.34.10 and above
 - BAT v0.9.4.1 and above


## Obtaining KLFitter

The KLFitter source can be obtained from this
[KLFitter git repository](https://gitlab.cern.ch/KLFitter/KLFitter/). A list of
all tagged releases can be found under 
[KLFitter/tags](https://gitlab.cern.ch/KLFitter/KLFitter/tags). To download the
source code, you can use git via

```
$ git clone https://gitlab.cern.ch/KLFitter/KLFitter.git
$ cd KLFitter
# For checking out a tagged release:
$ git checkout VERSION_TAG
```

## Installation via cmake

Before proceeding please make sure that you have a valid installation of ROOT
on your system and that cmake is correctly configured (version 3.1 or higher is
required). As cmake does fully automatised configuration, it is possible to
build KLFitter and the BAT library simultaneously. For doing so, please type
in your KLFitter directory:

```
mkdir build && cd build
cmake -DBUILTIN_BAT=TRUE ..
make -j
```

This will build both BAT and KLFitter in the `build` subfolder of the KLFitter
directory. If you already have an existing installation of BAT, you can also
link the KLFitter library against that version. For this, make sure that the
environment variable `$BATINSTALLDIR` points to the BAT installation directory
and is exported. Then the KLFitter cmake procedure is

```
mkdir build && cd build
cmake ..
make -j
```

## Installation via Makefile

to be written ...