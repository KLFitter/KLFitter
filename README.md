# KLFitter – The Kinematic Likelihood Fitter

[![pipeline status](https://gitlab.cern.ch/KLFitter/KLFitter/badges/master/pipeline.svg)](https://gitlab.cern.ch/KLFitter/KLFitter/commits/master)

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
   [COPYING](COPYING) and [COPYING.lesser](COPYING.LESSER) for the GNU General
   Public License and the additional terms of the GNU Lesser General Public
   License, respectively.


### Dependencies

KLFitter depends on the ROOT and BAT libraries. For information about ROOT,
please consult the [ROOT webpage](https://root.cern.ch/). BAT releases and
information about BAT and its installation can be found on the
[library's webpage](http://www.mppmu.mpg.de/bat/). The following versions have
been tested with KLFitter and are working:
 - ROOT v5.34.10 and above
 - BAT v0.9.4.1


### Obtaining KLFitter

The KLFitter source code can be obtained from this repository 
([https://gitlab.cern.ch/KLFitter/KLFitter/]). A list of all releases of
KLFitter can be found under 
[KLFitter/tags](https://gitlab.cern.ch/KLFitter/KLFitter/tags). To download the
source code, you can use the git clone command:

```
$ git clone https://gitlab.cern.ch/KLFitter/KLFitter.git
$ cd KLFitter
  (for checking out a specific release)
$ git checkout VERSION_TAG
```

### Installation via cmake (recommended)

Before proceeding please make sure that you have a valid installation of ROOT
on your system and that cmake is correctly configured (version 3.1 or higher is
required). As cmake does fully automatised configuration, it is possible to
build KLFitter and the BAT library simultaneously. BAT will be downloaded during
the cmake build process. For doing so, please change into the KLFitter
directory and then do:

```
mkdir build && cd build
cmake -DBUILTIN_BAT=TRUE ..
make -j
```

This will build both BAT and KLFitter in the _build_ subfolder of the KLFitter
directory. If you already have an existing installation of BAT, you can also
link the KLFitter library against that version (and not download BAT during the
cmake build process). For this, make sure that the environment variable 
`$BATINSTALLDIR` points to the BAT installation directory and is exported. This
variable will be used by the _FindBAT.cmake_ that locates the BAT library. With
a local version of BAT, the `-DBUILTIN_BAT` flag of the cmake command can be 
omitted, i.e. the KLFitter cmake command becomes `cmake ..`.

If you encounter problems with the cmake configuration, these could be possible
reasons:
- If ROOT cannot be found by cmake, include the ROOT binary path to the `PATH`
variable and make sure it is exported and can be picked up by cmake.
- If you use a local version of BAT (no download via cmake) and it cannot be
found, make sure to set the `$BATINSTALLDIR` variable and export it. This
variable is used by cmake to locate the library.
- If linkage against ROOT or BAT fails: are those libraries compiled with the
same version of the compiler? Make sure that this is consistent. If you set up
a custom compiler for the build process, maybe cmake doesn't pick this compiler
up correctly. Do `whereis c++` and `which c++` point to the same binary? To
tell cmake to use the latter one, you can export the compiler location:

```
export CXX=`which c++`
export CC=`which gcc`
```


### Installation via Makefile (deprecated)

The repository also comes with a _Makefile_, although using cmake to build
KLFitter is the recommended procedure. Compilation via Makefile assumes that
you have a working installation of ROOT and BAT. The location of the BAT library
is determined with the `$BATINSTALLDIR` variable, so make sure to set and export
it. Then switch to the KLFitter repository directory and call

```
make -j && make -j install
```

to compile KLFitter. The latter command will create a subdirectory _build_ and
copy the header and library files into that directory for external usage.
