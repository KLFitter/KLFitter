>>>
Please note that this README file was written in 2011 and essentially has not
been updated since. Expect a major update to happen soon.
>>>


# The (K)inematic (L)ikelihood Fitter


At the time, this README was written, these were the main authors of KLFitter:

| Author           | Affiliation              | Email                                        |
| ---------------- | ------------------------ | -------------------------------------------- |
| Kevin Kroeninger | University of Goettingen | kevin.kroeninger *AT* phys.uni-goettingen.de | 
| Johannes Erdmann | University of Goettingen | johannes.erdmann *AT* phys.uni-goettingen.de |
| Olaf Nackenhorst | University of Goettingen | olaf.nackenhorst *AT* phys.uni-goettingen.de |


## Outline

1. Purpose / motivation
2. Installation
3. Class structure
4. The fitting procedure
5. Using the KLFitter
6. How to run KLFitter - an example 


## 1. Purpose / motivation

The K(inematic)L(ikelihood)Fitter is a library for kinematic fitting
using a likelihood approach. It is primarily developed for the case of
top quark reconstruction, but it can be easily modified to fit other
processes. KLFitter is experiment independent in a sense that
different experiments can be parameterized.


## 2.Installation

### 2.1 Dependencies

KLFitter depends on ROOT and the BAT library. BAT can be obtained
from the following web page:

http://www.mppmu.mpg.de/bat/ 

Please read the instructions on
how to install and run BAT. 

### 2.2 Installation

The installation procedure is simple: 

Check out the branch `master` or a tagged release of KLFitter from the [KLFitter
  git repository](https://gitlab.cern.ch/KLFitter/KLFitter/) by typing:

```
$ git clone https://gitlab.cern.ch/KLFitter/KLFitter.git
$ cd KLFitter
# For checking out a tagged release:
$ git checkout VERSION_TAG
```

Make sure that the environment variable `$BATINSTALLDIR` is correctly defined.
  In step 2.1, you should have installed BAT to a specific location, so export
  that location with

```
$ export BATINSTALLDIR=/path/to/BAT/install
```

Compile the KLFitter library - this will produce the standard library (classes
  in the folder 'lib') you can use in your analysis. Compilation is
  straightforward and can be done with

```
$ make -j
```

Additionally, in case you want to install KLFitter to a designated directory,
  you can create a _clean_ KLFitter install version that only includes the
  header and library files:
  
```
$ make install
```

The files will be copied into the folder 'dest-tmp'.

KLFitter also comes with lots of examples. These example implementations and
instructions on how to use them, please refer to the
[KLFitterExtras](https://gitlab.cern.ch/KLFitter/KLFitterExtras) repository.


## 3. Class structure

The KLFitter package is build in a modular way so as to allow the
implementation of different physics processes and parameterizations of
the detector response. All classes live in the "KLFitter" namespace. The
class structure is as follows:

The central class is "Fitter". It contains objects which describe the
detector (an instance of "DetectorBase"), the input particles (an
instance of "Particles") and the likelihood (an instance of
"LikelihoodBase"). An additional object, an instance of
"Permutations", helps to manage the possible permutations of jets to
quarks, and leptons.

The DetectorBase class summarizes the detector information. It
contains several objects (instances of "ResolutionBase") describing
the parameterizastion of, e.g., energy resolutions for different
particle types. The description of a detector should be done by
creating a class which inherits from DetectorBase. The resolution
objects can be created by additional classes which inherit from
ResolutionBase. An example exists for a dummy detector with Gaussian
energy resolutions.

The Particle class holds containers for different types of particles,
partons (a common class for all quarks and gluons), electron, muons, taus,
neutrinos, photons and boson (a common class for photons and W/Z
bosons). Parton-level or reconstructed particles can be described
here. 

The actual fitting procedure is done in the LikelihoodBase class. This
class inherits from the BAT BCModel class. In order to implement a new
process, the user has to create a class which inherits from
LikelihoodBase. The user has to define the parameters of the fit and
the likelihood itself. One example is the LikelihoodTopLeptonJets class for l+jets ttbar channel.


## 4. The fitting procedure

### 4.1 The likelihood

The KLFitter is based on a likelihood approach. The likelihood
function has to be specified separately for each process. It can
consist of any function defined by the user. Typically, it is the
product of some parameterization of the energy resolution of the
measured jets and leptons, and functions describing the physics
process, such as Breit-Wigner functions, matrix elements, etc. 

#In case of semileptonic ttbar events, the likelihood is the product of
#four quark energy resolution functions, four quark eta/phi resolution 
#functions one lepton energy resolution function, neutrino pT resolution 
#function and the Breit-Wigner functions of the W boson masses. Optionally, 
#the Breit-Wigner functions of the top quarks can be used either with a free
#or a fixed top pole mass parameter. 

In case of semileptonic ttbar events, the likelihood is the product of
four quark energy resolution functions, one lepton energy resolution function,
neutrino pT resolution function and the Breit-Wigner functions of the W boson
masses. Optionally, the Breit-Wigner functions of the top quarks can be used
either with a free or a fixed top pole mass parameter. 

In this example, the parameters which are estimated are: 
a) the energies of the four jets
b) the energy of the charged lepton
c) the x,y,z-components of the neutrino momentum
d) optionally, the top pole mass parameter

This likelihood is implemented in the class LikelihoodTopLeptonJets.

There is an extended version implemented in the class
LikelihoodTopLeptonJets_JetAngles, which also fits the angular
variables of the jets. Resolutions functions are multiplied to the
likelihood for eta and phi of the jets and there are, hence, eight
additional parameters which are estimated in the fit.
Please note that this likelihood is about a factor of two slower
compared to the simpler LikelihoodTopLeptonJets, because a lot of
trignometric functions need to be evaluated in each call of the
likelihood.

Other parameters can easily be added. 

### 4.2 Optimization

Optimization is done using the methods provided by BAT. By default,
the optimization method used is TMinuit.

### 4.3 Marginalization

Marginalization is done using the methods provided by BAT. By default,
this method is Markov Chain Monte Carlo. 

### 4.4 Default analysis

In the default analysis a Markov chain is run to coarsly sample the
parameter space and locate a global maximum. Minuit is then used to
find the maximum with a greater precision. 


## 5. Using the KLFitter

### 5.1 Setting the measured quantities

The measured quantities are the 4-vectors of four jets and one charged
lepton, as well as the missing transverse energy. An object of type 
Particles is passed to the Fitter object. The input quantities can 
either be set by defining TLorentzVectors and adding them to Particles 
object, or by using an interface. An example for a Root interface exists 
which reads the data from a flat Root tree.

### 5.2 Combinatorics

The association of jets with partons and of reconstructed and parton
level leptons leads to several possible combinations. The class
Permutations calculates a table with all possible permutations. A
specific combination (index) can be set by the user,

int KLFitter::Permutations::SetPermutation(int index). 

In case interchanging two (or three) particles leaves the likelihood
invariant two (three) indices can be defined which are interchangeable
(for example for hadronically decaying W bosons). The corresponding
permutations are removed from the table. The indices can be set via

int InvariantPartonPermutations(int index1, int index2, int index3 = -1). 

### 5.3 Output

The output in the provided examples comes in form of a ROOT file which contains the following
trees: 

* `TreeMeasured`: the measured particles.
* `TreeSelected`: the particles selected for the fitting.
* `TreeModel`: the model particles (e.g., including also the top quark and other
  intermediate particles).
* `TreeTruth`: the Monte Carlo truth information (if avialable).
* `TreeMap`: maps containing the mapping between the measured and the selected
  particles.

The TreeModel contains the results from the fit. The most important
stored variables are:

* `EventNumber`: the event number.
* `N_permutations`: the number of permutations.
* `best_permutation[N_permutations]`: a list of indices ordered by the
  EventProbability. The best permutation has the index best_permutation[0], etc.
* `EventProbability[N_permutations]`: the event probability for each
  permutation.
* `LogLikelihood[N_permutations]`: the log Likelihood for each permutation.
* `ConvergenceStatusBit[N_permutations]`: the convergence status bit word for
  each permutation.
* `x_E[N_permutations]`, `x_pt`: x is a particle known in the model (i.e.,
  initial state or intermediate particle, etc.). The `_E`, ... describe the
  kinematics for each permutation.
* `par_x[N_permutations]`: x is the best fit parameter for each permutation.
* `parerr_x[N_permutations]`: x is the estimated uncertainty on the best fit
  parameter for each permutation.

## 6. How to run KLFitter - an example

A example on how to use the KLFitter library is provided: ExampleCode.C.
The example is extensively commented. However, it is kept "minimal" in order to
briefly show the necessary steps to get the KLFitter running in your analysis.
It will hence not run out-of-the-box itself (for working example consult the example
directory). So, just go ahead and try it in your own analysis code.
