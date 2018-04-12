# Example of a KLFitter implementation

### Introduction

The purpose of the provided example is to demonstrate how to implement KLFitter
in your analysis. The example code and the input files are designed to provide a
minimalistic input and implementation. The example also shows how to retrieve
the output from KLFitter after the fitting. The ttbar lepton+jets likelihood,
implemented in the class
[LikelihoodTopLeptonJets](include/KLFitter/LikelihoodTopLeptonJets.h) has been
chosen for the example, but from the implementation it should be straightforward
to extend the implementation to other available likelihoods. Detailed
explanations of all likelihoods can be found in the [main KLFitter
documentation](doc/WhatIsKLF.md).


### How to run the example

To run the example, compile the KLFitter library first. For build and
installation instructions, please refer to the [README file](README.md). Then
run the executable, with the path to the KLFitter source directory as an
additional argument. For example, if you followed the build instructions of the
README, the cmake command will store the executable in the `bin` subdirectory.
Then call

```
$ ./bin/example-top-ljets.exe ../
```

The location of the KLFitter source directory is needed to locate the input file
for the example, stored under [data/examples](data/examples), and the transfer
functions, stored under [data/transferfunctions](data/transferfunctions).


### Implementation

The implementation of the example consists of three different files:

* An [input ROOT file](data/examples/top-ljets-input.root) with a few hundred
ttbar lepton+jets events. The file contains all necessary variables to run the
KLFitter lepton+jets likelihood (described in detail in the [main
documentation](doc/WhatIsKLF.md)).
* A [header file](util/TreeReaderTopLJets.h) for the `TreeReaderTopLJets` class
that provides an interface to read the variables from the input file easily.
This class is used to loop over all input events and load the event information.
* An [implementation file](util/example-top-ljets.cxx) that provides a
minimalistic implementation of KLFitter in a toy-analysis. The file contains
detailed comments about every command and substep to use KLFitter.


## Input file

The input file for the provided example (data/examples/top-ljets-input.root) is
a _HepSim sample_ taken [from the
homepage](http://atlaswww.hep.anl.gov/hepsim/info.php?item=142) of the HepSim
project [1]. The sample simulates ttbar+jet processes in proton-proton
collisions at a centre-of-mass energy of 13 TeV. The processes are generated
using the Madgraph matrix element generator interfaced with Herwig6. The
simulated detector corresponds to the [Snowmass
detector](https://arxiv.org/abs/1309.1057) which reflects the best performance
from future ATLAS and CMS detectors.

> [1] S.V.Chekanov. "HepSim: a repository with predictions for high-energy
> physics experiments", [arXiv:1403.1886](https://arxiv.org/abs/1403.1886)

The original Delphes format of the input ROOT files is transformed into a
standard ROOT tree, the events are skimmed and only variables relevant for the
KLFitter reconstruction are kept. The event skimming sets the following
requirements for events to be kept:
* at least 4 jets with a transverse momentum of at least 25 GeV. The eta of the
  jets is required to fulfil |eta| < 2.5.
* exactly one charged lepton (i.e. electron or muon) with a transverse momentum
  of at least 25 GeV and |eta| < 2.5.


### Input file variables

The input ROOT file contains only one ROOT tree, `nominal`. The branches in this
tree are described below:

* `lepton_pt`: Float. Measured transverse momentum of the charged lepton
  (electron or muon), in GeV.
* `lepton_eta`: Float. Measured eta of the charged lepton (electron or muon).
* `lepton_cl_eta`: Float. Measured eta of the charged lepton (electron or muon)
  as measured by the calorimeters ("cluster eta") - this is needed for electrons
  only.
* `lepton_phi`: Float. Measured phi of the charged lepton (electron or muon).
* `lepton_e`: Float. Measured energy of the charged lepton (electron or muon),
  in GeV.
* `met_met`: Float. Magnitude of the measured missing transverse momentum, in
  GeV.
* `met_phi`: Float. Phi component of the missing transverse momentum.
* `sumet`: Float. Scalar sum of the transverse energy of and event.
* `lepton_is_e`: Char. A flag that carries information if the lepton in an event
  is an electron.
* `lepton_is_mu`: Char. A flag that carries information if the lepton in an
  event is a muon. This in principle is redundant if `lepton_is_e`, but can be
  used for sanity checks.
* `jet_pt`: A vector a floats. Each element corresponds to the measured
  transverse momentum of a reconstructed jet in an event, in GeV.
* `jet_eta`: A vector a floats. Each element corresponds to the measured eta of
  a reconstructed jet in an event.
* `jet_phi`: A vector a floats. Each element corresponds to the measured phi of
  a reconstructed jet in an event.
* `jet_e`: A vector a floats. Each element corresponds to the measured energy of
  a reconstructed jet in an event, in GeV.
* `jet_btag_weight`: A vector a floats. Each element corresponds to the
  b-tagging weight of the corresponding reconstructed jet. A dummy values (1)
  are used in this example.
* `jet_has_btag`: A vector of chars. Each element represents a flag if the
  corresponding reconstructed jet is b-tagged or not.


## Output of the example

The example code prints some important variables for the first processed event.
The output variables are also stored in form of a ROOT file,
`top-ljets-output.root`, which contains the following branches:

* `klf_bhad_pt`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted pT of the b-quark assigned
  by KLFitter to come from top quark that decays hadronically.
* `klf_bhad_eta`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted eta of the b-quark assigned
  by KLFitter to come from the top quark that decays hadronically.
* `klf_bhad_phi`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted phi of the b-quark assigned
  by KLFitter to come from the top quark that decays hadronically.
* `klf_bhad_e`: A vector of floats. Each element corresponds to one permutation.
  Value in an element represents fitted energy of the b-quark assigned by
  KLFitter to come from the top quark that decays hadronically.
* `klf_bhad_jet_index`: A vector of unsigned integers. Each element corresponds
  to one permutation. Value in an element represents position in a jet vector of
  the b-quark assigned by KLFitter to come from the top quark that decays
  hadronically.
* `klf_blep_pt`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted pT of the b-quark assigned
  by KLFitter to come from the top quark that decays semi-leptonically.
* `klf_blep_eta`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted eta of b-quark assigned by
  KLFitter to come from the top quark that decays semi-leptonically.
* `klf_blep_phi`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted phi of the b-quark assigned
  by KLFitter to come from the top quark that decays semi-leptonically.
* `klf_blep_e`: A vector of floats. Each element corresponds to one permutation.
  Value in an element represents fitted energy of the b-quark assigned by
  KLFitter to come from the top quark that decays semi-leptonically.
* `klf_blep_jet_index`: A vector of unsigned integers. Each element corresponds
  to one permutation. Value in an element represents position in a jet vector of
  the b-quark assigned by KLFitter to come from the top quark that decays
  semi-leptonically.
* `klf_lquark1_pt`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted pT of the first quark
  assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark1_eta`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted eta of the first quark
  assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark1_phi`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted phi of the first quark
  assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark1_e`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted energy of the first quark
  assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark1_jet_index`: A vector of unsigned integers. Each element
  corresponds to one permutation. Value in an element represents position in a
  jet vector of the first quark assigned by KLFitter to come from the W boson
  that decays hadronically.
* `klf_lquark2_pt`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted pT of the second quark
  assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_eta`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted eta of the second quark
  assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_phi`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted phi of the second quark
  assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_e`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted energy of the second quark
  assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_jet_index`: A vector of unsigned integers. Each element
  corresponds to one permutation. Value in an element represents position in a
  jet vector of the second quark assigned by KLFitter to come from the W boson
  that decays hadronically.
* `klf_lepton_pt`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted pT of the lepton from the
  leptonically decaying W boson.
* `klf_lepton_eta`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted eta of the lepton from the
  leptonically decaying W boson.
* `klf_lepton_phi`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted phi of the lepton from the
  leptonically decaying W boson.
* `klf_lepton_e`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted energy of the lepton from
  the leptonically decaying W boson.
* `klf_neutrino_pt`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted pT of the neutrino from the
  leptonically decaying W boson.
* `klf_neutrino_eta`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted eta of the neutrino from
  the leptonically decaying W boson.
* `klf_neutrino_phi`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted phi of the neutrino from
  the leptonically decaying W boson.
* `klf_neutrino_e`: A vector of floats. Each element corresponds to one
  permutation. Value in an element represents fitted energy of the neutrino from
  the leptonically decaying W boson.
* `klf_loglikelihood`: A vector of doubles. Each element corresponds to one
  permutation. Value in an element represents logarithm of the KLFitter
  likelihood.
* `klf_event_probability`: A vector of doubles. Each element corresponds to one
  permutation. Value in an element represents event probability for given
  permutation.
* `klf_fit_minuit_did_not_converge`: A vector of chars. Each element corresponds
  to one permutation. Value in an element represents a flag about convergence of
  the fit.
* `klf_fit_aborted_to_nan`: A vector of chars. Each element corresponds to one
  permutation. Value in an element represents a flag storing information about
  fit being aborted due to NaN value.
* `klf_fit_parameter_at_limit`: A vector of chars. Each element corresponds to
  one permutation. Value in an element represents a flag storing information
  about fit being at the parameter(s) limits.
* `klf_fit_invalid_transfer_function`: A vector of chars. Each element
  corresponds to one permutation. Value in an element represents a flag storing
  information about fit being at the region where transfer functions are not
  valid.


## How to identify permutations with the highest event probability

If you want to use values corresponding to permutations with, for example, the
highest __event probability_, check the `klf_event_probability` vector and find
the position of the element with the highest value. Store the index of this
permutation in the vector. Then take the values from all variable vectors at
exactly the same position. These values correspond to the permutation with the
highest __event probability_.
