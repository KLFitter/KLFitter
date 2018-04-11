# Example of a KLFitter implementation

### Introduction
The purpose of the provided example is to demonstrate how to implement KLFitter in your analysis. The example code and the input files are designed to provide minimalistic input and the implementation. The example also shows how to retrieve the output from the KLFitter. Lepton+jets likelihood has been chosen in the example, but from the implementation it should be straightforward to extend the implementation to other available likelihoods described in the KLFitter [documentation](doc/WhatIsKLF.md).

### How to run the example
To run the example, compile the code then run the executable and provide the path to the KLFitter working directory.
E.g. `./example-top-ljets.exe ./`

### Implementation
The implementation of the example utilises three different files: An input ROOT file with the necessary information needed to run KLFitter lepton+jets likelihood (described in detail [here](doc/WhatIsKLF.md)), a header [file](util/TreeReaderTopLJets.h) that provides an interface between the input ROOT file and KLFitter implementation, i.e. allows to read the variables from the ROOT file and finally a source [file](util/example-top-ljets.cxx) that provides a minimalistic implementation of the KLFitter in a toy-analysis. The implementation file includes detailed comments about every step of the KLFitter implementation. 

### Input file
The input file for the provided example (data/examples/top-ljets-input.root) is a [HepSim sample](http://atlaswww.hep.anl.gov/hepsim/info.php?item=142). The samples simulates ttbar+jet processes in proton-proton collisions at center-of-mass energy 13 TeV. The processes are generated using Madgraph matrix element generator interfaced with Herwig6. Simulated detector corresponds to [snowmass detector](https://arxiv.org/abs/1309.1057) which reflects the best performance from future ATLAS and CMS detectors.

The original delphes format of the input ROOT files is skimmed and only variables relevant to KLFitter are kept. Furthermore, additional event selection is applied:
* at least 4 jets with pT > 25 GeV and |eta| < 2.5 are required
* exactly one charged lepton with pT > 25 GeV and |eta| < 2.5 is required

#### Input file variables
The input ROOT file contains only one ROOT tree, `nominal`. The branches in this tree are described below:

* `lepton_pt`: Float. Measured transverse momentum of the charged lepton (electron or muon), in GeV.
* `lepton_eta`: Float. Measured eta of the charged lepton (electron or muon).
* `lepton_cl_eta`: Float. Measured eta of the charged lepton (electron or muon) as measured by the calorimeters ("cluster eta") - this is needed for electrons only.
* `lepton_phi`: Float. Measured phi of the charged lepton (electron or muon).
* `lepton_e`: Float. Measured energy of the charged lepton (electron or muon), in GeV.
* `met_met`: Float. Magnitude of the measured missing transverse momentum, in GeV.
* `met_phi`: Float. Phi component of the missing transverse momentum.
* `sumet`: Float. Scalar sum of the transverse energy of and event.
* `lepton_is_e`: Char. A flag that carries information if the lepton in an event is an electron.
* `lepton_is_mu`: Char. A flag that carries information if the lepton in an event is a muon. This in principle is redundant if `lepton_is_e`, but can be used for sanity checks.
* `jet_pt`: A vector a floats. Each element corresponds to the measured transverse momentum of a reconstructed jet in an event, in GeV.
* `jet_eta`: A vector a floats. Each element corresponds to the measured eta of a reconstructed jet in an event.
* `jet_phi`: A vector a floats. Each element corresponds to the measured phi of a reconstructed jet in an event.
* `jet_e`: A vector a floats. Each element corresponds to the measured energy of a reconstructed jet in an event, in GeV.
* `jet_btag_weight`: A vector a floats. Each element corresponds to the b-tagging weight of the corresponding reconstructed jet. A dummy values (1) are used in this example.
* `jet_has_btag`: A vector of chars. Each element represents a flag if the corresponding reconstructed jet is b-tagged or not.


### Output of the example
The example code prints some important variables for the first processed event.
The output in the provided examples comes in form of a ROOT file `ljets_output.root` which
contains the following branches:

* `klf_bhad_pt`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted pT of the b-quark assigned by KLFitter to come from top quark that decays hadronically.
* `klf_bhad_eta`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted eta of the b-quark assigned by KLFitter to come from the top quark that decays hadronically.
* `klf_bhad_phi`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted phi of the b-quark assigned by KLFitter to come from the top quark that decays hadronically.
* `klf_bhad_e`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted energy of the b-quark assigned by KLFitter to come from the top quark that decays hadronically.
* `klf_bhad_jet_index`: A vector of unsigned integers. Each element corresponds to one permutation. Value in an element represents position in a jet vector of the b-quark assigned by KLFitter to come from the top quark that decays hadronically.
* `klf_blep_pt`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted pT of the b-quark assigned by KLFitter to come from the top quark that decays semi-leptonically.
* `klf_blep_eta`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted eta of b-quark assigned by KLFitter to come from the top quark that decays semi-leptonically.
* `klf_blep_phi`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted phi of the b-quark assigned by KLFitter to come from the top quark that decays semi-leptonically.
* `klf_blep_e`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted energy of the b-quark assigned by KLFitter to come from the top quark that decays semi-leptonically.
* `klf_blep_jet_index`: A vector of unsigned integers. Each element corresponds to one permutation. Value in an element represents position in a jet vector of the b-quark assigned by KLFitter to come from the top quark that decays semi-leptonically.
* `klf_lquark1_pt`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted pT of the first quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark1_eta`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted eta of the first quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark1_phi`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted phi of the first quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark1_e`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted energy of the first quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark1_jet_index`: A vector of unsigned integers. Each element corresponds to one permutation. Value in an element represents position in a jet vector of the first quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_pt`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted pT of the second quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_eta`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted eta of the second quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_phi`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted phi of the second quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_e`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted energy of the second quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lquark2_jet_index`: A vector of unsigned integers. Each element corresponds to one permutation. Value in an element represents position in a jet vector of the second quark assigned by KLFitter to come from the W boson that decays hadronically.
* `klf_lepton_pt`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted pT of the lepton from the leptonically decaying W boson.
* `klf_lepton_eta`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted eta of the lepton from the leptonically decaying W boson.
* `klf_lepton_phi`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted phi of the lepton from the leptonically decaying W boson.
* `klf_lepton_e`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted energy of the lepton from the leptonically decaying W boson.
* `klf_neutrino_pt`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted pT of the neutrino from the leptonically decaying W boson.
* `klf_neutrino_eta`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted eta of the neutrino from the leptonically decaying W boson.
* `klf_neutrino_phi`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted phi of the neutrino from the leptonically decaying W boson.
* `klf_neutrino_e`: A vector of floats. Each element corresponds to one permutation. Value in an element represents fitted energy of the neutrino from the leptonically decaying W boson.
* `klf_loglikelihood`: A vector of doubles. Each element corresponds to one permutation. Value in an element represents logarithm of the KLFitter likelihood.
* `klf_event_probability`: A vector of doubles. Each element corresponds to one permutation. Value in an element represents event probability for given permutation.
* `klf_fit_minuit_did_not_converge`: A vector of chars. Each element corresponds to one permutation. Value in an element represents a flag about convergence of the fit.
* `klf_fit_aborted_to_nan`: A vector of chars. Each element corresponds to one permutation. Value in an element represents a flag storing information about fit being aborted due to NaN value.
* `klf_fit_parameter_at_limit`: A vector of chars. Each element corresponds to one permutation. Value in an element represents a flag storing information about fit being at the parameter(s) limits.
* `klf_fit_invalid_transfer_function`: A vector of chars. Each element corresponds to one permutation. Value in an element represents a flag storing information about fit being at the region where transfer functions are not valid.

### How to identify permutations with the highest Event probability
If you want to e.g. use values corresponding to permutations with the highest `Event probability` check the `klf_event_probability` vector and find the position of the element with the highest value. Then take the values from each vector at exactly the same position. These values correspond to the permutation with the highest `Event probability`.
