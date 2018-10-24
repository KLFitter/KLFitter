/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */

// Include the interface for ROOT file reading
#include "TreeReaderTopLJets.h"

// c++ includes
#include <iostream>
#include <vector>

// KLFitter includes
#include "KLFitter/DetectorSnowmass.h"
#include "KLFitter/Fitter.h"
#include "KLFitter/LikelihoodTopLeptonJets.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"

// ROOT includes
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "ERROR: Expecting 1 argument but " << argc - 1;
    std::cerr << " arguments provided. Exiting." << std::endl;
    return 1;
  }

  std::string base_dir{argv[1]};

  // Open the input ROOT file and check it afterwards.
  auto f_input = TFile::Open((base_dir + "/data/examples/top-ljets-input.root").c_str(), "READ");
  if (f_input == nullptr) {
    std::cerr << "ERROR: Cannot open the input ROOT file. Aborting." << std::endl;
    return 1;
  }

  // Create and connect the instance of TTree for reading. Name
  // of the input tree is "nominal". Check afterwards if
  // everything went fine.
  auto t = static_cast<TTree*>(f_input->Get("nominal"));
  if (t == nullptr) {
    std::cerr << "ERROR: Cannot read 'nominal' tree from the input file. Aborting" << std::endl;
    return 1;
  }

  // Store the number of entries of the input file.
  const int nEntries = t->GetEntries();

  // Create an instance of the TreeReaderTopLJets that stores all
  // event information (i.e. all variables).
  TreeReaderTopLJets event(t);

  // Create an instance of a fitter.
  KLFitter::Fitter fitter{};

  // Create a detector - this is needed for the transfer
  // functions. If setting the detector fails, abort.
  KLFitter::DetectorSnowmass detector{(base_dir + "/data/transferfunctions/snowmass").c_str()};
  if (!fitter.SetDetector(&detector)) {
    std::cerr << "ERROR: Failed to set detector! Aborting" << std::endl;
    return 1;
  }

  // Create an instance of the lepton+jets likelihood.
  KLFitter::LikelihoodTopLeptonJets likelihood{};

  // Set the likelihood properties.
  likelihood.PhysicsConstants()->SetMassTop(172.5);  // mass in GeV
  // Other b-tagging mode options are e.g. 'kNotag' and 'kVeto'.
  // For all possible options, refer to the documentation.
  likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kWorkingPoint);
  likelihood.SetFlagTopMassFixed(true);

  // Set the likelihood in the fitter instance.
  if (!fitter.SetLikelihood(&likelihood)) {
    std::cerr << "ERROR: Failed to set likelihood. Aborting." << std::endl;
    return 1;
  }

  bool isFirst(true);

  // Open the output ROOT file.
  auto f_output = TFile::Open("top-ljets-output.root", "RECREATE");
  TTree out_tree{"KLFitter_output", "Example tree with KLFitter output"};

  // Declare variables to store event information.
  std::vector<float> klf_bhad_pt;
  std::vector<float> klf_bhad_eta;
  std::vector<float> klf_bhad_phi;
  std::vector<float> klf_bhad_e;
  std::vector<unsigned int> klf_bhad_jet_index;
  std::vector<float> klf_blep_pt;
  std::vector<float> klf_blep_eta;
  std::vector<float> klf_blep_phi;
  std::vector<float> klf_blep_e;
  std::vector<unsigned int> klf_blep_jet_index;
  std::vector<float> klf_lquark1_pt;
  std::vector<float> klf_lquark1_eta;
  std::vector<float> klf_lquark1_phi;
  std::vector<float> klf_lquark1_e;
  std::vector<unsigned int> klf_lquark1_jet_index;
  std::vector<float> klf_lquark2_pt;
  std::vector<float> klf_lquark2_eta;
  std::vector<float> klf_lquark2_phi;
  std::vector<float> klf_lquark2_e;
  std::vector<unsigned int> klf_lquark2_jet_index;
  std::vector<float> klf_lepton_pt;
  std::vector<float> klf_lepton_eta;
  std::vector<float> klf_lepton_phi;
  std::vector<float> klf_lepton_e;
  std::vector<float> klf_neutrino_pt;
  std::vector<float> klf_neutrino_eta;
  std::vector<float> klf_neutrino_phi;
  std::vector<float> klf_neutrino_e;
  std::vector<double> klf_loglikelihood;
  std::vector<double> klf_event_probability;
  std::vector<char> klf_fit_minuit_did_not_converge;
  std::vector<char> klf_fit_aborted_to_nan;
  std::vector<char> klf_fit_parameter_at_limit;
  std::vector<char> klf_fit_invalid_transfer_function;

  // Prepare the output variables as branches of the tree.
  out_tree.Branch("klf_bhad_pt", &klf_bhad_pt);
  out_tree.Branch("klf_bhad_eta", &klf_bhad_eta);
  out_tree.Branch("klf_bhad_phi", &klf_bhad_phi);
  out_tree.Branch("klf_bhad_e", &klf_bhad_e);
  out_tree.Branch("klf_bhad_jet_index", &klf_bhad_jet_index);
  out_tree.Branch("klf_blep_pt", &klf_blep_pt);
  out_tree.Branch("klf_blep_eta", &klf_blep_eta);
  out_tree.Branch("klf_blep_phi", &klf_blep_phi);
  out_tree.Branch("klf_blep_e", &klf_blep_e);
  out_tree.Branch("klf_blep_jet_index", &klf_blep_jet_index);
  out_tree.Branch("klf_lquark1_pt", &klf_lquark1_pt);
  out_tree.Branch("klf_lquark1_eta", &klf_lquark1_eta);
  out_tree.Branch("klf_lquark1_phi", &klf_lquark1_phi);
  out_tree.Branch("klf_lquark1_e", &klf_lquark1_e);
  out_tree.Branch("klf_lquark1_jet_index", &klf_lquark1_jet_index);
  out_tree.Branch("klf_lquark2_pt", &klf_lquark2_pt);
  out_tree.Branch("klf_lquark2_eta", &klf_lquark2_eta);
  out_tree.Branch("klf_lquark2_phi", &klf_lquark2_phi);
  out_tree.Branch("klf_lquark2_e", &klf_lquark2_e);
  out_tree.Branch("klf_lquark2_jet_index", &klf_lquark2_jet_index);
  out_tree.Branch("klf_lepton_pt", &klf_lepton_pt);
  out_tree.Branch("klf_lepton_eta", &klf_lepton_eta);
  out_tree.Branch("klf_lepton_phi", &klf_lepton_phi);
  out_tree.Branch("klf_lepton_e", &klf_lepton_e);
  out_tree.Branch("klf_neutrino_pt", &klf_neutrino_pt);
  out_tree.Branch("klf_neutrino_eta", &klf_neutrino_eta);
  out_tree.Branch("klf_neutrino_phi", &klf_neutrino_phi);
  out_tree.Branch("klf_neutrino_e", &klf_neutrino_e);
  out_tree.Branch("klf_loglikelihood", &klf_loglikelihood);
  out_tree.Branch("klf_event_probability", &klf_event_probability);
  out_tree.Branch("klf_fit_minuit_did_not_converge", &klf_fit_minuit_did_not_converge);
  out_tree.Branch("klf_fit_aborted_to_nan", &klf_fit_aborted_to_nan);
  out_tree.Branch("klf_fit_parameter_at_limit", &klf_fit_parameter_at_limit);
  out_tree.Branch("klf_fit_invalid_transfer_function", &klf_fit_invalid_transfer_function);

  // Loop over all events in the input tree.
  std::cout << "Started looping over " << nEntries << " entries" << std::endl;
  for (int ievent = 0; ievent < nEntries; ievent++) {
    if (ievent %100 == 0) {
      std::cout << "Processing event: " << ievent << " out of " << nEntries << " events" << std::endl;
    }

    // Read the information from the input ROOT file.
    event.GetEntry(ievent);

    // Clear all vectors for of the variables.
    klf_bhad_pt.clear();
    klf_bhad_eta.clear();
    klf_bhad_phi.clear();
    klf_bhad_e.clear();
    klf_bhad_jet_index.clear();
    klf_blep_pt.clear();
    klf_blep_eta.clear();
    klf_blep_phi.clear();
    klf_blep_e.clear();
    klf_blep_jet_index.clear();
    klf_lquark1_pt.clear();
    klf_lquark1_eta.clear();
    klf_lquark1_phi.clear();
    klf_lquark1_e.clear();
    klf_lquark1_jet_index.clear();
    klf_lquark2_pt.clear();
    klf_lquark2_eta.clear();
    klf_lquark2_phi.clear();
    klf_lquark2_e.clear();
    klf_lquark2_jet_index.clear();
    klf_lepton_pt.clear();
    klf_lepton_eta.clear();
    klf_lepton_phi.clear();
    klf_lepton_e.clear();
    klf_neutrino_pt.clear();
    klf_neutrino_eta.clear();
    klf_neutrino_phi.clear();
    klf_neutrino_e.clear();
    klf_loglikelihood.clear();
    klf_event_probability.clear();
    klf_fit_minuit_did_not_converge.clear();
    klf_fit_aborted_to_nan.clear();
    klf_fit_parameter_at_limit.clear();
    klf_fit_invalid_transfer_function.clear();

    // Create an instance of particles. You need to make sure that:
    // - the particles are in the range allowed by the transfer
    //   functions (eta and pt)
    // - the energies and momenta are in GeV
    //
    // Be aware that *all* particles you're adding are considered
    // in the fit (many particles lead to many permutations to be
    // considered and hence a long running time and not
    // necessarily good fitting results due to the many available
    // permutations)
    KLFitter::ParticleCollection particles{};

    // Add leptons. Depending on the two event variables
    // "lepton_is_e" and "lepton_is_mu", either an electron or a
    // lepton is added. Also set the lepton type as a parameter
    // of the likelihood.
    TLorentzVector lepton_p4;
    lepton_p4.SetPtEtaPhiE(event.lepton_pt, event.lepton_eta, event.lepton_phi, event.lepton_e);
    if (event.lepton_is_e) {
      likelihood.SetLeptonType(KLFitter::LikelihoodTopLeptonJets::kElectron);
      KLFitter::Particles::Electron el{"electron", lepton_p4};
      el.SetDetEta(event.lepton_cl_eta);
      particles.AddParticle(el);
    } else if (event.lepton_is_mu) {
      likelihood.SetLeptonType(KLFitter::LikelihoodTopLeptonJets::kMuon);
      KLFitter::Particles::Muon mu{"muon", lepton_p4};
      mu.SetDetEta(event.lepton_eta);
      particles.AddParticle(mu);
    } else {
      std::cerr << "WARNING: Event has no electrons or muons. Skipping." << std::endl;
      continue;
    }

    // Add jets - the input file already required at least 4 jets
    // per event, out of which at least 1 is b-tagged.
    for (unsigned int ijet = 0; ijet < 4; ijet++) {
      KLFitter::Particles::Parton parton{"parton" + std::to_string(ijet), TLorentzVector{}};
      parton.GetP4().SetPtEtaPhiE(event.jet_pt->at(ijet), event.jet_eta->at(ijet),
          event.jet_phi->at(ijet), event.jet_e->at(ijet));
      parton.SetDetEta(event.jet_eta->at(ijet));                // jet eta
      parton.SetIdentifier(ijet);                               // index of the jet to identify it
      parton.SetIsBTagged(event.jet_has_btag->at(ijet));        // Is the jet btagged?
      parton.SetBTagEfficiency(0.6);                            // tagging efficiency required for kWorkingPoint
      parton.SetBTagRejection(145.);                            // 1./tagging inefficiency required for kWorkingPoint
      parton.SetBTagWeight(event.jet_btag_weight->at(ijet));    // btag discriminant weight
      particles.AddParticle(parton);
    }

    // Add particles to the likelihood.
    if (!fitter.SetParticles(&particles)) {
      std::cerr << "ERROR: Failed to add particles to KLFitter. Aborting." << std::endl;
      return 1;
    }

    // Add MET information
    const double met_ex = event.met_met*cos(event.met_phi);
    const double met_ey = event.met_met*sin(event.met_phi);
    if (!fitter.SetET_miss_XY_SumET(met_ex, met_ey, event.sumet)) {
      std::cerr << "ERROR: Failed to add MET to fitter. Aborting." << std::endl;
      return 1;
    }

    // Loop over all permutations.
    const int nperm = fitter.Permutations()->NPermutations();
    for (int iperm  = 0; iperm < nperm; iperm++) {
      // Do the fitting magic.
      fitter.Fit(iperm);

      // Read the output and convergence status of the fit.
      unsigned int ConvergenceStatusBitWord = fitter.ConvergenceStatus();
      bool MinuitDidNotConverge = (ConvergenceStatusBitWord & fitter.MinuitDidNotConvergeMask) != 0;
      bool FitAbortedDueToNaN = (ConvergenceStatusBitWord & fitter.FitAbortedDueToNaNMask) != 0;
      bool AtLeastOneFitParameterAtItsLimit = (ConvergenceStatusBitWord & fitter.AtLeastOneFitParameterAtItsLimitMask) != 0;
      bool InvalidTransferFunctionAtConvergence = (ConvergenceStatusBitWord & fitter.InvalidTransferFunctionAtConvergenceMask) != 0;

      // Get log likelihood and event probability values. Note
      // that the event probablity is _not_ normalized.
      double likelihood = fitter.Likelihood()->LogLikelihood(fitter.Likelihood()->GetBestFitParameters());
      double event_probability = std::exp(fitter.Likelihood()->LogEventProbability());

      // Get the values of all fitted variables.
      auto modelParticles = fitter.Likelihood()->ParticlesModel();
      auto permutedParticles = fitter.Likelihood()->PParticlesPermuted();

      // Hadronic b quark.
      float bhad_pt = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 0)->Pt();
      float bhad_eta = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 0)->Eta();
      float bhad_phi = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 0)->Phi();
      float bhad_e = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 0)->E();
      unsigned int bhad_index = (*permutedParticles)->partons.at(0).GetIdentifier();

      // Leptonic b quark.
      float blep_pt = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 1)->Pt();
      float blep_eta = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 1)->Eta();
      float blep_phi = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 1)->Phi();
      float blep_e = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 1)->E();
      unsigned int blep_index = (*permutedParticles)->partons.at(1).GetIdentifier();

      // Light quark 1.
      float lquark1_pt = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 2)->Pt();
      float lquark1_eta = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 2)->Eta();
      float lquark1_phi = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 2)->Phi();
      float lquark1_e = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 2)->E();
      unsigned int lquark1_index = (*permutedParticles)->partons.at(2).GetIdentifier();

      // Light quark 2.
      float lquark2_pt = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 3)->Pt();
      float lquark2_eta = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 3)->Eta();
      float lquark2_phi = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 3)->Phi();
      float lquark2_e = modelParticles->GetP4(KLFitter::Particles::Type::kParton, 3)->E();
      unsigned int lquark2_index = (*permutedParticles)->partons.at(3).GetIdentifier();

      float lepton_pt = -9999;
      float lepton_eta = -9999;
      float lepton_phi = -9999;
      float lepton_e = -9999;

      // Always check for lepton type or the code will crash.
      if (event.lepton_is_e) {
        lepton_pt = modelParticles->GetP4(KLFitter::Particles::Type::kElectron, 0)->Pt();
        lepton_eta = modelParticles->GetP4(KLFitter::Particles::Type::kElectron, 0)->Eta();
        lepton_phi = modelParticles->GetP4(KLFitter::Particles::Type::kElectron, 0)->Phi();
        lepton_e = modelParticles->GetP4(KLFitter::Particles::Type::kElectron, 0)->E();
      } else if (event.lepton_is_mu) {
        lepton_pt = modelParticles->GetP4(KLFitter::Particles::Type::kMuon, 0)->Pt();
        lepton_eta = modelParticles->GetP4(KLFitter::Particles::Type::kMuon, 0)->Eta();
        lepton_phi = modelParticles->GetP4(KLFitter::Particles::Type::kMuon, 0)->Phi();
        lepton_e = modelParticles->GetP4(KLFitter::Particles::Type::kMuon, 0)->E();
      }

      // Neutrino parameters.
      float neutrino_pt = modelParticles->GetP4(KLFitter::Particles::Type::kNeutrino, 0)->Pt();
      float neutrino_eta = modelParticles->GetP4(KLFitter::Particles::Type::kNeutrino, 0)->Eta();
      float neutrino_phi = modelParticles->GetP4(KLFitter::Particles::Type::kNeutrino, 0)->Phi();
      float neutrino_e = modelParticles->GetP4(KLFitter::Particles::Type::kNeutrino, 0)->E();

      // Fill the vectors with the output variables. Note: it's
      // better/safer to store booleans as chars in ROOT files.
      klf_fit_minuit_did_not_converge.emplace_back(static_cast<char>(MinuitDidNotConverge));
      klf_fit_aborted_to_nan.emplace_back(static_cast<char>(FitAbortedDueToNaN));
      klf_fit_parameter_at_limit.emplace_back(static_cast<char>(AtLeastOneFitParameterAtItsLimit));
      klf_fit_invalid_transfer_function.emplace_back(static_cast<char>(InvalidTransferFunctionAtConvergence));

      klf_bhad_pt.emplace_back(bhad_pt);
      klf_bhad_eta.emplace_back(bhad_eta);
      klf_bhad_phi.emplace_back(bhad_phi);
      klf_bhad_e.emplace_back(bhad_e);
      klf_bhad_jet_index.emplace_back(bhad_index);
      klf_blep_pt.emplace_back(blep_pt);
      klf_blep_eta.emplace_back(blep_eta);
      klf_blep_phi.emplace_back(blep_phi);
      klf_blep_e.emplace_back(blep_e);
      klf_blep_jet_index.emplace_back(blep_index);
      klf_lquark1_pt.emplace_back(lquark1_pt);
      klf_lquark1_eta.emplace_back(lquark1_eta);
      klf_lquark1_phi.emplace_back(lquark1_phi);
      klf_lquark1_e.emplace_back(lquark1_e);
      klf_lquark1_jet_index.emplace_back(lquark1_index);
      klf_lquark2_pt.emplace_back(lquark2_pt);
      klf_lquark2_eta.emplace_back(lquark2_eta);
      klf_lquark2_phi.emplace_back(lquark2_phi);
      klf_lquark2_e.emplace_back(lquark2_e);
      klf_lquark2_jet_index.emplace_back(lquark2_index);
      klf_lepton_pt.emplace_back(lepton_pt);
      klf_lepton_eta.emplace_back(lepton_eta);
      klf_lepton_phi.emplace_back(lepton_phi);
      klf_lepton_e.emplace_back(lepton_e);
      klf_neutrino_pt.emplace_back(neutrino_pt);
      klf_neutrino_eta.emplace_back(neutrino_eta);
      klf_neutrino_phi.emplace_back(neutrino_phi);
      klf_neutrino_e.emplace_back(neutrino_e);
      klf_loglikelihood.emplace_back(likelihood);
      klf_event_probability.emplace_back(event_probability);

      // Print some values for the first event.
      if (isFirst) {
        printf("----------------------------------------------------------------------------------------------\n");
        printf("----------------------------------------Permutation %2i----------------------------------------\n", iperm);
        printf("----------------------------------------------------------------------------------------------\n");
        printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
        printf("Jet index         | %16i | %17i | %16i | %15i |\n",
        bhad_index, blep_index, lquark1_index, lquark2_index);
        printf("----------------------------------------------------------------------------------------------\n");
        printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
        printf("Output Energies   | %16.2f | %17.2f | %16.2f | %15.2f |\n",
        bhad_e, blep_e, lquark1_e, lquark2_e);
        printf("----------------------------------------------------------------------------------------------\n");
        printf("                  | lepton energy    | neutrino pz       | loglikelihood    |  probability    |\n");
        printf("Other values      | %16.2f | %17.2f | %16.2f | %15.2e |\n",
        lepton_e, neutrino_pt, likelihood, event_probability);
        printf("----------------------------------------------------------------------------------------------\n");
        printf("                  | Minuit Not Conv. | Fit Aborted: NaN  | >=1 Par at Limit | Invalid TF@Conv.|\n");
        printf("Status Code       | %16i | %17i | %16i | %15i |\n",
        MinuitDidNotConverge, FitAbortedDueToNaN, AtLeastOneFitParameterAtItsLimit, InvalidTransferFunctionAtConvergence);
      }
    }
    if (isFirst) isFirst = false;
    out_tree.Fill();
  }

  // Go to the output file and write it.
  f_output->cd();
  std::cout << std::endl << "Writing into output root file: " << "top-ljets-output.root" << std::endl << std::endl;
  out_tree.Write();

  // Close both input and output ROOT files.
  f_output->Close();
  f_input->Close();

  // Return 0 after a success run.
  return 0;
}
