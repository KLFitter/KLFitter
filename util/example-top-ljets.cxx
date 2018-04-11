// Include the interface for ROOT file reading
#include "TreeReaderTopLJets.h"

// KLFitter includes
#include "KLFitter/Fitter.h"
#include "KLFitter/DetectorAtlas_8TeV.h"
#include "KLFitter/LikelihoodTopLeptonJets.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/Permutations.h"

// ROOT includes
#include "TFile.h"

// c++ includes
#include <iostream>
#include <memory>

int main(int argc, char *argv[]){
  if (argc != 2){
    std::cerr << "ERROR: Expecting 1 argument but " << argc-1 << " arguments provided. Exiting." << std::endl;
    return 1;
  }

  std::string base_dir = argv[1];

  // Open the input ROOT file
  auto f_input = TFile::Open((base_dir + "/data/examples/top-ljets-input.root").c_str(), "READ");

  //Check if we opened the correct file
  if (f_input == nullptr){
    std::cerr << "ERROR: Cannot open the input ROOT file. Aborting." << std::endl;
    return 1;
  }

  // Create and connect the instance of TTree for reading
  // Name of the input tree is "nominal"
  TTree *t = (TTree*)f_input->Get("nominal");

  // Check if everything went fine
  if (t == nullptr){
    std::cerr << "ERROR: Cannot read 'nominal' tree from the input file. Aborting" << std::endl;
    return 1;
  }

  // Store the number of entried in the input file
  const int nEntries = t->GetEntries();

  // Create the instance of object that stores the event information
  TreeReaderTopLJets event(t);

  // Create new fitter
  KLFitter::Fitter fitter{};

  // Create detector - this is needed for transfer functions
  // The energy in the name of the class is irrelevant
  // it is just a dummy name
  KLFitter::DetectorAtlas_8TeV detector{(base_dir + "/data/transferfunctions/8TeV/ttbar/mc12_LCJets_v1").c_str()}; // 8 TeV TF

  if (!fitter.SetDetector(&detector)){
    std::cerr << "ERROR: Failed to set detector! Aborting" << std::endl;
    return 1;
  }

  // Create instance of likelihood
  KLFitter::LikelihoodTopLeptonJets likelihood{};

  // Set likelihood properties
  likelihood.PhysicsConstants()->SetMassTop(172.5); // in GeV
  likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kWorkingPoint); // kNotag, kVeto, kVetoLight, kVetoBoth, kWorkingPoint
  likelihood.SetFlagTopMassFixed(true);
  likelihood.SetLeptonType(KLFitter::LikelihoodTopLeptonJets::LeptonType::kMuon);

  // Set the likelihood
  if (!fitter.SetLikelihood(&likelihood)){
    std::cerr << "ERROR: Failed to set likelihood. Aborting." << std::endl;
    return 1;
  }

  bool isFirst(true);

  // open output ROOT file
  auto f_output = TFile::Open("top-ljets-output.root","RECREATE");
  TTree *out_tree = new TTree("KLFitter_output", "Example tree with KLFitter output");

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

  // prepare the output variables in a tree
  out_tree->Branch("klf_bhad_pt",&klf_bhad_pt);
  out_tree->Branch("klf_bhad_eta",&klf_bhad_eta);
  out_tree->Branch("klf_bhad_phi",&klf_bhad_phi);
  out_tree->Branch("klf_bhad_e",&klf_bhad_e);
  out_tree->Branch("klf_bhad_jet_index",&klf_bhad_jet_index);
  out_tree->Branch("klf_blep_pt",&klf_blep_pt);
  out_tree->Branch("klf_blep_eta",&klf_blep_eta);
  out_tree->Branch("klf_blep_phi",&klf_blep_phi);
  out_tree->Branch("klf_blep_e",&klf_blep_e);
  out_tree->Branch("klf_blep_jet_index",&klf_blep_jet_index);
  out_tree->Branch("klf_lquark1_pt",&klf_lquark1_pt);
  out_tree->Branch("klf_lquark1_eta",&klf_lquark1_eta);
  out_tree->Branch("klf_lquark1_phi",&klf_lquark1_phi);
  out_tree->Branch("klf_lquark1_e",&klf_lquark1_e);
  out_tree->Branch("klf_lquark1_jet_index",&klf_lquark1_jet_index);
  out_tree->Branch("klf_lquark2_pt",&klf_lquark2_pt);
  out_tree->Branch("klf_lquark2_eta",&klf_lquark2_eta);
  out_tree->Branch("klf_lquark2_phi",&klf_lquark2_phi);
  out_tree->Branch("klf_lquark2_e",&klf_lquark2_e);
  out_tree->Branch("klf_lquark2_jet_index",&klf_lquark2_jet_index);
  out_tree->Branch("klf_lepton_pt",&klf_lepton_pt);
  out_tree->Branch("klf_lepton_eta",&klf_lepton_eta);
  out_tree->Branch("klf_lepton_phi",&klf_lepton_phi);
  out_tree->Branch("klf_lepton_e",&klf_lepton_e);
  out_tree->Branch("klf_neutrino_pt",&klf_neutrino_pt);
  out_tree->Branch("klf_neutrino_eta",&klf_neutrino_eta);
  out_tree->Branch("klf_neutrino_phi",&klf_neutrino_phi);
  out_tree->Branch("klf_neutrino_e",&klf_neutrino_e);
  out_tree->Branch("klf_loglikelihood",&klf_loglikelihood);
  out_tree->Branch("klf_event_probability",&klf_event_probability);
  out_tree->Branch("klf_fit_minuit_did_not_converge",&klf_fit_minuit_did_not_converge);
  out_tree->Branch("klf_fit_aborted_to_nan",&klf_fit_aborted_to_nan);
  out_tree->Branch("klf_fit_parameter_at_limit",&klf_fit_parameter_at_limit);
  out_tree->Branch("klf_fit_invalid_transfer_function",&klf_fit_invalid_transfer_function);

  // loop over events
  std::cout << "Started looping over " << nEntries << " entries" << std::endl;
  for (int ievent = 0; ievent < nEntries; ievent++){
    if (ievent %100 == 0){
      std::cout << "Processing event: " << ievent << " out of " << nEntries << " events" << std::endl;
    }
    //Read the information from the input ROOT file
    event.GetEntry(ievent);

    // clear all vectors for output file
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

    // Create instance of particles

    // here, you need to make sure that
      // - the particles are in the range allowed by the transfer functions (eta and pt)
      // - the energies and momenta are in GeV
      // - be aware that *all* particles you're adding are considered in the fit
      //   (many particles lead to many permutations to be considered and hence a long
      //   running time and not necessarily good fitting results due to the many available
      //   permutations)
      // the arguments taken py AddParticle() are
      // - TLorentzVector of the physics 4-momentum
      // - detector eta for the evaluation of the transfer functions (for muons: just use the physics eta)
      // - type of particle
      // - an optional name of the particle (pass empty string in case you don't want to give your particle a name)
      // - index of the particle in your original collection (for convenience)
    KLFitter::Particles particles{};

    // Add leptons
    TLorentzVector lepton;
    lepton.SetPtEtaPhiE(event.lepton_pt, event.lepton_eta, event.lepton_phi, event.lepton_e);
    if (event.lepton_is_e){
        likelihood.SetLeptonType(1); // set lepton type to electron
      particles.AddParticle(&lepton, event.lepton_eta, KLFitter::Particles::kElectron);
    } else if (event.lepton_is_mu){
        likelihood.SetLeptonType(2); // set lepton type to muon
      particles.AddParticle(&lepton, event.lepton_eta, KLFitter::Particles::kMuon);
    } else {
      std::cerr << "WARNING: Event has no electrons nor muons. Skipping." << std::endl;
      continue;
    }

    // Add jets - the input file required at least 4 jets per event and at least 1 btagged
    for (unsigned int ijet = 0; ijet < 4; ijet++){ // we only want to consider 4 jets
      TLorentzVector jet;
      jet.SetPtEtaPhiE(event.jet_pt->at(ijet), event.jet_eta->at(ijet),
          event.jet_phi->at(ijet), event.jet_e->at(ijet));
      // arguments are as follows:
      //  1) TLorentzVector of jet
      //  2) jet eta
      //  3) KLFitter particle type. kParton for jets
      //  4) Internal name
      //  5) Index of the jet
      //  6) Is the jet btagged?
      //  7) tagging effciency required for kWorkingPoint
      //  8) 1./tagging inefficiency required for kWorkingPoint
      //  9) true flavour type
      //  10) btag weight
      particles.AddParticle(&jet, event.jet_eta->at(ijet), KLFitter::Particles::kParton, "", ijet,
          (int)event.jet_has_btag->at(ijet), 0.6, 145., KLFitter::Particles::kNone, event.jet_btag_weight->at(ijet));
    }

    // Add particles to the likelihood

    //fitter.SetParticles(p.get());
    if (!fitter.SetParticles(&particles)){
      std::cerr << "ERROR: Failed to add particles to KLFitter. Aborting." << std::endl;
      return 1;
    }
    // Add MET information
    const double met_ex = event.met_met*cos(event.met_phi);
    const double met_ey = event.met_met*sin(event.met_phi);


    if (!fitter.SetET_miss_XY_SumET(met_ex, met_ey, event.sumet)){
      std::cerr << "ERROR: Failed to add MET to fitter. Aborting." << std::endl;
      return 1;
    }

    // Loop over all permutations
    const int nperm = fitter.Permutations()->NPermutations();
    for (int iperm  = 0; iperm < nperm; iperm++){
      // Do the fitting magic
      fitter.Fit(iperm);

      // Read the output of the fit
      //
      // Minuit converge status
      unsigned int ConvergenceStatusBitWord = fitter.ConvergenceStatus();
      bool MinuitDidNotConverge = (ConvergenceStatusBitWord &
        fitter.MinuitDidNotConvergeMask) != 0;
      bool FitAbortedDueToNaN = (ConvergenceStatusBitWord &
        fitter.FitAbortedDueToNaNMask) != 0;
      bool AtLeastOneFitParameterAtItsLimit = (ConvergenceStatusBitWord &
        fitter.AtLeastOneFitParameterAtItsLimitMask) != 0;
      bool InvalidTransferFunctionAtConvergence = (ConvergenceStatusBitWord &
        fitter.InvalidTransferFunctionAtConvergenceMask) != 0;

      //Get likelihood and event probability values
      double likelihood = fitter.Likelihood()->LogLikelihood(
        fitter.Likelihood()->GetBestFitParameters());

      // probbaility is NOT normalized!
      double event_probability = std::exp(fitter.Likelihood()->LogEventProbability());

      // Get the values of fitted variables
      auto modelParticles = fitter.Likelihood()->ParticlesModel();
      auto permutedParticles = fitter.Likelihood()->PParticlesPermuted();

      // bhad quark
      float bhad_pt = modelParticles->Parton(0)->Pt();
      float bhad_eta = modelParticles->Parton(0)->Eta();
      float bhad_phi = modelParticles->Parton(0)->Phi();
      float bhad_e = modelParticles->Parton(0)->E();
      unsigned int bhad_index = (*permutedParticles)->JetIndex(0);

      // blep quark
      float blep_pt = modelParticles->Parton(1)->Pt();
      float blep_eta = modelParticles->Parton(1)->Eta();
      float blep_phi = modelParticles->Parton(1)->Phi();
      float blep_e = modelParticles->Parton(1)->E();
      unsigned int blep_index = (*permutedParticles)->JetIndex(1);

      //light quark 1
      float lquark1_pt = modelParticles->Parton(2)->Pt();
      float lquark1_eta = modelParticles->Parton(2)->Eta();
      float lquark1_phi = modelParticles->Parton(2)->Phi();
      float lquark1_e = modelParticles->Parton(2)->E();
      unsigned int lquark1_index = (*permutedParticles)->JetIndex(2);

      //light quark 2
      float lquark2_pt = modelParticles->Parton(3)->Pt();
      float lquark2_eta = modelParticles->Parton(3)->Eta();
      float lquark2_phi = modelParticles->Parton(3)->Phi();
      float lquark2_e = modelParticles->Parton(3)->E();
      unsigned int lquark2_index = (*permutedParticles)->JetIndex(3);

      float lepton_pt = -9999;
      float lepton_eta = -9999;
      float lepton_phi = -9999;
      float lepton_e = -9999;

      // always check for lepton type or the code will crash
      if (event.lepton_is_e){
        // electron
        lepton_pt = modelParticles->Electron(0)->Pt();
        lepton_eta = modelParticles->Electron(0)->Eta();
        lepton_phi = modelParticles->Electron(0)->Phi();
        lepton_e = modelParticles->Electron(0)->E();
      } else if (event.lepton_is_mu){
        // muon
        lepton_pt = modelParticles->Muon(0)->Pt();
        lepton_eta = modelParticles->Muon(0)->Eta();
        lepton_phi = modelParticles->Muon(0)->Phi();
        lepton_e = modelParticles->Muon(0)->E();
      }

      // neutrino
      float neutrino_pt = modelParticles->Neutrino(0)->Pt();
      float neutrino_eta = modelParticles->Neutrino(0)->Eta();
      float neutrino_phi = modelParticles->Neutrino(0)->Phi();
      float neutrino_e = modelParticles->Neutrino(0)->E();

      // fill the vectors
      // it's better/safer to store booleans as chars in ROOT files
      klf_fit_minuit_did_not_converge.emplace_back(static_cast<char> (MinuitDidNotConverge) );
      klf_fit_aborted_to_nan.emplace_back(static_cast<char> (FitAbortedDueToNaN) );
      klf_fit_parameter_at_limit.emplace_back(static_cast<char> (AtLeastOneFitParameterAtItsLimit) );
      klf_fit_invalid_transfer_function.emplace_back(static_cast<char> (InvalidTransferFunctionAtConvergence) );

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

      // print some values for the first event
      if (isFirst){
        printf("----------------------------------------------------------------------------------------------\n");
        printf("----------------------------------------Permutation %2i----------------------------------------\n",iperm);
        printf("----------------------------------------------------------------------------------------------\n");
        printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
        printf("Jet index         | %16i | %17i | %16i | %15i |\n",
        bhad_index, blep_index, lquark1_index, lquark2_index );
        printf("----------------------------------------------------------------------------------------------\n");
        printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
        printf("Output Energies   | %16.2f | %17.2f | %16.2f | %15.2f |\n",
        bhad_e, blep_e, lquark1_e, lquark2_e );
        printf("----------------------------------------------------------------------------------------------\n");
        printf("                  | lepton energy    | neutrino pz       | loglikelihood    |  probability    |\n");
        printf("Other values      | %16.2f | %17.2f | %16.2f | %15.2e |\n",
        lepton_e, neutrino_pt, likelihood, event_probability );
        printf("----------------------------------------------------------------------------------------------\n");
        printf("                  | Minuit Not Conv. | Fit Aborted: NaN  | >=1 Par at Limit | Invalid TF@Conv.|\n");
        printf("Status Code       | %16i | %17i | %16i | %15i |\n",
        MinuitDidNotConverge, FitAbortedDueToNaN, AtLeastOneFitParameterAtItsLimit, InvalidTransferFunctionAtConvergence);
      }
    }
    if (isFirst) isFirst = false;
    out_tree->Fill();
  }

  // Close the ROOT files

  // go to output file
  f_output->cd();
  // write to output file
  std::cout << std::endl << "Writing into output root file: " << "ljets_output.root" << std::endl << std::endl;
  out_tree->Write();

  // close ROOT files
  f_output->Close();
  f_input->Close();

  // Successful run
  return 0;
}
