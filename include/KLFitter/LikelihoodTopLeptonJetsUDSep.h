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

#ifndef KLFITTER_LIKELIHOODTOPLEPTONJETSUDSEP_H_
#define KLFITTER_LIKELIHOODTOPLEPTONJETSUDSEP_H_

#include <iostream>

#include "KLFitter/LikelihoodTopLeptonJets.h"

class TH1F;
class TH2F;
class TLorentzVector;

// ---------------------------------------------------------

namespace KLFitter {
class ResolutionBase;

/**
  * \class KLFitter::LikelihoodTopLeptonJetsUDSep
  * \brief A class implementing a likelihood for the ttbar lepton+jets channel.
  *
  * This class represents a likelihood for the ttbar into lepton+jets.
  */
class LikelihoodTopLeptonJetsUDSep : public LikelihoodTopLeptonJets {
 public:
  /// Enumerate for light-jet reweighting methods
  enum LJetSeparationMethod {
    kNone,           ///< description here
    kPermReweight,   ///< description here
    kPermReweight2D  ///< description here
  };

  /// \name Constructors and destructors
  /// @{

  /**
   * The default constructor.
   */
  LikelihoodTopLeptonJetsUDSep();

  /// The (defaulted) destructor.
  ~LikelihoodTopLeptonJetsUDSep();

  /// @}
  /// \name Member functions (Get)
  /// @{

  /// Returns the probability of a jet to have the tag weight and pT of an b type jet.
  double BJetProb(double tagweight, double pt);

  /// Returns the probability of a jet to have the pT of an b type jet.
  double BJetPt(double pt);

  /// Returns the probability of a jet to have the tag weight of an b type jet.
  double BJetTagWeight(double tagweight);

  /// Returns the probability of a jet to have the tag weight and pT of an down type jet.
  double DownJetProb(double tagweight, double pt);

  /// Returns the probability of a jet to have the pT of an down type jet.
  double DownJetPt(double pt);

  /// Returns the probability of a jet to have the tag weight of an down type jet.
  double DownJetTagWeight(double tagweight);

  /// Returns the probability of a jet to have the tag weight and pT of an up type jet.
  double UpJetProb(double tagweight, double pt);

  /// Returns the probability of a jet to have the pT of an up type jet.
  double UpJetPt(double pt);

  /// Returns the probability of a jet to have the tag weight of an up type jet.
  double UpJetTagWeight(double tagweight);

  /// @}
  /// \name Member functions (Set)
  /// @{

  /// Set a flag. If flag is true the permutations are reweighted with the pT and tag weight probabilities.
  void SetLJetSeparationMethod(LJetSeparationMethod flag) { m_ljet_separation_method = flag; }

  /**
   * Set histogram for tag weight distribution of b jets.
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetBJet2DWeightHisto(TH2F* hist) { m_bjet_2d_weight_histo = hist; return 1; }

  /**
   * Set histogram for pT distribution of b jets (reco level).
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetBJetPtHisto(TH1F* hist) { m_bjet_pt_histo = hist; return 1; }

  /**
   * Set histogram for tag weight distribution of b jets.
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetBJetTagWeightHisto(TH1F* hist) { m_bjet_tag_weight_histo = hist; return 1; }

  /**
   * Set histogram for tag weight distribution of down type jets.
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetDownJet2DWeightHisto(TH2F* hist) { m_down_jet_2d_weight_histo = hist; return 1; }

  /**
   * Set histogram for pT distribution of down jets (reco level).
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetDownJetPtHisto(TH1F* hist) { m_down_jet_pt_histo = hist; return 1; }

  /**
   * Set histogram for tag weight distribution of down type jets.
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetDownJetTagWeightHisto(TH1F* hist) { m_down_jet_tag_weight_histo = hist; return 1; }

  /**
   * Set histogram for tag weight distribution of up type jets.
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetUpJet2DWeightHisto(TH2F* hist) { m_up_jet_2d_weight_histo = hist; return 1; }

  /**
   * Set histogram for pT distribution of up jets (reco level).
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetUpJetPtHisto(TH1F* hist) { m_up_jet_pt_histo = hist; return 1; }

  /**
   * Set histogram for tag weight distribution of up type jets.
   * @param hist Pointer to histogram.
   * @return An error flag.
   */
  int SetUpJetTagWeightHisto(TH1F* hist) { m_up_jet_tag_weight_histo = hist; return 1; }

  /// @}
  /// \name Member functions (BAT)
  /// @{

  /**
   * Define the parameters of the fit.
   */
  void DefineParameters() override;

  /**
   * Return the log of the event probability fof the current
   * combination
   * @return The event probability
   */
  double LogEventProbability() override;

  /**
   * Return the contribution from b tagging to the log of the
   * event probability for the current combination
   * @return The event probability contribution
   */
  double LogEventProbabilityBTag() override;

  /// @}
  /// \name Member functions (misc)
  /// @{

  /**
   * Check if the permutation is LH invariant.
   * @return Permutation of the invariant partner, -1 if there is no one.
   */
  int LHInvariantPermutationPartner(int iperm, int nperms, int *switchpar1, int *switchpar2) override;

  /**
   * Return the contribution from pT and b tag weight probability (by LJetSeparationMethod)
   * to the log of the event probability for the current combination
   * @return The event probability contribution
   */
  double LogEventProbabilityLJetReweight();

  /// @}

 protected:
  /**
   * Define the model particles
   * @return An error code.
   */
  int DefineModelParticles() override;

  /**
   * Remove forbidden particle permutations.
   * @return An error code.
   */
  int RemoveForbiddenParticlePermutations() override { return 1; }

  /**
   * Remove invariant particle permutations.
   * @return An error code.
   */
  int RemoveInvariantParticlePermutations() override;

  /// @{
  /// \name Member attributes

  /// A flag for using an additional reweighting of the permutations with the pT
  /// and tag weight probability (default: false);
  LJetSeparationMethod m_ljet_separation_method;

  /// A pointer to the 2d histogram "tag weight vs. pT for bQuarks"
  TH2F* m_bjet_2d_weight_histo;

  /// A pointer to the histogram of the down b jet pT distribution.
  TH1F* m_bjet_pt_histo;

  /// A pointer to the histogram of the up b tag weight distribution.
  TH1F* m_bjet_tag_weight_histo;

  /// A pointer to the 2d histogram "tag weight vs. pT for downQuarks"
  TH2F* m_down_jet_2d_weight_histo;

  /// A pointer to the histogram of the down type jet pT distribution.
  TH1F* m_down_jet_pt_histo;

  /// A pointer to the histogram of the down quark tag weight distribution.
  TH1F* m_down_jet_tag_weight_histo;

  /// A pointer to the 2d histogram "tag weight vs. pT for upQuarks"
  TH2F* m_up_jet_2d_weight_histo;

  /// A pointer to the histogram of the up type jet pT distribution.
  TH1F* m_up_jet_pt_histo;

  /// A pointer to the histogram of the up quark tag weight distribution.
  TH1F* m_up_jet_tag_weight_histo;

  /// @}
};
}  // namespace KLFitter

#endif  // KLFITTER_LIKELIHOODTOPLEPTONJETSUDSEP_H_
