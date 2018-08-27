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

#ifndef KLFITTER_DETECTORATLAS_7TEV_H_
#define KLFITTER_DETECTORATLAS_7TEV_H_

#include <memory>
#include <string>

#include "KLFitter/DetectorBase.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
 * A class for describing of the ATLAS detector. This class
 * holds the description of the ATLAS detector.
 */
class DetectorAtlas_7TeV : public DetectorBase {
 public:
  /**
   * The default constructor.
   * @param folder The folder with transfer function parameters.
   */
  explicit DetectorAtlas_7TeV(std::string folder = "");

  /// The (defaulted) destructor.
  ~DetectorAtlas_7TeV();

  /** \name Member functions (Get)  */
  /** @{ */

  /**
   * Return the energy resolution of light jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  ResolutionBase* ResEnergyLightJet(double eta = 0.) override;

  /**
   * Return the energy resolution of b jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  ResolutionBase* ResEnergyBJet(double eta = 0.) override;

  /**
   * Return the energy resolution of gluon jets.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  ResolutionBase* ResEnergyGluonJet(double eta = 0.) override;

  /**
   * Return the energy resolution of electrons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  ResolutionBase* ResEnergyElectron(double eta = 0.) override;

  /**
   * Return the energy resolution of muons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  ResolutionBase* ResEnergyMuon(double eta = 0.) override;

  /**
   * Return the energy resolution of photons.
   * @param eta The eta of the particle.
   * @return A pointer to the energy resolution object.
   */
  ResolutionBase* ResEnergyPhoton(double eta = 0.) override;

  /**
   * Return the missing ET resolution.
   * @return A pointer to the missing ET resolution.
   */
  ResolutionBase* ResMissingET() override;

  /**
   * Return the eta resolution of light jets.
   * @param eta The eta of the particle.
   * @return A pointer to the eta resolution object.
   */
  ResolutionBase* ResEtaLightJet(double eta = 0.) override;

  /**
   * Return the eta resolution of b jets.
   * @param eta The eta of the particle.
   * @return A pointer to the eta resolution object.
   */
  ResolutionBase* ResEtaBJet(double eta = 0.) override;

  /**
   * Return the phi resolution of light jets.
   * @param eta The phi of the particle.
   * @return A pointer to the phi resolution object.
   */
  ResolutionBase* ResPhiLightJet(double eta = 0.) override;

  /**
   * Return the phi resolution of b jets.
   * @param eta The phi of the particle.
   * @return A pointer to the phi resolution object.
   */
  ResolutionBase* ResPhiBJet(double eta = 0.) override;

  /** @} */

 private:
  /// The energy resolution of light jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta1;
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta2;
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta3;
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta4;
  std::unique_ptr<ResolutionBase> m_res_energy_light_jet_eta5;

  /// The energy resolution of b jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_energy_bjet_eta1;
  std::unique_ptr<ResolutionBase> m_res_energy_bjet_eta2;
  std::unique_ptr<ResolutionBase> m_res_energy_bjet_eta3;
  std::unique_ptr<ResolutionBase> m_res_energy_bjet_eta4;
  std::unique_ptr<ResolutionBase> m_res_energy_bjet_eta5;

  /// The energy resolution of gluon jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_energy_gluon_jet_eta1;
  std::unique_ptr<ResolutionBase> m_res_energy_gluon_jet_eta2;
  std::unique_ptr<ResolutionBase> m_res_energy_gluon_jet_eta3;
  std::unique_ptr<ResolutionBase> m_res_energy_gluon_jet_eta4;

  /// The energy resolution of electrons for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_energy_electron_eta1;
  std::unique_ptr<ResolutionBase> m_res_energy_electron_eta2;
  std::unique_ptr<ResolutionBase> m_res_energy_electron_eta3;
  std::unique_ptr<ResolutionBase> m_res_energy_electron_eta4;

  /// The energy resolution of muons for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_energy_muon_eta1;
  std::unique_ptr<ResolutionBase> m_res_energy_muon_eta2;
  std::unique_ptr<ResolutionBase> m_res_energy_muon_eta3;

  /// The energy resolution of photons for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_energy_photon_eta1;
  std::unique_ptr<ResolutionBase> m_res_energy_photon_eta2;
  std::unique_ptr<ResolutionBase> m_res_energy_photon_eta3;
  std::unique_ptr<ResolutionBase> m_res_energy_photon_eta4;

  /// The eta resolution of light jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_eta_light_jet_eta1;
  std::unique_ptr<ResolutionBase> m_res_eta_light_jet_eta2;
  std::unique_ptr<ResolutionBase> m_res_eta_light_jet_eta3;
  std::unique_ptr<ResolutionBase> m_res_eta_light_jet_eta4;

  /// The eta resolution of b jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_eta_bjet_eta1;
  std::unique_ptr<ResolutionBase> m_res_eta_bjet_eta2;
  std::unique_ptr<ResolutionBase> m_res_eta_bjet_eta3;
  std::unique_ptr<ResolutionBase> m_res_eta_bjet_eta4;

  /// The phi resolution of light jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_phi_light_jet_eta1;
  std::unique_ptr<ResolutionBase> m_res_phi_light_jet_eta2;
  std::unique_ptr<ResolutionBase> m_res_phi_light_jet_eta3;
  std::unique_ptr<ResolutionBase> m_res_phi_light_jet_eta4;

  /// The phi resolution of b jets for different eta regions.
  std::unique_ptr<ResolutionBase> m_res_phi_bjet_eta1;
  std::unique_ptr<ResolutionBase> m_res_phi_bjet_eta2;
  std::unique_ptr<ResolutionBase> m_res_phi_bjet_eta3;
  std::unique_ptr<ResolutionBase> m_res_phi_bjet_eta4;

  /// Missing ET resolution in x and y
  std::unique_ptr<ResolutionBase> m_res_missing_ET;

  /// The eta binning for jets
  const double m_jet_eta_bin_1{0.8};
  const double m_jet_eta_bin_2{1.37};
  const double m_jet_eta_bin_3{1.52};
  const double m_jet_eta_bin_4{2.5};
  const double m_jet_eta_bin_5{4.5};

  /// The eta binning for electrons
  const double m_electron_eta_bin_1{0.8};
  const double m_electron_eta_bin_2{1.37};
  const double m_electron_eta_bin_3{1.52};
  const double m_electron_eta_bin_4{2.5};

  /// The eta binning for muons
  const double m_muon_eta_bin_1{1.11};
  const double m_muon_eta_bin_2{1.25};
  const double m_muon_eta_bin_3{2.5};

  /// The eta binning for muons
  const double m_photon_eta_bin_1{1.11};
  const double m_photon_eta_bin_2{1.25};
  const double m_photon_eta_bin_3{2.5};
  const double m_photon_eta_bin_4{3.0};
};
}  // namespace KLFitter

#endif  // KLFITTER_DETECTORATLAS_7TEV_H_
