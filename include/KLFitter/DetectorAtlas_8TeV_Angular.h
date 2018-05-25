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

#ifndef KLFITTER_DETECTORATLAS_8TEV_ANGULAR_H_
#define KLFITTER_DETECTORATLAS_8TEV_ANGULAR_H_

#include <memory>
#include <string>

#include "KLFitter/DetectorAtlas_8TeV.h"

// ---------------------------------------------------------

namespace KLFitter {
/**
 * \brief Extension of DetectorAtlas_8TeV with angular transfer
 * functions.

 * This version of the detector adds additional functionality for
 * eta and phi resolutions.
 */
class DetectorAtlas_8TeV_Angular : public DetectorAtlas_8TeV {
 public:
  /// \name Constructors and destructors
  /// @{

  /**
   * The default constructor, which extends the constructor of
   * DetectorAtlas_8TeV() by eta and phi resolutions.
   * @param folder The folder with transfer function parameters.
   */
  explicit DetectorAtlas_8TeV_Angular(std::string folder = "");

  /// The (defaulted) destructor.
  ~DetectorAtlas_8TeV_Angular();

  /// @}
  /// \name Member functions (Get)
  /// @{

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

  /// @}

 private:
  /// \name Eta resolutions for light and b jets.
  /// @{

  std::unique_ptr<ResolutionBase> fResEtaLightJet_eta1;
  std::unique_ptr<ResolutionBase> fResEtaLightJet_eta2;
  std::unique_ptr<ResolutionBase> fResEtaLightJet_eta3;
  std::unique_ptr<ResolutionBase> fResEtaLightJet_eta4;

  std::unique_ptr<ResolutionBase> fResEtaBJet_eta1;
  std::unique_ptr<ResolutionBase> fResEtaBJet_eta2;
  std::unique_ptr<ResolutionBase> fResEtaBJet_eta3;
  std::unique_ptr<ResolutionBase> fResEtaBJet_eta4;

  /// @}
  /// \name Phi resolutions for light and b jets.
  /// @{

  std::unique_ptr<ResolutionBase> fResPhiLightJet_eta1;
  std::unique_ptr<ResolutionBase> fResPhiLightJet_eta2;
  std::unique_ptr<ResolutionBase> fResPhiLightJet_eta3;
  std::unique_ptr<ResolutionBase> fResPhiLightJet_eta4;

  std::unique_ptr<ResolutionBase> fResPhiBJet_eta1;
  std::unique_ptr<ResolutionBase> fResPhiBJet_eta2;
  std::unique_ptr<ResolutionBase> fResPhiBJet_eta3;
  std::unique_ptr<ResolutionBase> fResPhiBJet_eta4;

  /// @}
};
}  // namespace KLFitter

#endif  // KLFITTER_DETECTORATLAS_8TEV_ANGULAR_H_
