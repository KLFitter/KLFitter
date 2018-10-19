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

#ifndef KLFITTER_TEST_PARTICLES_JET_H_
#define KLFITTER_TEST_PARTICLES_JET_H_

#include "gtest/gtest.h"

#include "KLFitter/Particles/Jet.h"

TEST(TestParticleJet, GetName) {
  KLFitter::Particles::Jet e{"test_name", TLorentzVector{}};
  EXPECT_EQ("test_name", e.GetName());
}

TEST(TestParticleJet, ConstructAndGetFourVector) {
  TLorentzVector p4{15, 23.4, 27, 3};

  // Test whether four vector is correctly stored at construction.
  KLFitter::Particles::Jet e{"test_name", p4};
  EXPECT_EQ(p4, e.GetP4());
  EXPECT_FLOAT_EQ(23.4, e.GetP4().Y());

  // Now test whether Jet::SetP4() works.
  p4.SetX(17.342);
  p4.SetY(12.232);
  e.SetP4(p4);
  EXPECT_EQ(p4, e.GetP4());
  EXPECT_FLOAT_EQ(12.232, e.GetP4().Y());
}

TEST(TestParticleJet, SetAndGetIdentifier) {
  KLFitter::Particles::Jet e{"", TLorentzVector{}};
  unsigned int id = 25;
  e.SetIdentifier(id);
  EXPECT_EQ(id, e.GetIdentifier());
  id = 729;
  e.SetIdentifier(id);
  EXPECT_EQ(id, e.GetIdentifier());
}

TEST(TestParticleJet, SetAndGetDetEta) {
  KLFitter::Particles::Jet e{"", TLorentzVector{}};
  double eta = 23.523;
  e.SetDetEta(eta);
  EXPECT_EQ(eta, e.GetDetEta());
}

TEST(TestParticleJet, SetAndGetBTagEfficiency) {
  KLFitter::Particles::Jet e{"", TLorentzVector{}};
  double eff = 82.363;
  e.SetBTagEfficiency(eff);
  EXPECT_DOUBLE_EQ(eff, e.GetBTagEfficiency());
}

TEST(TestParticleJet, SetAndGetBTagRejection) {
  KLFitter::Particles::Jet e{"", TLorentzVector{}};
  double eff = 12.363;
  e.SetBTagRejection(eff);
  EXPECT_DOUBLE_EQ(eff, e.GetBTagRejection());
}

TEST(TestParticleJet, SetAndGetBTagWeight) {
  KLFitter::Particles::Jet e{"", TLorentzVector{}};
  EXPECT_EQ(false, e.GetBTagWeightIsSet());
  double weight = -.3921;
  e.SetBTagWeight(weight);
  EXPECT_EQ(true, e.GetBTagWeightIsSet());
  EXPECT_DOUBLE_EQ(weight, e.GetBTagWeight());
}

TEST(TestParticleJet, SetAndGetBTagWeightIsSet) {
  KLFitter::Particles::Jet e{"", TLorentzVector{}};
  EXPECT_EQ(false, e.GetBTagWeightIsSet());
  e.SetBTagWeightIsSet(true);
  EXPECT_EQ(true, e.GetBTagWeightIsSet());
  e.SetBTagWeightIsSet(false);
  EXPECT_EQ(false, e.GetBTagWeightIsSet());
}

TEST(TestParticleJet, SetAndGetIsBTagged) {
  KLFitter::Particles::Jet e{"", TLorentzVector{}};
  EXPECT_EQ(false, e.GetIsBTagged());
  e.SetIsBTagged(true);
  EXPECT_EQ(true, e.GetIsBTagged());
  e.SetIsBTagged(false);
  EXPECT_EQ(false, e.GetIsBTagged());
}

TEST(TestParticleJet, SetAndGetTrueFlavor) {
  KLFitter::Particles::Jet e{"", TLorentzVector{}};
  EXPECT_EQ(KLFitter::Particles::JetTrueFlavor::kNone, e.GetTrueFlavor());
  e.SetTrueFlavor(KLFitter::Particles::JetTrueFlavor::kB);
  EXPECT_EQ(KLFitter::Particles::JetTrueFlavor::kB, e.GetTrueFlavor());
  e.SetTrueFlavor(KLFitter::Particles::JetTrueFlavor::kLight);
  EXPECT_EQ(KLFitter::Particles::JetTrueFlavor::kLight, e.GetTrueFlavor());
}

#endif  // KLFITTER_TEST_PARTICLES_JET_H_
