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

#include "gtest/gtest.h"

#include "KLFitter/Particles/Parton.h"
#include "TLorentzVector.h"

namespace KLFitter {
TEST(TestParticleParton, GetName) {
  KLFitter::Particles::Parton e{"test_name", TLorentzVector{}};
  EXPECT_EQ("test_name", e.GetName());
}

TEST(TestParticleParton, ConstructAndGetFourVector) {
  TLorentzVector p4{15, 23.4, 27, 3};

  // Test whether four vector is correctly stored at construction.
  KLFitter::Particles::Parton e{"test_name", p4};
  EXPECT_EQ(p4, e.GetP4());
  EXPECT_FLOAT_EQ(23.4, e.GetP4().Y());

  // Now test whether Parton::SetP4() works.
  p4.SetX(17.342);
  p4.SetY(12.232);
  e.SetP4(p4);
  EXPECT_EQ(p4, e.GetP4());
  EXPECT_FLOAT_EQ(12.232, e.GetP4().Y());
}

TEST(TestParticleParton, SetAndGetIdentifier) {
  KLFitter::Particles::Parton e{"", TLorentzVector{}};
  int id = 25;
  e.SetIdentifier(id);
  EXPECT_EQ(id, e.GetIdentifier());
  id = 729;
  e.SetIdentifier(id);
  EXPECT_EQ(id, e.GetIdentifier());
}

TEST(TestParticleParton, SetAndGetDetEta) {
  KLFitter::Particles::Parton e{"", TLorentzVector{}};
  double eta = 23.523;
  e.SetDetEta(eta);
  EXPECT_EQ(eta, e.GetDetEta());
}

TEST(TestParticleParton, SetAndGetBTagEfficiency) {
  KLFitter::Particles::Parton e{"", TLorentzVector{}};
  double eff = 82.363;
  e.SetBTagEfficiency(eff);
  EXPECT_DOUBLE_EQ(eff, e.GetBTagEfficiency());
}

TEST(TestParticleParton, SetAndGetBTagRejection) {
  KLFitter::Particles::Parton e{"", TLorentzVector{}};
  double eff = 12.363;
  e.SetBTagRejection(eff);
  EXPECT_DOUBLE_EQ(eff, e.GetBTagRejection());
}

TEST(TestParticleParton, SetAndGetBTagWeight) {
  KLFitter::Particles::Parton e{"", TLorentzVector{}};
  EXPECT_FALSE(e.GetBTagWeightIsSet());
  double weight = -.3921;
  e.SetBTagWeight(weight);
  EXPECT_TRUE(e.GetBTagWeightIsSet());
  EXPECT_DOUBLE_EQ(weight, e.GetBTagWeight());
}

TEST(TestParticleParton, SetAndGetBTagWeightIsSet) {
  KLFitter::Particles::Parton e{"", TLorentzVector{}};
  EXPECT_FALSE( e.GetBTagWeightIsSet());
  e.SetBTagWeightIsSet(true);
  EXPECT_TRUE(e.GetBTagWeightIsSet());
  e.SetBTagWeightIsSet(false);
  EXPECT_FALSE(e.GetBTagWeightIsSet());
}

TEST(TestParticleParton, SetAndGetIsBTagged) {
  KLFitter::Particles::Parton e{"", TLorentzVector{}};
  EXPECT_FALSE(e.GetIsBTagged());
  e.SetIsBTagged(true);
  EXPECT_TRUE(e.GetIsBTagged());
  e.SetIsBTagged(false);
  EXPECT_FALSE(e.GetIsBTagged());
}

TEST(TestParticleParton, SetAndGetTrueFlavor) {
  KLFitter::Particles::Parton e{"", TLorentzVector{}};
  EXPECT_EQ(KLFitter::Particles::PartonTrueFlavor::kNone, e.GetTrueFlavor());
  e.SetTrueFlavor(KLFitter::Particles::PartonTrueFlavor::kB);
  EXPECT_EQ(KLFitter::Particles::PartonTrueFlavor::kB, e.GetTrueFlavor());
  e.SetTrueFlavor(KLFitter::Particles::PartonTrueFlavor::kLight);
  EXPECT_EQ(KLFitter::Particles::PartonTrueFlavor::kLight, e.GetTrueFlavor());
}
}  // namespace KLFitter
