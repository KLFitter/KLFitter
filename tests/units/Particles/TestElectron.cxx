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

#include "KLFitter/Particles/Electron.h"
#include "TLorentzVector.h"

namespace KLFitter {
TEST(TestParticleElectron, GetName) {
  KLFitter::Particles::Electron e{"test_name", TLorentzVector{}};
  EXPECT_EQ("test_name", e.GetName());
}

TEST(TestParticleElectron, ConstructAndGetFourVector) {
  TLorentzVector p4{15, 23.4, 27, 3};

  // Test whether four vector is correctly stored at construction.
  KLFitter::Particles::Electron e{"test_name", p4};
  EXPECT_EQ(p4, e.GetP4());
  EXPECT_FLOAT_EQ(23.4, e.GetP4().Y());

  // Now test whether Electron::SetP4() works.
  p4.SetX(17.342);
  p4.SetY(12.232);
  e.SetP4(p4);
  EXPECT_EQ(p4, e.GetP4());
  EXPECT_FLOAT_EQ(12.232, e.GetP4().Y());
}

TEST(TestParticleElectron, SetAndGetIdentifier) {
  KLFitter::Particles::Electron e{"", TLorentzVector{}};
  int id = 25;
  e.SetIdentifier(id);
  EXPECT_EQ(id, e.GetIdentifier());
  id = 729;
  e.SetIdentifier(id);
  EXPECT_EQ(id, e.GetIdentifier());
}

TEST(TestParticleElectron, SetAndGetDetEta) {
  KLFitter::Particles::Electron e{"", TLorentzVector{}};
  double eta = 23.523;
  e.SetDetEta(eta);
  EXPECT_EQ(eta, e.GetDetEta());
}

TEST(TestParticleElectron, SetAndGetCharge) {
  KLFitter::Particles::Electron e{"", TLorentzVector{}};
  float charge{-12.3};
  e.SetCharge(charge);
  EXPECT_EQ(charge, e.GetCharge());
}
}  // namespace KLFitter 
