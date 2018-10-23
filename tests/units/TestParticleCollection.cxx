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

#include "KLFitter/ParticleCollection.h"
#include "KLFitter/Particles/Boson.h"
#include "KLFitter/Particles/Electron.h"
#include "KLFitter/Particles/Parton.h"
#include "KLFitter/Particles/Muon.h"
#include "KLFitter/Particles/Neutrino.h"
#include "KLFitter/Particles/Photon.h"
#include "KLFitter/Particles/Tau.h"
#include "KLFitter/Particles/Track.h"
#include "TLorentzVector.h"

#include <vector>

namespace KLFitter {
TEST(TestParticleCollection, CopyAssign) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Electron p{"test_p", vec};
  coll.electrons.emplace_back(p);

  // Now make a copy of the above object.
  ParticleCollection new_coll{coll};
  EXPECT_EQ(new_coll.electrons.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(new_coll.electrons.at(0).GetP4().X(), p.GetP4().X());
  EXPECT_EQ(new_coll.electrons.size(), (size_t)1);
  EXPECT_EQ(new_coll.muons.size(), (size_t)0);
}

TEST(TestParticleCollection, AddParton) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Parton p{"test_p", vec};
  coll.AddParticle(p);

  // Test whether adding the particle worked.
  EXPECT_EQ(coll.partons.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(coll.partons.at(0).GetP4().X(), p.GetP4().X());

  // This should throw an exception, because you can't add
  // particles with identical names.
  EXPECT_THROW(coll.AddParticle(p), std::invalid_argument);
}

TEST(TestParticleCollection, AddElectron) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Electron p{"test_p", vec};
  coll.AddParticle(p);

  // Test whether adding the particle worked.
  EXPECT_EQ(coll.electrons.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(coll.electrons.at(0).GetP4().X(), p.GetP4().X());

  // This should throw an exception, because you can't add
  // particles with identical names.
  EXPECT_THROW(coll.AddParticle(p), std::invalid_argument);
}

TEST(TestParticleCollection, AddMuon) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Muon p{"test_p", vec};
  coll.AddParticle(p);

  // Test whether adding the particle worked.
  EXPECT_EQ(coll.muons.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(coll.muons.at(0).GetP4().X(), p.GetP4().X());

  // This should throw an exception, because you can't add
  // particles with identical names.
  EXPECT_THROW(coll.AddParticle(p), std::invalid_argument);
}

TEST(TestParticleCollection, AddPhoton) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Photon p{"test_p", vec};
  coll.AddParticle(p);

  // Test whether adding the particle worked.
  EXPECT_EQ(coll.photons.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(coll.photons.at(0).GetP4().X(), p.GetP4().X());

  // This should throw an exception, because you can't add
  // particles with identical names.
  EXPECT_THROW(coll.AddParticle(p), std::invalid_argument);
}

TEST(TestParticleCollection, AddTau) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Tau p{"test_p", vec};
  coll.AddParticle(p);

  // Test whether adding the particle worked.
  EXPECT_EQ(coll.taus.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(coll.taus.at(0).GetP4().X(), p.GetP4().X());

  // This should throw an exception, because you can't add
  // particles with identical names.
  EXPECT_THROW(coll.AddParticle(p), std::invalid_argument);
}

TEST(TestParticleCollection, AddNeutrino) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Neutrino p{"test_p", vec};
  coll.AddParticle(p);

  // Test whether adding the particle worked.
  EXPECT_EQ(coll.neutrinos.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(coll.neutrinos.at(0).GetP4().X(), p.GetP4().X());

  // This should throw an exception, because you can't add
  // particles with identical names.
  EXPECT_THROW(coll.AddParticle(p), std::invalid_argument);
}

TEST(TestParticleCollection, AddBoson) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Boson p{"test_p", vec};
  coll.AddParticle(p);

  // Test whether adding the particle worked.
  EXPECT_EQ(coll.bosons.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(coll.bosons.at(0).GetP4().X(), p.GetP4().X());

  // This should throw an exception, because you can't add
  // particles with identical names.
  EXPECT_THROW(coll.AddParticle(p), std::invalid_argument);
}

TEST(TestParticleCollection, AddTrack) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Track p{"test_p", vec};
  coll.AddParticle(p);

  // Test whether adding the particle worked.
  EXPECT_EQ(coll.tracks.at(0).GetName(), p.GetName());
  EXPECT_FLOAT_EQ(coll.tracks.at(0).GetP4().X(), p.GetP4().X());

  // This should throw an exception, because you can't add
  // particles with identical names.
  EXPECT_THROW(coll.AddParticle(p), std::invalid_argument);
}

TEST(TestParticleCollection, RemoveParticle) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Parton p{"test_p", vec};
  coll.AddParticle(p);

  // First make sure it's actually there.
  EXPECT_EQ(coll.partons.at(0).GetP4(), p.GetP4());

  // Now remove it and try to access it again.
  coll.RemoveParticle(Particles::Type::kParton, 0);
  EXPECT_THROW(coll.partons.at(0).GetP4(), std::out_of_range);
}

TEST(TestParticleCollection, FindArbitraryParticle) {
  ParticleCollection coll{};
  Particles::Electron electron{"test_electron", TLorentzVector{}};
  electron.SetCharge(-1);
  coll.AddParticle(electron);
  Particles::Parton parton{"test_parton", TLorentzVector{}};
  parton.SetBTagWeight(0.375);
  coll.AddParticle(parton);

  // Now try to find these two particles.
  auto result = dynamic_cast<const Particles::Electron*>(coll.FindParticle("test_electron"));
  EXPECT_EQ(electron.GetCharge(), result->GetCharge());
  auto result2 = dynamic_cast<const Particles::Parton*>(coll.FindParticle("test_parton"));
  EXPECT_EQ(parton.GetBTagWeight(), result2->GetBTagWeight());

  // Now try to find one that doesn't exist.
  EXPECT_EQ(coll.FindParticle("test"), nullptr);
}

TEST(TestParticleCollection, FindParticleOfType) {
  ParticleCollection coll{};
  Particles::Electron electron{"test_electron", TLorentzVector{}};
  electron.SetCharge(-1);
  coll.AddParticle(electron);

  // Now try to find this particle.
  auto uncasted = coll.FindParticle(Particles::Type::kElectron, "test_electron");
  auto result = dynamic_cast<const Particles::Electron*>(uncasted);
  EXPECT_EQ(electron.GetCharge(), result->GetCharge());

  // Now try to find one that doesn't exist.
  EXPECT_EQ(coll.FindParticle(Particles::Type::kBoson, "boson"), nullptr);
}

TEST(TestParticleCollection, GetFourVectors) {
  ParticleCollection coll{};
  TLorentzVector vec{15.3, 283.123, -3.123, 2391.2};
  Particles::Electron p{"test_p", vec};
  coll.electrons.emplace_back(p);

  // Now retrieve the four-vector, make sure it is the same.
  const auto&& p4 = coll.GetP4(Particles::Type::kElectron, 0);
  EXPECT_EQ(*p4, vec);

  // Now do the same with the constant four-momentum.
  const auto&& const_coll = static_cast<const ParticleCollection>(coll);
  const auto&& const_p4 = const_coll.GetP4(Particles::Type::kElectron, 0);
  EXPECT_EQ(*const_p4, vec);

  // Verify that the first p4 is not constant.
  p4->SetX(13.8);
}

TEST(TestParticleCollection, NumberOfParticles) {
  ParticleCollection coll{};
  Particles::Electron el1{"el1", TLorentzVector{}};
  coll.electrons.emplace_back(el1);
  Particles::Electron el2{"el2", TLorentzVector{}};
  coll.electrons.emplace_back(el2);
  Particles::Parton p1{"p1", TLorentzVector{}};
  coll.partons.emplace_back(p1);

  // Verify that all particle numbers are correct.
  EXPECT_EQ(coll.NParticles(), (size_t)3);
  EXPECT_EQ(coll.NParticles(Particles::Type::kElectron), (size_t)2);
  EXPECT_EQ(coll.NParticles(Particles::Type::kParton), (size_t)1);
  EXPECT_EQ(coll.NParticles(Particles::Type::kBoson), (size_t)0);

  // Remove one to see if the results are adjusted.
  coll.RemoveParticle(Particles::Type::kParton, 0);
  EXPECT_EQ(coll.NParticles(), (size_t)2);
  EXPECT_EQ(coll.NParticles(Particles::Type::kParton), (size_t)0);
}

TEST(TestParticleCollection, NumberOfBTags) {
  ParticleCollection coll{};
  Particles::Parton p1{"p1", TLorentzVector{}};
  p1.SetIsBTagged(true);
  coll.partons.emplace_back(p1);
  Particles::Parton p2{"p2", TLorentzVector{}};
  p2.SetIsBTagged(true);
  coll.partons.emplace_back(p2);

  // Also add one without b-tag to test the default value.
  Particles::Parton p3{"p3", TLorentzVector{}};
  p3.SetIsBTagged(false);
  coll.partons.emplace_back(p3);

  // Verify the correct number of b-tags is set.
  EXPECT_EQ(coll.NBTags(), (size_t)2);

  // Remove one b-tag and test if the number is adjusted.
  coll.partons.at(1).SetIsBTagged(false);
  EXPECT_EQ(coll.NBTags(), (size_t)1);
}
}  // namespace KLFitter
