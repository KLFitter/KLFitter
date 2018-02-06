#include "KLFitter/LikelihoodTopLeptonJets.h"

#include <iostream>

int main() {
  const float met{26125.748};
  const float met_PHI{0.3639200};

  const std::string lep_type{"Muon"};
  const float lep_pt{30501.886};
  const float lep_eta{0.4483959};
  const float lep_phi{2.9649317};
  const float lep_E{33620.113};

  const float jet1_pt{133569.53};
  const float jet1_eta{0.2231264};
  const float jet1_phi{1.7798618};
  const float jet1_E{137562.92};
  const float jet1_btag_weight{0.6868029};
  const bool jet1_has_btag{false};

  const float jet2_pt{77834.281};
  const float jet2_eta{0.8158330};
  const float jet2_phi{-1.533635};
  const float jet2_E{105723.34};
  const float jet2_btag_weight{-0.869940};
  const bool jet2_has_btag{false};

  const float jet3_pt{49327.293};
  const float jet3_eta{1.9828589};
  const float jet3_phi{-1.878274};
  const float jet3_E{182640.06};
  const float jet3_btag_weight{0.9999086};
  const bool jet3_has_btag{true};

  const float jet4_pt{43140.816};
  const float jet4_eta{0.4029131};
  const float jet4_phi{-0.472721};
  const float jet4_E{47186.804};
  const float jet4_btag_weight{-0.223728};
  const bool jet4_has_btag{false};

  KLFitter::LikelihoodTopLeptonJets lh{};

  std::cout << "Hello World" << std::endl;
  return 0;
}
