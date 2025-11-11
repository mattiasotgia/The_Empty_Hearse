#pragma once 

#include "sbnana/CAFAna/Core/Tree.h"
#include "helper.h"

using namespace huntingMichelsVars;

void huntingMichels () {

  std::map<std::string, std::string> availableLoaders = {
    {"nominal", "msotgia_v09_89_01_01p03_bdtCheatingMichel_respun_nominal"},
    {"cheated", "msotgia_v09_89_01_01p03_bdtCheatingMichel_respun_cheated"},
    {"onlyMva", "msotgia_v09_89_01_01p03_bdtCheatingMichel_respun_onlyMva"}
  };

  std::vector<std::string> runningLoaders = {"nominal", "onlyMva"};

  std::vector<std::unique_ptr<ana::Tree>> treesElectrons;

  for (auto const& runningLoader: runningLoaders) {

    ana::SpectrumLoader loader(availableLoaders.at(runningLoader));

    treesElectrons.emplace_back(std::make_unique<ana::Tree>( 
      runningLoader.c_str(),
      std::vector<std::string>{"E", "score", "completeness", "purity"},
      loader,
      std::vector<ana::SpillMultiVar>{E, score, completeness, purity},
      ana::kNoSpillCut,
      true
    ));

    loader.Go();
  }

  std::unique_ptr<TFile> treeWriter(new TFile("huntingMichels.root", "RECREATE"));

  treeWriter->mkdir("michels");
  for (auto const& tree: treesElectrons) {
    tree->SaveTo(treeWriter->GetDirectory("michels"));
  }

  return;
}
