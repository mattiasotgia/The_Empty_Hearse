#pragma once

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"

#include "TProfile.h"

#include <vector>     
#include <tuple>
#include <set>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

#define ELECTRON 11
#define PLANE 2

namespace huntingMichelsVars {

  std::vector<std::pair<double, double>> getEnergyTrackscore(const caf::SRSpillProxy* sr) {

    std::vector<std::pair<double, double>> returnStuff;

    for (auto const& true_particle: sr->true_particles) {

      bool isElectron = (ELECTRON == std::abs(true_particle.pdg));

      if (not isElectron) continue;
      if (not true_particle.contained) continue;
      if (true_particle.cryostat == -1) continue;

      int trueG4ID = true_particle.G4ID;

      // Loop through reco (all) to match with best match
      for (auto const& slice: sr->slc) {
        for (auto const& pfp: slice.reco.pfp) {
          if (pfp.trk.truth.bestmatch.G4ID == trueG4ID) {
            // particle match!
            returnStuff.emplace_back(std::make_pair(true_particle.plane[true_particle.cryostat][PLANE].visE, pfp.trackScore));
          }
        }
      }
    }

    return returnStuff;
  
  }

  const ana::SpillMultiVar E ([](const caf::SRSpillProxy* sr) -> std::vector<double> {
    std::vector<double> returnE;
    std::vector<std::pair<double, double>> commonReturn = getEnergyTrackscore(sr);
    for (auto const& [E, score]: commonReturn) {
	    returnE.emplace_back(E);
    }
    return returnE;
  });

  const ana::SpillMultiVar score ([](const caf::SRSpillProxy* sr) -> std::vector<double> {
    std::vector<double> returnScore;
    std::vector<std::pair<double, double>> commonReturn = getEnergyTrackscore(sr);
    for (auto const& [E, score]: commonReturn) {
	    returnScore.emplace_back(score);
    }
    return returnScore;
  });

} // huntingMichels
