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

  struct electronVariables {
    double E, score, completeness, purity;
  };

  std::vector<electronVariables> getEnergyTrackscore(const caf::SRSpillProxy* sr) {

    std::vector<electronVariables> returnStuff;

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
            returnStuff.emplace_back(electronVariables({
              // true_particle.plane[true_particle.cryostat][PLANE].visE,
              true_particle.genE,
              pfp.trackScore, 
              pfp.trk.truth.bestmatch.hit_completeness, 
              pfp.trk.truth.bestmatch.hit_purity
            }));
          }
        }
      }
    }

    return returnStuff;
  
  }

  const ana::SpillMultiVar E ([](const caf::SRSpillProxy* sr) -> std::vector<double> {
    std::vector<double> returnE;
    std::vector<electronVariables> commonReturn = getEnergyTrackscore(sr);
    for (auto const& variable: commonReturn) {
	    returnE.emplace_back(variable.E);
    }
    return returnE;
  });

  const ana::SpillMultiVar score ([](const caf::SRSpillProxy* sr) -> std::vector<double> {
    std::vector<double> returnScore;
    std::vector<electronVariables> commonReturn = getEnergyTrackscore(sr);
    for (auto const& variable: commonReturn) {
	    returnScore.emplace_back(variable.score);
    }
    return returnScore;
  });

  const ana::SpillMultiVar completeness ([](const caf::SRSpillProxy* sr) -> std::vector<double> {
    std::vector<double> returnScore;
    std::vector<electronVariables> commonReturn = getEnergyTrackscore(sr);
    for (auto const& variable: commonReturn) {
	    returnScore.emplace_back(variable.completeness);
    }
    return returnScore;
  });

  const ana::SpillMultiVar purity ([](const caf::SRSpillProxy* sr) -> std::vector<double> {
    std::vector<double> returnScore;
    std::vector<electronVariables> commonReturn = getEnergyTrackscore(sr);
    for (auto const& variable: commonReturn) {
	    returnScore.emplace_back(variable.purity);
    }
    return returnScore;
  });

} // huntingMichels
