
#ifndef SLICE_HELPER_H
#define SLICE_HELPER_H

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TMath.h"

#include "helper.h"
#include "selection.h"
#include "pdg_helper.h"

#include <iostream>
#include <string>
#include <vector>

namespace vars {
    namespace slice {
        const ana::Var slice_efficiency ([](const caf::SRSliceProxy* slice) -> double {
            return slice->tmatch.eff;
        });

        const ana::Var slice_purity ([](const caf::SRSliceProxy* slice) -> double {
            return slice->tmatch.pur;
        });
    }

    namespace slice_proton {
        const ana::MultiVar slice_efficiency ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> slcEff;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && vars::pdg::zero_vtXDist(slice, iPfp)) {
                    slcEff.push_back(slice->tmatch.eff);
                }
                iPfp++;
            }
            return slcEff;
        });

        const ana::MultiVar slice_purity ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> slcPur;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && vars::pdg::zero_vtXDist(slice, iPfp)) {
                    slcPur.push_back(slice->tmatch.pur);
                }
                iPfp++;
            }
            return slcPur;
        });
    }
}

#endif
