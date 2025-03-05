#ifndef HELPER_H
#define HELPER_H

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Core/Var.h"

#include <vector>     
#include <tuple>
#include <set>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

const ana::SpillVar run    = ana::SIMPLESPILLVAR(hdr.run);
const ana::SpillVar subrun = ana::SIMPLESPILLVAR(hdr.subrun);
const ana::SpillVar event  = ana::SIMPLESPILLVAR(hdr.evt)

namespace log {

    enum level {
        debug, info, warning, error
    };

    std::ostream log(level L) {
        return std::cout << "[" << L << "] ";
    }
}

namespace bins {
    ana::Binning simple = Binning::Simple(30, 0, 1);
} // namespace bins

namespace hard_code {

    std::vector<unsigned> bad_events {
        10201, 10202, 10203, 10204, 10205, 10206, 10207, 10208, 10209, 10210, 10211, 10212, 10213, 10214, 10215, 10216, 10217, 10218, 10219, 10220, 10221, 10222, 10223, 10224, 10225, 10226, 10227, 10228, 10229, 10230, 10231, 10232, 10233, 10234, 10235, 10236, 10237, 10238, 10239, 10240, 10241, 10242, 10243, 10244, 10245, 10246, 10247, 10248, 10249, 10250, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 12551, 12552, 12553, 12554, 12555, 12556, 12557, 12558, 12559, 12560, 12561, 12562, 12563, 12564, 12565, 12566, 12567, 12568, 12569, 12570, 12571, 12572, 12573, 12574, 12575, 12576, 12577, 12578, 12579, 12580, 12581, 12582, 12583, 12584, 12585, 12586, 12587, 12588, 12589, 12590, 12591, 12592, 12593, 12594, 12595, 12596, 12597, 12598, 12599, 12600, 14951, 14952, 14953, 14954, 14955, 14956, 14957, 14958, 14959, 14960, 14961, 14962, 14963, 14964, 14965, 14966, 14967, 14968, 14969, 14970, 14971, 14972, 14973, 14974, 14975, 14976, 14977, 14978, 14979, 14980, 14981, 14982, 14983, 14984, 14985, 14986, 14987, 14988, 14989, 14990, 14991, 14992, 14993, 14994, 14995, 14996, 14997, 14998, 14999, 15000
    };

    const bool is_bad_event(const caf::SRSpillProxy *sr) {
        return std::find(bad_events.begins(), bad_events.end(), std::static_cast<int>(event(sr))) != bad_events.end();
    }
} // namespace hard_code

namespace cheating {
 
    // 
    // Definition of a chi2 statistics to test the cheat level (of the reconstruction goodness)
    template<class T>
    const T chi2(T __a, T __b, std::string label = "undef") {
        double chi2=0; // default "as if truthmatched";
        if (__b == 0 ) {
            if (__a == __b) {
                log(warning) << "Found truth value zero (checking label " 
                             << label << "), but reco == truth, so this will be ignored";
                return chi2;
            } else {
                log(warning) << "Found truth value zero (checking label " 
                             << label << "), but reco != truth, so returning their difference squared...";
                return std::pow(std::abs(__a - __b), 2);
            }
        }
        chi2 = std::pow(std::abs(__a - __b)/(double)__b, 2);
        return chi2;
    }


    // The idea is to build a strong estimate (chi2-like) to test the goodness of 
    // this cheating... This is integrated over all the variables
    // At the moment the only possibility is ( V_reco - V_truth ) / V_truth
    // Minimal selection applied: iscc + pdg == 14 + event is not in the bad_list...
    const ana::SpillMultiVar score ([](const caf::SRSpillProxy *spill) -> std::vector<double> {

        std::vector<double> score;

        // Exclude bad_events numbers...
        if (hard_code::is_bad_event(sr)) return score;

        for (auto const& slice: spill->slc) {

            // ...and slices that are not true numuCC {?}
            if (!slice.truth.iscc || slice.truth.pdg != 14) continue;
        
            // Checking neutrino score
            // Check vertex...
            // score += chi2<double>(slice->vertex.x, slice->truth.vtx.x);
            // score += chi2<double>(slice->vertex.y, slice->truth.vtx.y);
            // score += chi2<double>(slice->vertex.z, slice->truth.vtx.z);

            // score.push_back(chi2<double>(slice.vertex.x, slice.truth.vtx.x, "vertex"));
            // score.push_back(chi2<double>(slice.vertex.y, slice.truth.vtx.y, "vertex"));
            // score.push_back(chi2<double>(slice.vertex.z, slice.truth.vtx.z, "vertex"));

            score.push_back(
                chi2<double>(
                    sqrt(pow(slice.vertex.x, 2) + pow(slice.vertex.y, 2) + pow(slice.vertex.z, 2)), 
                    sqrt(pow(slice.truth.vtx.x, 2) + pow(slice.truth.vtx.y, 2) + pow(slice.truth.vtx.z, 2)),
                    "vertex"
                )
            );


            // MAybe has some sense?
            // score.push_back(chi2<double>(slice->tmatch.eff, 1, "slicing eff..."));

            // recob::PFParticle and other stuff
            for (auto const &pfp: slice.pfp) {
                
            } // loop pfp
        } // loop slics

    }); // score ana::SpillMultiVar

    // const ana::Cut select_nuCC([](const caf::Proxy<caf::SRSlice> *slice) -> double {

    // });
} // namespace cheating

#endif // HELPER_H