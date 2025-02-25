
#ifndef SELECTION_DEF_H
#define SELECTION_DEF_H

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
const ana::SpillVar event  = ana::SIMPLESPILLVAR(hdr.evt);

const ana::Cut is_numu_CC([](const caf::SRSliceProxy *slc) -> bool {
    return slc->truth.iscc;
});

const ana::Cut is_clear_cosmic([](const caf::SRSliceProxy *slc) -> bool {
    return slc->is_clear_cosmic;
});

std::string which_cryo(const double &x) { 
    return (x<0 ? "cryo E (cryo 0)" : "cryo W (cryo 1)"); 
}

const ana::SpillVar test([](const caf::SRSpillProxy *sr) -> double {

    int cc_cosmic = 0;
    for (const auto &islc: sr->slc) {
        if (islc.is_clear_cosmic) 
            cc_cosmic++;
    }

    std::cout << "\nIn event " << run(sr) << ":" << subrun(sr) << ":" << event(sr) << std::endl;
    std::cout << "Found # " << sr->nslc << " slice(s), of which "
              << "there are # " << cc_cosmic << " clear cosmic(s)" << std::endl;
    std::cout << "Running throught pfps in slice" << std::endl;

    for (auto const &slice: sr->slc) {
        std::cout << "--> slice # " << slice.producer  << ", which " 
                  << (slice.truth.iscc ? "is" : "is NOT") << " true nuCC" 
                  << ", with slice.primary.size() = " << slice.primary.size() << std::endl;
        if (slice.is_clear_cosmic) std::cout << "    is_clear_cosmic = true" << std::endl;
        if (!slice.is_clear_cosmic) std::cout << "    is_clear_cosmic = false" << std::endl;
        for (auto const &pfp: slice.reco.pfp) {
            std::cout << "  --> track in  " << which_cryo(pfp.trk.start.x) << ", with lenght " << pfp.trk.len << std::endl;
            std::cout << "  --> shower in " << which_cryo(pfp.shw.start.x) << ", with lenght " << pfp.shw.len << std::endl;
        }

    }
    // std::cout << " Slice ! " << std::endl;
    // if (slc->is_clear_cosmic)
    //     std::cout << " Found clear cosmic " << std::endl;
    return 0;
});

#endif 