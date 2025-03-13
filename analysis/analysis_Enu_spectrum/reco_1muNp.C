
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

#include "helper.h"
#include "selection.h"

template <class T> 
struct plot1D {
    std::string name;
    const ana::Binning bin;
    T var;
};

struct spectra {
    ana::Spectrum *spectra;
    std::string name;
};


void reco_1muNp() {
    
    // REMARK: not all events are in both processing, some are only in the 
    // not cheated (17700 evt), nd not in the cheated sample (17500)...
    // These are accounted for in the hard_code namespace (hard_code::is_bad_event...)
    ana::SpectrumLoader loader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    ana::SpectrumLoader loader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_stage1tocaf_flatcafs");



};

#endif