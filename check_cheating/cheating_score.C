
#ifndef CHEATING_SCORE_C
#define CHEATING_SCORE_C

#include "helper.h"

template <class T> struct plot1D {
    std::string name;
    const ana::Binning bin;
    ana::SpectrumLoader loader;
    const T var;
};

struct spectra {
    Spectrum *spectra;
    std::string name;
};

void cheating_score() {

    // REMARK: not all events are in both processing, some are only in the 
    // not cheated (17700 evt), nd not in the cheated sample (17500)...
    // These are accounted for in the hard_code namespace (hard_code::is_bad_event...)
    ana::SpectrumLoader dataloader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_stage1tocaf_flatcafs");
    ana::SpectrumLoader dataloader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_stage1tocaf_flatcafs");
    
    std::vector<plot1D<ana::SpillMultiVar>> plots {
        {"cheating_checks", bins::simple, dataloader_cheated, cheating::score}, 
        {"non_cheated_checks", bins::simple, dataloader_non_cheated, cheating::score}
    };

    std::vector<spectra> cheating_spectra;

    for (auto const& plot: plots) {
        cheating_spectra.push_back({
            new ana::Spectrum(plot.name, plot.bin, plot.loader, plot.var, ana::kNoSpillCut), plot.name
        });
    }

    return;
}

#endif // CHEATING_SCORE_C
    