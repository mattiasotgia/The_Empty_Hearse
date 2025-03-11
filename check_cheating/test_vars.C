#include "helper.h"

#include "TCanvas.h"
#include "TH1D.h"

#include "TStyle.h"
#include "TLegend.h"

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

const ana::SpillCut valid_events ([](const caf::SRSpillProxy *spill) -> bool {

    // std::vector<unsigned> good_events = {
    //     5, 6, 7, 11, 13, 17, 19, 21, 22, 23, 24, 25, 29, 30, 31, 35, 36, 39, 41
    // };
    // 41: vertex is out of the TPC

    std::vector<unsigned> good_events = {
        5, 6, 7, 11, 13, 17, 19, 21, 22, 23, 24, 25, 29, 30, 31, 35, 36, 39
    };

    // Those are two events 21: 1µ1p0π; 35: 1µNpMπ (N=2, M=3) not interesting for analysis
    std::vector<unsigned> golden_events = {
        19, 
        21, 
        35
    };

    // return std::find(good_events.begin(), good_events.end(), static_cast<int>(event(spill))) != good_events.end();
    return std::find(golden_events.begin(), golden_events.end(), static_cast<int>(event(spill))) != golden_events.end();
});

void test_vars() {

    // ana::SpectrumLoader loader("/exp/icarus/data/users/msotgia/cheating-tests/test-single-event/caf/stage1_allpass.flat.caf.root");
    // ana::SpectrumLoader loader("/exp/icarus/data/users/msotgia/cheating-tests/test-single-event/caf/stage1_cheat_upto_2d.flat.caf.root");
    // ana::SpectrumLoader loader("/exp/icarus/data/users/msotgia/cheating-tests/test-single-event/caf/stage1_cheat_upto_vertex.flat.caf.root");
    // ana::SpectrumLoader loader("/exp/icarus/data/users/msotgia/cheating-tests/test-single-event/caf/stage1_cheat_upto_3d.flat.caf.root");
    ana::SpectrumLoader loader("/exp/icarus/data/users/msotgia/cheating-tests/test-single-event/caf/stage1_cheat_upto_neutrino.flat.caf.root");
    // ana::SpectrumLoader loader("/exp/icarus/data/users/msotgia/cheating-tests/test-single-event/caf/stage1_cheat_upto_bdt.flat.caf.root");

    
    ana::Tree("test", {"test_var"}, loader, {cheating::test_variables}, cheating::cut_bad_events && valid_events); 
    loader.Go();

    return;
}   
