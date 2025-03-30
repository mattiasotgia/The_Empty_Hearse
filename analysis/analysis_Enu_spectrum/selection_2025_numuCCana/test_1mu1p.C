#include "helper_selection.h"
#include "sbnana/CAFAna/Core/Tree.h"

const ana::SpillCut valid_events ([](const caf::SRSpillProxy *spill) -> bool {

    std::vector<unsigned> golden_events = {15428};

    // return std::find(good_events.begin(), good_events.end(), static_cast<int>(event(spill))) != good_events.end();
    return std::find(golden_events.begin(), golden_events.end(), static_cast<int>(event(spill))) != golden_events.end();
});

const ana::SpillVar spill_dE ([](const caff::SRSpillProxy *spill) -> double {
    int ipfp_muon = -1;
    double value = -9999;
    for (auto const& slice: spill->slc) {
        if (!automatic_selection_1muNp(spill, slice, 10, 100))
            continue;

        ipfp_muon = automatic_selection_mu_index (spill, slice, 10, 100);
        if (ipfp_muon == -1)
            continue;
        
        if (classification_type (spill, slice) != 2) 
            continue;
        
        // COMPUTE VALUE HERE
        double reco_E = Neutrino_energy_reco_Np (slice, ipfp_muon, 10);
        double true_E = slice.truth.E;
        if (std::isnan(true_E) || std::isnan(reco_E))
            continue;
        
        value = reco_E - true_E
    }
});

void test_1muNp () {
    ana::SpectrumLoader loader_nominal("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    
}