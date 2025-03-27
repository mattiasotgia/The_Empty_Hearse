
#ifndef EVENT_BY_EVENT_C
#define EVENT_BY_EVENT_C

#include "sbnana/CAFAna/Core/Tree.h"

#include "selection.h"
#include "helper.h" 


using logger::log;
using level_t = logger::level;
using cut_type_t = var_utils::cut_type;

const ana::SpillVar check_from_slice (
    const ana::Cut &reco_cut = ana::kNoCut,  
    cut_type_t what_to_cut_on = cut_type_t::RECO, 
    const ana::Cut &truth_cut = ana::kNoCut,
    std::function<particle_data::int_type(
        const caf::SRSpillProxy*, const caf::SRSliceProxy*
    )> classification = var_utils::classification_type
) {
    return ana::SpillVar([=](const caf::SRSpillProxy *spill) -> double {
        int selected_slices = 0;
        double slice_value = -9999;

        // if (what_to_cut_on == cut_type_t::MC_1muNp || what_to_cut_on == cut_type_t::MC_1mu1p) {
        //     for (auto const& nu: spill->mc.nu) {
        //         if (what_to_cut_on == cut_type_t::MC_1muNp && 
        //             classification_type_MC(spill, nu) != particle_data::int_type::true_visible_1muNp) continue;
        //         if (what_to_cut_on == cut_type_t::MC_1mu1p && 
        //             classification_type_MC(spill, nu) != particle_data::int_type::true_visible_1mu1p) continue;
                
                
        //     } // loop spill->nu
        // } // MC classification :)
        
        std::ofstream file;
        file.open("writer_log.txt", std::ios_base::app);

        file << "run.event = " << run(spill) << "." << event(spill) << ":" << std::endl;
        std::size_t islc = -1;
        for (auto const& slice: spill->slc) {
            islc++;
            file << "--> slice #: " << islc  
                 << " (efficiency = " << slice.tmatch.eff << ", purity = " << slice.tmatch.pur << ")" 
                 << std::flush;

            if (what_to_cut_on == cut_type_t::RECO && !reco_cut(&slice)) {
                file << " cut on cut_type_t::RECO && reco_cut" << std::endl;
                continue; 
            }
            if (what_to_cut_on == cut_type_t::TRUE_1muNp && (
                classification(spill, &slice) != particle_data::int_type::true_visible_1muNp || !truth_cut(&slice)
            )) {
                file << " cut on cut_type_t::TRUE_1muNp && (particle_data::int_type::true_visible_1muNp + truth_cut)" << std::endl;
                continue;
            }
            if (what_to_cut_on == cut_type_t::TRUE_1mu1p && (
                classification(spill, &slice) != particle_data::int_type::true_visible_1mu1p || !truth_cut(&slice)
            )) {
                file << " cut on cut_type_t::TRUE_1mu1p && (particle_data::int_type::true_visible_1mu1p + truth_cut)" << std::endl;
                continue;
            }
            if (what_to_cut_on == cut_type_t::BOTH_1muNp && (
                !reco_cut(&slice) || !truth_cut(&slice) || 
                classification(spill, &slice) != particle_data::int_type::true_visible_1muNp
            )) {
                file << " cut on cut_type_t::BOTH_1muNp && (particle_data::int_type::true_visible_1muNp + truth_cut + reco_cut)" << std::endl;
                continue;
            }
            if (what_to_cut_on == cut_type_t::BOTH_1mu1p && (
                !reco_cut(&slice) || !truth_cut(&slice) || 
                classification(spill, &slice) != particle_data::int_type::true_visible_1mu1p
            )) {
                file << " cut on cut_type_t::BOTH_1mu1p && (particle_data::int_type::true_visible_1mu1p + truth_cut + reco_cut)" << std::endl;
                continue;
            }
            std::cout << "--> slice #: " << islc  
                      << " (efficiency = " << slice.tmatch.eff << ", purity = " << slice.tmatch.pur << ")" 
                      << " this slice is SELECTED as 1µNp (" 
                      << "run.event = " << run(spill) << "." << event(spill) << ")" << std::endl;

            file << " this slice is SELECTED as 1µNp" << std::endl; 
            
            selected_slices ++ ;
        } // loop spill->slc
    
        if (selected_slices > 1) 
            log(level_t::error) << "Something wrong with run:event = " 
                                << run(spill) << ":" << event(spill) 
                                << " => found " << selected_slices
                                << " slice(s) 1µNp"
                                << std::endl;
        return slice_value;
    });
}

const ana::Cut def_cut = (
    // ana::kNoCut
    // cuts::reco::slice_1muNp             &&
    // cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_at_least_mu       // &&
    // cuts::reco::slice_mu_50_length      &&
    // cuts::reco::slice_mu_in_contained   &&
    // cuts::reco::slice_mu_not_crossing   &&
    // cuts::reco::slice_barycenter        &&
    // cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_numuCC           &&
    cuts::truth::slice_vtx_in_FV
);

const cut_type_t local_cut_type = cut_type_t::BOTH_1muNp;

const ana::SpillCut valid_events ([](const caf::SRSpillProxy *spill) -> bool {

    std::vector<unsigned> golden_events = {15428};

    // return std::find(good_events.begin(), good_events.end(), static_cast<int>(event(spill))) != good_events.end();
    return std::find(golden_events.begin(), golden_events.end(), static_cast<int>(event(spill))) != golden_events.end();
});

const ana::SpillVar check_events = check_from_slice (def_cut, local_cut_type, def_cut_truth);

void event_by_event () {

    // ana::SpectrumLoader loader_cheated("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-150.?-side-by-side/*/stage1_cheated*flat.caf.root");
    // ana::SpectrumLoader loader_nominal("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-150.?-side-by-side/*/stage1_nominal*flat.caf.root");

    // ana::SpectrumLoader loader_cheated("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-150.2-side-by-side/stage1/stage1_cheated.flat.caf.root");
    // ana::SpectrumLoader loader_nominal("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-150.2-side-by-side/stage1/stage1_nominal.flat.caf.root");

    // ana::SpectrumLoader loader_cheated("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-150.3-side-by-side/caf/stage1_cheated.flat.caf.root");
    // ana::SpectrumLoader loader_nominal("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-150.3-side-by-side/caf/stage1_nominal.flat.caf.root");

    ana::SpectrumLoader loader_nominal("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");

    // std::unique_ptr<ana::Tree> events_cheated(new ana::Tree("events_cheated", {"label"}, loader_cheated, {check_events}, ana::kNoSpillCut)); 
    // std::unique_ptr<ana::Tree> events_nominal(new ana::Tree("events_nominal", {"label"}, loader_nominal, {check_events}, ana::kNoSpillCut)); 

    // std::unique_ptr<ana::Tree> events_cheated(new ana::Tree("events_cheated", {"label"}, loader_cheated, {cheating::test_variables}, ana::kNoSpillCut && valid_events)); 
    std::unique_ptr<ana::Tree> events_nominal(new ana::Tree("events_nominal", {"label"}, loader_nominal, {cheating::test_variables}, ana::kNoSpillCut && valid_events)); 

    std::ofstream file;
    file.open("writer_log.txt", std::ios_base::app);
    // file << "________________ CHEATED DATA ________________" << std::endl;
    // loader_cheated.Go();
    file << "________________ NOMINAL DATA ________________" << std::endl;
    loader_nominal.Go();

    return;
} // void event_by_event

#endif // EVENT_BY_EVENT_C
