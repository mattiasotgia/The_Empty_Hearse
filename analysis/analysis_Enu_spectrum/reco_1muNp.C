
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

#include "helper.h"
#include "selection.h"

#include "sbnana/CAFAna/Core/Tree.h"

const ana::Cut def_cut = (
    // ana::kNoCut
    cuts::reco::slice_1muNp             &&
    cuts::reco::slice_at_least_mu       &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_mu_50_length      &&
    cuts::reco::slice_mu_in_contained   &&
    cuts::reco::slice_mu_not_crossing   &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

using level_t = logger::level;

const ana::SpillVar spill_reco_E                = var_utils::make_spill_from_slice (vars::reco::slice_neutrino_energy_1muNp,    def_cut);
const ana::SpillVar spill_reco_pT               = var_utils::make_spill_from_slice (vars::reco::slice_neutrino_pT_1muNp,        def_cut);
const ana::SpillVar spill_vertex_difference     = var_utils::make_spill_from_slice (vars::reco::slice_vertex_difference,        def_cut);
const ana::SpillVar spill_muon_purity           = var_utils::make_spill_from_slice (vars::reco::slice_muon_hit_purity,          def_cut);
const ana::SpillVar spill_muon_completeness     = var_utils::make_spill_from_slice (vars::reco::slice_muon_hit_completeness,    def_cut);

void reco_1muNp() {
    
    // REMARK: not all events are in both processing, some are only in the 
    // not cheated (17700 evt), nd not in the cheated sample (17500)...
    // These are accounted for in the hard_code namespace (hard_code::is_bad_event...,
    // accessed through cheating::cut_bad_events)

    ana::SpectrumLoader loader_cheated("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-slice0-problem/caf/stage1_full_cheat.flat.caf.root");
    // ana::SpectrumLoader loader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_stage1tocaf_flatcafs");
    std::unique_ptr<ana::Tree> cheated(new ana::Tree(
        "cheated", 
        {"reco_E", "muon_hit_completeness", "muon_hit_purity", "vertex_difference"}, 
        loader_cheated, 
        {spill_reco_E, spill_muon_completeness, spill_muon_purity, spill_vertex_difference}, 
        // {vars::reco::slice_neutrino_energy_1muNp},
        cheating::cut_bad_events
        // ana::kNoCut // def_cut
        // def_cut
    )); 
    loader_cheated.Go();
    
    ana::SpectrumLoader loader_non_cheated("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-slice0-problem/caf/stage1_no_cheat.flat.caf.root");
    // ana::SpectrumLoader loader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    std::unique_ptr<ana::Tree> non_cheated(new ana::Tree(
        "non_cheated", 
        {"reco_E", "muon_hit_completeness", "muon_hit_purity", "vertex_difference"}, 
        loader_non_cheated, 
        {spill_reco_E, spill_muon_completeness, spill_muon_purity, spill_vertex_difference}, 
        // {vars::reco::slice_neutrino_energy_1muNp},
	    cheating::cut_bad_events
        // ana::kNoCut // def_cut
        // def_cut
    )); 
    loader_non_cheated.Go();

    std::unique_ptr<TFile> file(new TFile("reco1muNp_alldata.root", "RECREATE"));
    std::cout << "__ WRITING :     cheated __ " << std::endl;
    cheated->SaveTo(file->mkdir("cheated"));
    std::cout << "__ WRITING : non cheated __ " << std::endl;
    non_cheated->SaveTo(file->mkdir("non_cheated"));
};

#endif