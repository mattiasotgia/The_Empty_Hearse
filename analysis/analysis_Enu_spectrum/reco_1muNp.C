
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

#include "helper.h"
#include "selection.h"

#include "sbnana/CAFAna/Core/Tree.h"

const ana::Cut def_cut = (
    // ana::kNoCut
    // cuts::reco::slice_1muNp             &&
    // cuts::reco::slice_at_least_mu       &&
    // cuts::reco::slice_vtx_in_FV         &&
    // cuts::reco::slice_mu_50_length      &&
    // cuts::reco::slice_mu_in_contained   &&
    // cuts::reco::slice_mu_not_crossing   &&
    // cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_numuCC           &&
    cuts::truth::slice_vtx_in_FV
);

const var_utils::cut_type local_cut_type = var_utils::cut_type::BOTH_1muNp;

using level_t = logger::level;

const ana::SpillVar spill_dE                    = var_utils::make_spill_from_slice (vars::truth::slice_neutrino_dE,             def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_reco_E                = var_utils::make_spill_from_slice (vars::reco::slice_neutrino_energy_1muNp,    def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_true_E                = var_utils::make_spill_from_slice (vars::truth::slice_neutrino_energy,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_reco_pT               = var_utils::make_spill_from_slice (vars::reco::slice_neutrino_pT_1muNp,        def_cut, local_cut_type, def_cut_truth);

// Checks on cheating
const ana::SpillVar spill_vertex_difference     = var_utils::make_spill_from_slice (vars::reco::slice_vertex_difference,        def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_vertex_difference_x   = var_utils::make_spill_from_slice (vars::reco::slice_vertex_difference_x,      def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_vertex_difference_y   = var_utils::make_spill_from_slice (vars::reco::slice_vertex_difference_y,      def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_vertex_difference_z   = var_utils::make_spill_from_slice (vars::reco::slice_vertex_difference_z,      def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_muon_purity           = var_utils::make_spill_from_slice (vars::reco::slice_muon_hit_purity,          def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_muon_completeness     = var_utils::make_spill_from_slice (vars::reco::slice_muon_hit_completeness,    def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_proton_purity         = var_utils::make_spill_from_slice (vars::reco::slice_proton_hit_purity,        def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_proton_completeness   = var_utils::make_spill_from_slice (vars::reco::slice_proton_hit_completeness,  def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_muon_true_L       = var_utils::make_spill_from_slice (vars::reco::slice_pid_muon_true_length,     def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_muon_reco_L       = var_utils::make_spill_from_slice (vars::reco::slice_pid_muon_reco_length,     def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_true_L     = var_utils::make_spill_from_slice (vars::reco::slice_pid_proton_true_length,   def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_reco_L     = var_utils::make_spill_from_slice (vars::reco::slice_pid_proton_reco_length,   def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_muon_momentum_rangeP  = var_utils::make_spill_from_slice (vars::reco::slice_muon_momentum_rangeP,     def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_proton_momentum_rangeP= var_utils::make_spill_from_slice (vars::reco::slice_proton_momentum_rangeP,   def_cut, local_cut_type, def_cut_truth);


void reco_1muNp() {
    
    // REMARK: not all events are in both processing, some are only in the 
    // not cheated (17700 evt), nd not in the cheated sample (17500)...
    // These are accounted for in the hard_code namespace (hard_code::is_bad_event...,
    // accessed through cheating::cut_bad_events)

    // ana::SpectrumLoader loader_cheated("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-slice0-problem/caf/stage1_full_cheat.flat.caf.root");
    ana::SpectrumLoader loader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_stage1tocaf_flatcafs");
    std::unique_ptr<ana::Tree> cheated(new ana::Tree(
        "cheated", 
        {"reco_E", "true_E", "delta_E", "reco_pT"}, 
        loader_cheated, 
        {spill_reco_E, spill_true_E, spill_dE, spill_reco_pT}, 
        cheating::cut_bad_events
    )); 

    std::unique_ptr<ana::Tree> cheated_checks(new ana::Tree(
        "cheated_checks", 
        {
            "vertex_difference", 
            "vertex_difference_x", 
            "vertex_difference_y", 
            "vertex_difference_z", 
            "muon_purity", 
            "muon_completeness", 
            "proton_purity", 
            "proton_completeness", 
            "pid_muon_true_L", 
            "pid_muon_reco_L", 
            "pid_proton_true_L", 
            "pid_proton_reco_L", 
            "muon_momentum_rangeP", 
            "proton_momentum_rangeP"
        }, 
        loader_cheated, 
        {
            spill_vertex_difference, 
            spill_vertex_difference_x, 
            spill_vertex_difference_y, 
            spill_vertex_difference_z, 
            spill_muon_purity, 
            spill_muon_completeness, 
            spill_proton_purity, 
            spill_proton_completeness, 
            spill_pid_muon_true_L, 
            spill_pid_muon_reco_L, 
            spill_pid_proton_true_L, 
            spill_pid_proton_reco_L, 
            spill_muon_momentum_rangeP, 
            spill_proton_momentum_rangeP
        }, 
        cheating::cut_bad_events
    )); 
    loader_cheated.Go();
    
    // ana::SpectrumLoader loader_non_cheated("/exp/icarus/data/users/msotgia/thesis/cheating-tests/test-slice0-problem/caf/stage1_no_cheat.flat.caf.root");
    ana::SpectrumLoader loader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    std::unique_ptr<ana::Tree> non_cheated(new ana::Tree(
        "non_cheated", 
        {"reco_E", "true_E", "delta_E", "reco_pT"}, 
        loader_non_cheated, 
        {spill_reco_E, spill_true_E, spill_dE, spill_reco_pT},
	    cheating::cut_bad_events
    )); 

    std::unique_ptr<ana::Tree> non_cheated_checks(new ana::Tree(
        "non_cheated_checks", 
        {
            "vertex_difference", 
            "vertex_difference_x", 
            "vertex_difference_y", 
            "vertex_difference_z", 
            "muon_purity", 
            "muon_completeness", 
            "proton_purity", 
            "proton_completeness", 
            "pid_muon_true_L", 
            "pid_muon_reco_L", 
            "pid_proton_true_L", 
            "pid_proton_reco_L", 
            "muon_momentum_rangeP", 
            "proton_momentum_rangeP"
        }, 
        loader_non_cheated, 
        {
            spill_vertex_difference, 
            spill_vertex_difference_x, 
            spill_vertex_difference_y, 
            spill_vertex_difference_z, 
            spill_muon_purity, 
            spill_muon_completeness, 
            spill_proton_purity, 
            spill_proton_completeness, 
            spill_pid_muon_true_L, 
            spill_pid_muon_reco_L, 
            spill_pid_proton_true_L, 
            spill_pid_proton_reco_L, 
            spill_muon_momentum_rangeP, 
            spill_proton_momentum_rangeP
        }, 
        cheating::cut_bad_events
    )); 
    loader_non_cheated.Go();

    std::unique_ptr<TFile> file(new TFile("reco1muNp_alldata.root", "RECREATE"));
    std::cout << "__ WRITING :     cheated __ " << std::endl;
    cheated->SaveTo(file->mkdir("cheated"));
    cheated_checks->SaveTo(file->mkdir("cheated_checks"));
    std::cout << "__ WRITING : non cheated __ " << std::endl;
    non_cheated_checks->SaveTo(file->mkdir("non_cheated_checks"));
};

#endif