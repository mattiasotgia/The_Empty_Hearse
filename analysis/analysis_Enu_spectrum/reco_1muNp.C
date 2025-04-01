
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

#include "helper.h"
#include "selection.h"
#include "selection_nopid.h"

#include "sbnana/CAFAna/Core/Tree.h"

#define SPILLVAR(_def, _var, _reco, _what, _true)                                                                               \
    var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(_def, _var, _reco, _what, _true)
#define SPILLMULTIVAR(_def, _var, _reco, _what, _true)                                                                          \
    var_utils::make_spill_from_slice<ana::SpillMultiVar, ana::MultiVar, std::vector<double>>(_def, _var, _reco, _what, _true)
#define NOPID_SPILLVAR(_var, _true, _type)                                                                                      \
    no_pid::make_spill_from_slice<const ana::SpillVar, double>(_var, _true, _type)
#define NOPID_SPILLMULTIVAR(_var, _true, _type)                                                                                 \
    no_pid::make_spill_from_slice<const ana::SpillMultiVar, std::vector<double>>(_var, _true, _type)

using level_t = logger::level;
using cut_type_t = var_utils::cut_type_t;



const ana::Cut def_cut = (
    // ana::kNoCut
    // (cuts::reco::slice_1mu1p || cuts::reco::slice_1muNp) &&
    cuts::reco::slice_1mu1p             &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_vtx_in_FV
);

const cut_type_t local_cut_type = cut_type_t::BOTH_1mu1p;
const particle_data::int_type_t type = particle_data::int_type_t::true_visible_1mu1p;

const ana::SpillVar spill_dE                            = SPILLVAR(-9999, vars::truth::slice_neutrino_dE,                   def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_reco_E                        = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,          def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_true_E                        = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,               def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_reco_pT                       = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,              def_cut, local_cut_type, def_cut_truth);

// Checks on cheating
const ana::SpillVar spill_vertex_difference             = SPILLVAR(-9999, vars::reco::slice_vertex_difference,              def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_vertex_difference_x           = SPILLVAR(-9999, vars::reco::slice_vertex_difference_x,            def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_vertex_difference_y           = SPILLVAR(-9999, vars::reco::slice_vertex_difference_y,            def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_vertex_difference_z           = SPILLVAR(-9999, vars::reco::slice_vertex_difference_z,            def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_muon_purity                   = SPILLVAR(-9999, vars::reco::slice_muon_hit_purity,                def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_muon_completeness             = SPILLVAR(-9999, vars::reco::slice_muon_hit_completeness,          def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_proton_purity                 = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity,              def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_proton_completeness           = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness,        def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_muon_L_reco_true_ratio    = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_reco_true_ratio,     def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_L_reco_true_ratio  = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio,   def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_muon_P_reco_true_ratio    = SPILLVAR(-9999, vars::reco::slice_muon_P_reco_true_ratio,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_P_reco_true_ratio  = SPILLVAR(-9999, vars::reco::slice_proton_P_reco_true_ratio,       def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_muon_true_L               = SPILLVAR(-9999, vars::reco::slice_pid_muon_true_length,           def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_muon_reco_L               = SPILLVAR(-9999, vars::reco::slice_pid_muon_reco_length,           def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_true_L             = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_reco_L             = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_muon_R                    = SPILLVAR(-9999, vars::reco::slice_muon_R,                         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_R                  = SPILLVAR(-9999, vars::reco::slice_leading_proton_R,               def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_muon_momentum_rangeP          = SPILLVAR(-9999, vars::reco::slice_muon_momentum_rangeP,           def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_proton_momentum_rangeP        = SPILLVAR(-9999, vars::reco::slice_proton_momentum_rangeP,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_CT3D_rangeP                   = SPILLVAR(-9999, vars::reco::slice_CT3D_rangeP_muon_proton,        def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_CT3D_trueP                    = SPILLVAR(-9999, vars::reco::slice_CT3D_trueP_muon_proton,         def_cut, local_cut_type, def_cut_truth);

// NO PID selection (MAIN)
const ana::SpillVar      nopid_spill_vertex_difference_3D             = NOPID_SPILLVAR(no_pid::vars::slice_vertex_delta_3D,                   def_cut_truth, type);
const ana::SpillVar      nopid_spill_vertex_difference_x              = NOPID_SPILLVAR(no_pid::vars::slice_vertex_delta_x,                    def_cut_truth, type);
const ana::SpillVar      nopid_spill_vertex_difference_y              = NOPID_SPILLVAR(no_pid::vars::slice_vertex_delta_y,                    def_cut_truth, type);
const ana::SpillVar      nopid_spill_vertex_difference_z              = NOPID_SPILLVAR(no_pid::vars::slice_vertex_delta_z,                    def_cut_truth, type);
const ana::SpillVar      nopid_spill_muon_hit_purity                  = NOPID_SPILLVAR(no_pid::vars::slice_muon_hit_purity,                   def_cut_truth, type);
const ana::SpillVar      nopid_spill_leading_proton_hit_purity        = NOPID_SPILLVAR(no_pid::vars::slice_leading_proton_hit_purity,         def_cut_truth, type);
const ana::SpillVar      nopid_spill_muon_hit_completeness            = NOPID_SPILLVAR(no_pid::vars::slice_muon_hit_completeness,             def_cut_truth, type);
const ana::SpillVar      nopid_spill_leading_proton_hit_completeness  = NOPID_SPILLVAR(no_pid::vars::slice_leading_proton_hit_completeness,   def_cut_truth, type);
const ana::SpillVar      nopid_spill_muon_L_reco_true_ratio           = NOPID_SPILLVAR(no_pid::vars::slice_muon_L_reco_true_ratio,            def_cut_truth, type);
const ana::SpillVar      nopid_spill_leading_proton_L_reco_true_ratio = NOPID_SPILLVAR(no_pid::vars::slice_leading_proton_L_reco_true_ratio,  def_cut_truth, type);
const ana::SpillVar      nopid_spill_muon_P_reco_true_ratio           = NOPID_SPILLVAR(no_pid::vars::slice_muon_P_reco_true_ratio,            def_cut_truth, type);
const ana::SpillVar      nopid_spill_leading_proton_P_reco_true_ratio = NOPID_SPILLVAR(no_pid::vars::slice_leading_proton_P_reco_true_ratio,  def_cut_truth, type);
const ana::SpillVar      nopid_spill_neutrino_reco_E                  = NOPID_SPILLVAR(no_pid::vars::slice_neutrino_reco_E,                   def_cut_truth, type);
const ana::SpillVar      nopid_spill_neutrino_true_E                  = NOPID_SPILLVAR(no_pid::vars::slice_neutrino_true_E,                   def_cut_truth, type);
const ana::SpillVar      nopid_spill_neutrino_reco_dE                 = NOPID_SPILLVAR(no_pid::vars::slice_neutrino_reco_dE,                  def_cut_truth, type);
const ana::SpillVar      nopid_spill_neutrino_reco_pT                 = NOPID_SPILLVAR(no_pid::vars::slice_neutrino_reco_pT,                  def_cut_truth, type);
const ana::SpillVar      nopid_spill_muon_R                           = NOPID_SPILLVAR(no_pid::vars::slice_muon_R,                            def_cut_truth, type);
const ana::SpillVar      nopid_spill_leading_proton_R                 = NOPID_SPILLVAR(no_pid::vars::slice_leading_proton_R,                  def_cut_truth, type);
const ana::SpillMultiVar nopid_spill_protons_hit_purity               = NOPID_SPILLMULTIVAR(no_pid::vars::slice_protons_hit_purity,           def_cut_truth, type);
const ana::SpillMultiVar nopid_spill_protons_hit_completeness         = NOPID_SPILLMULTIVAR(no_pid::vars::slice_protons_hit_completeness,     def_cut_truth, type);
const ana::SpillMultiVar nopid_spill_protons_L_reco_true_ratio        = NOPID_SPILLMULTIVAR(no_pid::vars::slice_protons_L_reco_true_ratio,    def_cut_truth, type);
const ana::SpillMultiVar nopid_spill_protons_P_reco_true_ratio        = NOPID_SPILLMULTIVAR(no_pid::vars::slice_protons_P_reco_true_ratio,    def_cut_truth, type);
const ana::SpillMultiVar nopid_spill_protons_R                        = NOPID_SPILLMULTIVAR(no_pid::vars::slice_protons_R,                    def_cut_truth, type);

void reco_1muNp() {
    
    // REMARK: not all events are in both processing, some are only in the 
    // not cheated (17700 evt), nd not in the cheated sample (17500)...
    // These are accounted for in the hard_code namespace (hard_code::is_bad_event...,
    // accessed through cheating::cut_bad_events)

    ana::SpectrumLoader loader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_stage1tocaf_flatcafs");
    std::unique_ptr<ana::Tree> cheated(new ana::Tree(
        "cheated", 
        {"reco_E", "true_E", "delta_E", "reco_pT"},                 ////< 4 labels 
        loader_cheated, 
        {spill_reco_E, spill_true_E, spill_dE, spill_reco_pT},      ////< 4 variables
        cheating::cut_bad_events && cuts::reco::spill_CRTPMTNeutrino
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
            "leading_proton_purity", 
            "leading_proton_completeness", 
            "muon_L_reco_true_ratio",
            "leading_proton_L_reco_true_ratio",
            "muon_P_reco_true_ratio",
            "leading_proton_P_reco_true_ratio",
            "muon_true_L", 
            "leading_proton_true_L", 
            "muon_R",
            "leading_proton_R",
            "muon_rangeP",
            "leading_proton_rangeP",
            "CT3D_muon_leading_proton_from_rangeP",
            "CT3D_muon_leading_proton_from_trueP"
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
            spill_pid_muon_L_reco_true_ratio,
            spill_pid_proton_L_reco_true_ratio,
            spill_pid_muon_P_reco_true_ratio,
            spill_pid_proton_P_reco_true_ratio,
            spill_pid_muon_true_L,  
            spill_pid_proton_true_L, 
            spill_pid_muon_R,
            spill_pid_proton_R,
            spill_muon_momentum_rangeP,
            spill_proton_momentum_rangeP,
            spill_CT3D_rangeP,
            spill_CT3D_trueP
        }, 
        cheating::cut_bad_events && cuts::reco::spill_CRTPMTNeutrino
    )); 
    
    std::unique_ptr<ana::Tree> cheated_nopid_checks(new ana::Tree(
        "cheated_nopid_checks",
        {
            "vertex_difference_3D",
            "vertex_difference_x",
            "vertex_difference_y",
            "vertex_difference_z",
            "muon_hit_purity",
            "leading_proton_hit_purity",
            "muon_hit_completeness",
            "leading_proton_hit_completeness",
            "muon_L_reco_true_ratio",
            "leading_proton_L_reco_true_ratio",
            "muon_P_reco_true_ratio",
            "leading_proton_P_reco_true_ratio",
            "neutrino_reco_E",
            "neutrino_true_E",
            "neutrino_reco_dE",
            "neutrino_reco_pT",
            "muon_R",
            "leading_proton_R"
        }, 
        loader_cheated,
        {
            nopid_spill_vertex_difference_3D,
            nopid_spill_vertex_difference_x,
            nopid_spill_vertex_difference_y,
            nopid_spill_vertex_difference_z,
            nopid_spill_muon_hit_purity,
            nopid_spill_leading_proton_hit_purity,
            nopid_spill_muon_hit_completeness,
            nopid_spill_leading_proton_hit_completeness,
            nopid_spill_muon_L_reco_true_ratio,
            nopid_spill_leading_proton_L_reco_true_ratio,
            nopid_spill_muon_P_reco_true_ratio,
            nopid_spill_leading_proton_P_reco_true_ratio,
            nopid_spill_neutrino_reco_E,
            nopid_spill_neutrino_true_E,
            nopid_spill_neutrino_reco_dE,
            nopid_spill_neutrino_reco_pT,
            nopid_spill_muon_R,
            nopid_spill_leading_proton_R
        },
        cheating::cut_bad_events && cuts::truth::spill_1mu1p_MC && cuts::truth::spill_sliceless
    ));

    std::unique_ptr<ana::Tree> cheated_nopid_checks_multip(new ana::Tree(
        "cheated_nopid_checks_multip",
        {
            "protons_hit_purity",
            "protons_hit_completeness",
            "protons_L_reco_true_ratio",
            "protons_P_reco_true_ratio",
            "protons_R"
        }, 
        loader_cheated,
        {
            nopid_spill_protons_hit_purity,
            nopid_spill_protons_hit_completeness,
            nopid_spill_protons_L_reco_true_ratio,
            nopid_spill_protons_P_reco_true_ratio,
            nopid_spill_protons_R
        },
        cheating::cut_bad_events && cuts::truth::spill_1mu1p_MC && cuts::truth::spill_sliceless
    ));
    loader_cheated.Go();
    
    ana::SpectrumLoader loader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    std::unique_ptr<ana::Tree> non_cheated(new ana::Tree(
        "non_cheated", 
        {"reco_E", "true_E", "delta_E", "reco_pT"}, 
        loader_non_cheated, 
        {spill_reco_E, spill_true_E, spill_dE, spill_reco_pT},
	    cheating::cut_bad_events && cuts::reco::spill_CRTPMTNeutrino
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
            "leading_proton_purity", 
            "leading_proton_completeness", 
            "muon_L_reco_true_ratio",
            "leading_proton_L_reco_true_ratio",
            "muon_P_reco_true_ratio",
            "leading_proton_P_reco_true_ratio",
            "muon_true_L", 
            "leading_proton_true_L", 
            "muon_R",
            "leading_proton_R",
            "muon_rangeP",
            "leading_proton_rangeP",
            "CT3D_muon_leading_proton_from_rangeP",
            "CT3D_muon_leading_proton_from_trueP"
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
            spill_pid_muon_L_reco_true_ratio,
            spill_pid_proton_L_reco_true_ratio,
            spill_pid_muon_P_reco_true_ratio,
            spill_pid_proton_P_reco_true_ratio,
            spill_pid_muon_true_L,  
            spill_pid_proton_true_L, 
            spill_pid_muon_R,
            spill_pid_proton_R,
            spill_muon_momentum_rangeP,
            spill_proton_momentum_rangeP,
            spill_CT3D_rangeP,
            spill_CT3D_trueP
        }, 
        cheating::cut_bad_events && cuts::reco::spill_CRTPMTNeutrino
    )); 

    std::unique_ptr<ana::Tree> non_cheated_nopid_checks(new ana::Tree(
        "non_cheated_nopid_checks",
        {
            "vertex_difference_3D",
            "vertex_difference_x",
            "vertex_difference_y",
            "vertex_difference_z",
            "muon_hit_purity",
            "leading_proton_hit_purity",
            "muon_hit_completeness",
            "leading_proton_hit_completeness",
            "muon_L_reco_true_ratio",
            "leading_proton_L_reco_true_ratio",
            "muon_P_reco_true_ratio",
            "leading_proton_P_reco_true_ratio",
            "neutrino_reco_E",
            "neutrino_true_E",
            "neutrino_reco_dE",
            "neutrino_reco_pT",
            "muon_R",
            "leading_proton_R"
        }, 
        loader_non_cheated,
        {
            nopid_spill_vertex_difference_3D,
            nopid_spill_vertex_difference_x,
            nopid_spill_vertex_difference_y,
            nopid_spill_vertex_difference_z,
            nopid_spill_muon_hit_purity,
            nopid_spill_leading_proton_hit_purity,
            nopid_spill_muon_hit_completeness,
            nopid_spill_leading_proton_hit_completeness,
            nopid_spill_muon_L_reco_true_ratio,
            nopid_spill_leading_proton_L_reco_true_ratio,
            nopid_spill_muon_P_reco_true_ratio,
            nopid_spill_leading_proton_P_reco_true_ratio,
            nopid_spill_neutrino_reco_E,
            nopid_spill_neutrino_true_E,
            nopid_spill_neutrino_reco_dE,
            nopid_spill_neutrino_reco_pT,
            nopid_spill_muon_R,
            nopid_spill_leading_proton_R
        },
        cheating::cut_bad_events && cuts::truth::spill_1mu1p_MC && cuts::truth::spill_sliceless
    ));

    std::unique_ptr<ana::Tree> non_cheated_nopid_checks_multip(new ana::Tree(
        "non_cheated_nopid_checks_multip",
        {
            "protons_hit_purity",
            "protons_hit_completeness",
            "protons_L_reco_true_ratio",
            "protons_P_reco_true_ratio",
            "protons_R"
        }, 
        loader_non_cheated,
        {
            nopid_spill_protons_hit_purity,
            nopid_spill_protons_hit_completeness,
            nopid_spill_protons_L_reco_true_ratio,
            nopid_spill_protons_P_reco_true_ratio,
            nopid_spill_protons_R
        },
        cheating::cut_bad_events && cuts::truth::spill_1mu1p_MC && cuts::truth::spill_sliceless
    ));
    loader_non_cheated.Go();

    std::unique_ptr<TFile> file_reco_1muNp(new TFile("reco1muNp_checks.root", "RECREATE"));
    std::cout << "__ WRITING :     cheated __ " << std::endl;
    cheated->SaveTo(file_reco_1muNp->mkdir("cheated"));
    cheated_checks->SaveTo(file_reco_1muNp->GetDirectory("cheated"));
    cheated_nopid_checks->SaveTo(file_reco_1muNp->mkdir("cheated_nopid"));
    cheated_nopid_checks_multip->SaveTo(file_reco_1muNp->GetDirectory("cheated_nopid"));
    std::cout << "__ WRITING : non cheated __ " << std::endl;
    non_cheated->SaveTo(file_reco_1muNp->mkdir("non_cheated"));
    non_cheated_checks->SaveTo(file_reco_1muNp->GetDirectory("non_cheated"));
    non_cheated_nopid_checks->SaveTo(file_reco_1muNp->mkdir("non_cheated_nopid"));
    non_cheated_nopid_checks_multip->SaveTo(file_reco_1muNp->GetDirectory("non_cheated_nopid"));
};

#endif 
