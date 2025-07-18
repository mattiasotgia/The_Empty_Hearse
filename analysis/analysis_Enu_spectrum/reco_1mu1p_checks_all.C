#pragma once

#include "helper.h"
#include "selection.h"

#include "sbnana/CAFAna/Core/Tree.h"

#define SPILLVAR(_def, _var, _reco, _what, _true)                                                                               \
    var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(_def, _var, _reco, _what, _true)
#define SPILLMULTIVAR(_def, _var, _reco, _what, _true)                                                                          \
    var_utils::make_spill_from_slice<ana::SpillMultiVar, ana::MultiVar, std::vector<double>>(_def, _var, _reco, _what, _true)


using level_t = logger::level;
using cut_type_t = var_utils::cut_type_t;

const ana::Cut def_cut = (
    cuts::reco::slice_1mu1p             &&
    cuts::reco::slice_at_least_mu       &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut s8 = (
    cuts::reco::slice_1mu1p             &&
    cuts::reco::slice_at_least_mu       &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut s1 = (ana::kNoCut);
const ana::Cut s2 = (cuts::reco::slice_barycenter);
const ana::Cut s3 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV);
const ana::Cut s4 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV && cuts::reco::slice_all_trk_contained);
const ana::Cut s5 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV && cuts::reco::slice_all_trk_contained && cuts::reco::slice_at_least_mu);
const ana::Cut s6 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV && cuts::reco::slice_all_trk_contained && cuts::reco::slice_at_least_mu && cuts::reco::slice_1p_only);
const ana::Cut s7 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV && cuts::reco::slice_all_trk_contained && cuts::reco::slice_at_least_mu && cuts::reco::slice_1p_only && cuts::reco::slice_0pi_only && cuts::reco::slice_0showers_only);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_vtx_in_FV
);

const cut_type_t reco_true_1u1p = cut_type_t::BOTH_1mu1p;
const cut_type_t true_1u1p = cut_type_t::TRUE_1mu1p;
const cut_type_t reco_1u1p = cut_type_t::RECO;
const particle_data::int_type_t type = particle_data::int_type_t::true_visible_1mu1p;

const cut_type_t local_cut_type = reco_true_1u1p;

const ana::SpillVar spill_dE                            = SPILLVAR(-9999, vars::truth::slice_neutrino_dE,                   def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_reco_E                        = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,          def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_true_E                        = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,               def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_reco_pT                       = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,              def_cut, local_cut_type, def_cut_truth);
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
const ana::SpillVar spill_pid_muon_L_ratio_proxy        = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_ratio_proxy,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_true_L             = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_reco_L             = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_muon_R                    = SPILLVAR(-9999, vars::reco::slice_muon_R,                         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_pid_proton_R                  = SPILLVAR(-9999, vars::reco::slice_leading_proton_R,               def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_muon_momentum_rangeP          = SPILLVAR(-9999, vars::reco::slice_muon_momentum_rangeP,           def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_proton_momentum_rangeP        = SPILLVAR(-9999, vars::reco::slice_proton_momentum_rangeP,         def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_CT3D_rangeP                   = SPILLVAR(-9999, vars::reco::slice_CT3D_rangeP_muon_proton,        def_cut, local_cut_type, def_cut_truth);
const ana::SpillVar spill_CT3D_trueP                    = SPILLVAR(-9999, vars::reco::slice_CT3D_trueP_muon_proton,         def_cut, local_cut_type, def_cut_truth);


void reco_1mu1p_checks_all() {

    ana::SpectrumLoader nominal_loader                     ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_nominal");
    ana::SpectrumLoader cheated_2d_loader                  ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_2d");
    ana::SpectrumLoader cheated_vtx_loader                 ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_vtx");
    ana::SpectrumLoader cheated_vtxSelection_loader        ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_vtxSelection");
    ana::SpectrumLoader cheated_3d_loader                  ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_3d");
    ana::SpectrumLoader cheated_nuH_loader                 ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_nuH");
    ana::SpectrumLoader cheated_mva_loader                 ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_mva");
    ana::SpectrumLoader cheated_2d_vtx_loader              ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_2d_vtx");
    ana::SpectrumLoader cheated_2d_vtx_3d_loader           ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_2d_vtx_3d");
    ana::SpectrumLoader cheated_2d_vtx_3d_nu_loader        ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_2d_vtx_3d_nu");
    ana::SpectrumLoader cheated_2d_vtx_3d_nu_mva_loader    ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_2d_vtx_3d_nu_mva");

    // ana::SpectrumLoader nominal_loader                  ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_nominal");
    // ana::SpectrumLoader cheated_Mva_loader              ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_mva");
    // ana::SpectrumLoader cheated_Vtx_loader              ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_vtx");
    // ana::SpectrumLoader cheated_VtxSelection_loader     ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_vtx_selection");    
    // ana::SpectrumLoader cheated_2D_loader               ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D");
    // ana::SpectrumLoader cheated_2D_Vtx_loader           ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D_vtx");
    // ana::SpectrumLoader cheated_2D_Vtx_3D_loader        ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D_vtx_3D");
    // ana::SpectrumLoader cheated_2D_Vtx_3D_Nu_loader     ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D_vtx_3D_Nu");
    // ana::SpectrumLoader cheated_2D_Vtx_3D_Nu_Mva_loader ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D_vtx_3D_Nu_mva");

    std::map<std::string, ana::SpectrumLoader*> loaders_available = {
        {"nominal",                  &nominal_loader},
        {"cheated_2d",               &cheated_2d_loader},
        {"cheated_vtx",              &cheated_vtx_loader},
        {"cheated_vtxSelection",     &cheated_vtxSelection_loader},
        {"cheated_3d",               &cheated_3d_loader},
        {"cheated_nuH",              &cheated_nuH_loader},
        {"cheated_mva",              &cheated_mva_loader},
        {"cheated_2d_vtx",           &cheated_2d_vtx_loader},
        {"cheated_2d_vtx_3d",        &cheated_2d_vtx_3d_loader},
        {"cheated_2d_vtx_3d_nu",     &cheated_2d_vtx_3d_nu_loader},
        {"cheated_2d_vtx_3d_nu_mva", &cheated_2d_vtx_3d_nu_mva_loader}
    };

    // Running all :)
    std::vector<std::string> running_loaders = {
        "nominal",
        "cheated_2d",
        "cheated_vtx",
        "cheated_vtxSelection",
        "cheated_3d",
        "cheated_nuH",
        "cheated_mva",
        "cheated_2d_vtx",
        "cheated_2d_vtx_3d",
        "cheated_2d_vtx_3d_nu",
        "cheated_2d_vtx_3d_nu_mva"
    };

    std::vector<std::unique_ptr<ana::Tree>> trees;

    for (auto const& running_loader: running_loaders) {

        ana::SpectrumLoader* loader = nullptr;
        try {
            loader = loaders_available.at(running_loader);
            // Proceed with using `loader`
        } catch (const std::out_of_range& e) {
            logger::log(level_t::warning) << "Error: Loader '" << running_loader << "' not found in loaders_available.\n";
            logger::log(level_t::warning) << "Exception: " << e.what() << '\n';
            // Handle the error (e.g., fallback, abort, log, etc.)
            continue;
        }
        // ana::SpectrumLoader& loader = *loaders_available.at(running_loader);

        trees.emplace_back(std::make_unique<ana::Tree>(
            running_loader.c_str(),
            std::vector<std::string>{
                "event",
                "dE",
                "reco_E",
                "true_E",
                "reco_pT",
                "vertex_difference",
                "vertex_difference_x",
                "vertex_difference_y",
                "vertex_difference_z",
                "muon_purity",
                "muon_completeness",
                "proton_purity",
                "proton_completeness",
                "pid_muon_L_reco_true_ratio",
                "pid_proton_L_reco_true_ratio",
                "pid_muon_P_reco_true_ratio",
                "pid_proton_P_reco_true_ratio",
                "pid_muon_true_L",
                "pid_muon_reco_L",
                "epid_muon_L_ratio_proxy",
                "pid_proton_true_L",
                "pid_proton_reco_L",
                "pid_muon_R",
                "pid_proton_R",
                "muon_momentum_rangeP",
                "proton_momentum_rangeP",
                "CT3D_rangeP",
                "CT3D_trueP"
            }, 
            *loader,
            std::vector<ana::SpillVar>{
                event,
                spill_dE,
                spill_reco_E,
                spill_true_E,
                spill_reco_pT,
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
                spill_pid_muon_reco_L,
                spill_pid_muon_L_ratio_proxy,
                spill_pid_proton_true_L,
                spill_pid_proton_reco_L,
                spill_pid_muon_R,
                spill_pid_proton_R,
                spill_muon_momentum_rangeP,
                spill_proton_momentum_rangeP,
                spill_CT3D_rangeP,
                spill_CT3D_trueP
            },
            cuts::reco::spill_CRTPMTNeutrino
        ));


        loader->Go();
    }

    std::unique_ptr<TFile> file_1mu1p(new TFile("2k_v2_variables_1u1p.root", "RECREATE"));
    file_1mu1p->mkdir("tests");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1mu1p->GetDirectory("tests"));


    return;
}
