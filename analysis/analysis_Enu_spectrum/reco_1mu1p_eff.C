
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

#include "helper.h"
#include "selection.h"

#include "sbnana/CAFAna/Core/Tree.h"
#include <map>

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

const ana::Cut def_cut_truth = (
    cuts::truth::slice_vtx_in_FV
);

const cut_type_t reco_true_1u1p = cut_type_t::BOTH_1mu1p;
const cut_type_t true_1u1p = cut_type_t::TRUE_1mu1p;
const cut_type_t reco_1u1p = cut_type_t::RECO;
const particle_data::int_type_t type = particle_data::int_type_t::true_visible_1mu1p;

// true + reco
const ana::SpillVar spill_reco_E_reco_true_cut  = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, reco_true_1u1p, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut  = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_true_1u1p, def_cut_truth);
const ana::SpillVar spill_reco_pT_reco_true_cut = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, reco_true_1u1p, def_cut_truth);

// true only
const ana::SpillVar spill_reco_E_true_cut       = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, true_1u1p, def_cut_truth);
const ana::SpillVar spill_true_E_true_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, true_1u1p, def_cut_truth);
const ana::SpillVar spill_reco_pT_true_cut      = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, true_1u1p, def_cut_truth);

// reco only
const ana::SpillVar spill_reco_E_reco_cut       = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, reco_1u1p, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_1u1p, def_cut_truth);
const ana::SpillVar spill_reco_pT_reco_cut      = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, reco_1u1p, def_cut_truth);

void reco_1mu1p_eff() {

    ana::SpectrumLoader nominal_loader                  ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_nominal");
    ana::SpectrumLoader cheated_Mva_loader              ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_mva");
    ana::SpectrumLoader cheated_Vtx_loader              ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_vtx");
    ana::SpectrumLoader cheated_VtxSelection_loader     ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_vtx_selection");    
    ana::SpectrumLoader cheated_2D_loader               ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D");
    ana::SpectrumLoader cheated_2D_Vtx_loader           ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D_vtx");
    ana::SpectrumLoader cheated_2D_Vtx_3D_loader        ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D_vtx_3D");
    ana::SpectrumLoader cheated_2D_Vtx_3D_Nu_loader     ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D_vtx_3D_Nu");
    ana::SpectrumLoader cheated_2D_Vtx_3D_Nu_Mva_loader ("msotgia_v09_89_01_01p03_down_vtx_both_ifdh_reco_cheated_2D_vtx_3D_Nu_mva");

    std::map<std::string, ana::SpectrumLoader*> loaders_available = {
        {"nominal_loader",                  &nominal_loader},
        {"cheated_Mva_loader",              &cheated_Mva_loader},
        {"cheated_Vtx_loader",              &cheated_Vtx_loader},
        {"cheated_VtxSelection_loader",     &cheated_VtxSelection_loader},
        {"cheated_2D_loader",               &cheated_2D_loader},
        {"cheated_2D_Vtx_loader",           &cheated_2D_Vtx_loader},
        {"cheated_2D_Vtx_3D_loader",        &cheated_2D_Vtx_3D_loader},
        {"cheated_2D_Vtx_3D_Nu_loader",     &cheated_2D_Vtx_3D_Nu_loader},
        {"cheated_2D_Vtx_3D_Nu_Mva_loader", &cheated_2D_Vtx_3D_Nu_Mva_loader}
    };

    // Running all :)
    std::vector<std::string> running_loaders = {
        "nominal_loader",
        "cheated_Mva_loader",
        "cheated_Vtx_loader",
        "cheated_VtxSelection_loader",
        "cheated_2D_loader",
        "cheated_2D_Vtx_loader",
        "cheated_2D_Vtx_3D_loader",
        "cheated_2D_Vtx_3D_Nu_loader",
        "cheated_2D_Vtx_3D_Nu_Mva_loader"
    };

    std::vector<std::unique_ptr<ana::Tree>> trees;

    for (auto const& running_loader: running_loaders) {
        
        ana::SpectrumLoader& loader = *loaders_available.at(running_loader);

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(), 
            std::vector<std::string>{"event", "reco_E", "true_E", "reco_pT"}, 
            loader,
            std::vector<ana::SpillVar>{event, spill_reco_E_reco_true_cut, spill_true_E_reco_true_cut, spill_reco_pT_reco_true_cut},
            cuts::reco::spill_CRTPMTNeutrino
        ));
        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_" + running_loader).c_str(), 
            std::vector<std::string>{"event", "reco_E", "true_E", "reco_pT"}, 
            loader,
            std::vector<ana::SpillVar>{event, spill_reco_E_reco_cut, spill_true_E_reco_cut, spill_reco_pT_reco_cut},
            cuts::reco::spill_CRTPMTNeutrino
        ));
        trees.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(), 
            std::vector<std::string>{"event", "reco_E", "true_E", "reco_pT"}, 
            loader,
            std::vector<ana::SpillVar>{event, spill_reco_E_true_cut, spill_true_E_true_cut, spill_reco_pT_true_cut},
            cuts::reco::spill_CRTPMTNeutrino
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1mu1p(new TFile("2k_efficiency_plot_1u1p.root", "RECREATE"));
    file_1mu1p->mkdir("efficiency_studies");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1mu1p->GetDirectory("efficiency_studies"));
};

#endif 
