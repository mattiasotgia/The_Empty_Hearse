
#pragma once

#include "helper.h"
#include "selection.h"
#include "slice_helper.h"

#include "sbnana/CAFAna/Core/Tree.h"
#include <map>

#define SPILLVAR(_def, _var, _reco, _what, _true)                                                                               \
    var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(_def, _var, _reco, _what, _true)
#define SPILLMULTIVAR(_def, _var, _reco, _what, _true)                                                                          \
    var_utils::make_spill_from_slice<ana::SpillMultiVar, ana::MultiVar, std::vector<double>>(_def, _var, _reco, _what, _true)

using level_t = logger::level;
using cut_type_t = var_utils::cut_type_t;

const ana::Cut def_cut = (
    (cuts::reco::slice_1mu1p || cuts::reco::slice_1muNp) &&
    cuts::reco::slice_at_least_mu       &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_vtx_in_FV
);

const cut_type_t reco_true_1uNp = cut_type_t::BOTH_1muNp;
const cut_type_t true_1uNp = cut_type_t::TRUE_1muNp;
const cut_type_t reco_1uNp = cut_type_t::RECO;

// true + reco
const ana::SpillVar spill_reco_E_reco_true_cut  = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut  = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_reco_pT_reco_true_cut = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_slcEff_reco_true_cut  = SPILLVAR(-9999, vars::slice::slice_efficiency,            def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_slcPur_reco_true_cut  = SPILLVAR(-9999, vars::slice::slice_purity,                def_cut, reco_true_1uNp, def_cut_truth);

// true only
const ana::SpillVar spill_reco_E_true_cut       = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_true_E_true_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_reco_pT_true_cut      = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar spill_slcEff_true_cut       = SPILLVAR(-9999, vars::slice::slice_efficiency,            def_cut, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_slcPur_true_cut       = SPILLVAR(-9999, vars::slice::slice_purity,                def_cut, true_1uNp, def_cut_truth); 

// reco only
const ana::SpillVar spill_reco_E_reco_cut       = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_reco_pT_reco_cut      = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_slcEff_reco_cut       = SPILLVAR(-9999, vars::slice::slice_efficiency,            def_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_slcPur_reco_cut       = SPILLVAR(-9999, vars::slice::slice_purity,                def_cut, reco_1uNp, def_cut_truth);

void CCNp_efficiency() {

    std::map<std::string, std::string> loaders_available = {
        {"nominal",                  "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_nominal"},
        {"cheated_2d",               "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d"},
        {"cheated_vtx",              "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_vtx"},
        {"cheated_vtxSelection",     "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_vtxSelection"},
        {"cheated_3d",               "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_3d"},
        {"cheated_nuH",              "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_nu"},
        {"cheated_mva",              "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_mva"},
        {"cheated_2d_vtx",           "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d_vtx"},
        {"cheated_2d_vtx_3d",        "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d_vtx_3d"},
        {"cheated_2d_vtx_3d_nu",     "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d_vtx_3d_nu"},
        {"cheated_2d_vtx_3d_nu_mva", "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d_vtx_3d_nu_mva"}
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
        
        ana::SpectrumLoader loader (loaders_available.at(running_loader).c_str());

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(), 
            std::vector<std::string>{"event", "reco_E", "true_E", "reco_pT", "sliceEfficiency", "slicePurity"}, 
            loader,
            std::vector<ana::SpillVar>{event, spill_reco_E_reco_true_cut, spill_true_E_reco_true_cut, spill_reco_pT_reco_true_cut, spill_slcEff_reco_true_cut,spill_slcPur_reco_true_cut},
            cuts::reco::spill_CRTPMTNeutrino
        ));
        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_" + running_loader).c_str(), 
            std::vector<std::string>{"event", "reco_E", "true_E", "reco_pT", "sliceEfficiency", "slicePurity"}, 
            loader,
            std::vector<ana::SpillVar>{event, spill_reco_E_reco_cut, spill_true_E_reco_cut, spill_reco_pT_reco_cut, spill_slcEff_reco_cut, spill_slcPur_reco_cut},
            cuts::reco::spill_CRTPMTNeutrino
        ));
        trees.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(), 
            std::vector<std::string>{"event", "reco_E", "true_E", "reco_pT", "sliceEfficiency", "slicePurity"}, 
            loader,
            std::vector<ana::SpillVar>{event, spill_reco_E_true_cut, spill_true_E_true_cut, spill_reco_pT_true_cut, spill_slcEff_true_cut, spill_slcPur_true_cut},
            cuts::reco::spill_CRTPMTNeutrino
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1muNp(new TFile("CCNp_efficiency.root", "RECREATE"));
    file_1muNp->mkdir("efficiency_studies");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1muNp->GetDirectory("efficiency_studies"));
};

