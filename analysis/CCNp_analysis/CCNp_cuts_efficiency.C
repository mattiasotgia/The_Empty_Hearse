
#pragma once

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
    (cuts::reco::slice_1mu1p || cuts::reco::slice_1muNp) &&
    cuts::reco::slice_at_least_mu       &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_vtx_in_FV
);

const ana::Cut s1 = (ana::kNoCut);
const ana::Cut s2 = (cuts::reco::slice_barycenter);
const ana::Cut s3 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV);
const ana::Cut s4 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV && cuts::reco::slice_all_trk_contained);
const ana::Cut s5 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV && cuts::reco::slice_all_trk_contained && cuts::reco::slice_at_least_mu);
const ana::Cut s6 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV && cuts::reco::slice_all_trk_contained && cuts::reco::slice_at_least_mu && cuts::reco::slice_Np_only);
const ana::Cut s7 = (cuts::reco::slice_barycenter && cuts::reco::slice_vtx_in_FV && cuts::reco::slice_all_trk_contained && cuts::reco::slice_at_least_mu && cuts::reco::slice_Np_only && cuts::reco::slice_0pi_only && cuts::reco::slice_0showers_only);

const cut_type_t reco_true_1uNp = cut_type_t::BOTH_1muNp;
const cut_type_t true_1uNp = cut_type_t::TRUE_1muNp;
const cut_type_t reco_1uNp = cut_type_t::RECO;


// true + reco
const ana::SpillVar spill_true_E_reco_true_cut_s1 = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s1, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut_s2 = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s2, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut_s3 = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s3, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut_s4 = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s4, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut_s5 = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s5, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut_s6 = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s6, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut_s7 = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s7, reco_true_1uNp, def_cut_truth);

// true only
const ana::SpillVar spill_true_E_true_cut_s1      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s1, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_true_E_true_cut_s2      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s2, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_true_E_true_cut_s3      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s3, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_true_E_true_cut_s4      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s4, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_true_E_true_cut_s5      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s5, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_true_E_true_cut_s6      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s6, true_1uNp, def_cut_truth); 
const ana::SpillVar spill_true_E_true_cut_s7      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s7, true_1uNp, def_cut_truth); 

// reco only
const ana::SpillVar spill_true_E_reco_cut_s1      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s1, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut_s2      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s2, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut_s3      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s3, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut_s4      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s4, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut_s5      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s5, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut_s6      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s6, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut_s7      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, s7, reco_1uNp, def_cut_truth);

void CCNp_cuts_efficiency() {

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
        
        ana::SpectrumLoader& loader = *loaders_available.at(running_loader);

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(), 
            std::vector<std::string>{
                "event",
                "s1",
                "s2",
                "s3",
                "s4",
                "s5",
                "s6",
                "s7"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                event,
                spill_true_E_reco_true_cut_s1,
                spill_true_E_reco_true_cut_s2,
                spill_true_E_reco_true_cut_s3,
                spill_true_E_reco_true_cut_s4,
                spill_true_E_reco_true_cut_s5,
                spill_true_E_reco_true_cut_s6,
                spill_true_E_reco_true_cut_s7
            },
            cuts::reco::spill_CRTPMTNeutrino
        ));
        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_" + running_loader).c_str(), 
            std::vector<std::string>{
                "event",
                "s1",
                "s2",
                "s3",
                "s4",
                "s5",
                "s6",
                "s7"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                event,
                spill_true_E_reco_cut_s1,
                spill_true_E_reco_cut_s2,
                spill_true_E_reco_cut_s3,
                spill_true_E_reco_cut_s4,
                spill_true_E_reco_cut_s5,
                spill_true_E_reco_cut_s6,
                spill_true_E_reco_cut_s7
            },
            cuts::reco::spill_CRTPMTNeutrino
        ));
        trees.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(), 
            std::vector<std::string>{
                "event",
                "s1",
                "s2",
                "s3",
                "s4",
                "s5",
                "s6",
                "s7"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                event,
                spill_true_E_true_cut_s1,
                spill_true_E_true_cut_s2,
                spill_true_E_true_cut_s3,
                spill_true_E_true_cut_s4,
                spill_true_E_true_cut_s5,
                spill_true_E_true_cut_s6,
                spill_true_E_true_cut_s7
            },
            cuts::reco::spill_CRTPMTNeutrino
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1muNp(new TFile("CCNp_efficiency_byCut.root", "RECREATE"));
    file_1muNp->mkdir("cuts");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1muNp->GetDirectory("cuts"));
};

