#pragma once 

#include "sbnana/CAFAna/Core/Tree.h"
#include "pdg_helper.h"

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


const ana::SpillVar reco_true_vtxPull_z = SPILLVAR(-9999, vars::reco::slice_vertex_difference_z, def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_vtxPull_y = SPILLVAR(-9999, vars::reco::slice_vertex_difference_y, def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_vtxPull_x = SPILLVAR(-9999, vars::reco::slice_vertex_difference_x, def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_vtxPull   = SPILLVAR(-9999, vars::reco::slice_vertex_difference,   def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar true_vtxPull_z      = SPILLVAR(-9999, vars::reco::slice_vertex_difference_z, def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_vtxPull_y      = SPILLVAR(-9999, vars::reco::slice_vertex_difference_y, def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_vtxPull_x      = SPILLVAR(-9999, vars::reco::slice_vertex_difference_x, def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_vtxPull        = SPILLVAR(-9999, vars::reco::slice_vertex_difference,   def_cut, true_1uNp, def_cut_truth);



const ana::SpillVar reco_true_muon_vtxDist          = SPILLVAR(-9999, vars::pdg::slice_muon_vtxDist,        def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muonPandoraPrimary    = SPILLVAR(-9999, vars::pdg::slice_muonPandoraPrimary,  def_cut, reco_true_1uNp, def_cut_truth);

const ana::SpillVar true_muon_vtxDist               = SPILLVAR(-9999, vars::pdg::slice_muon_vtxDist,        def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muonPandoraPrimary         = SPILLVAR(-9999, vars::pdg::slice_muonPandoraPrimary,  def_cut, true_1uNp, def_cut_truth);

const ana::SpillMultiVar reco_true_proton_vtxDist       = SPILLMULTIVAR({}, vars::pdg::slice_proton_vtxDist,        def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_protonPandoraPrimary = SPILLMULTIVAR({}, vars::pdg::slice_protonPandoraPrimary,  def_cut, reco_true_1uNp, def_cut_truth);

const ana::SpillMultiVar true_proton_vtxDist            = SPILLMULTIVAR({}, vars::pdg::slice_proton_vtxDist,        def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_protonPandoraPrimary      = SPILLMULTIVAR({}, vars::pdg::slice_protonPandoraPrimary,  def_cut, true_1uNp, def_cut_truth);

void CCNp_vertexStudy () {

    std::unique_ptr<TFile> file_1muNp(new TFile("CCNp_vertexStudy.root", "RECREATE"));

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

    std::vector<std::unique_ptr<ana::Tree>> trees_muons, trees_protons, trees_vertex;

    for (auto const& running_loader: running_loaders) {
        
        ana::SpectrumLoader loader (loaders_available.at(running_loader).c_str());

        trees_muons.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(),
            std::vector<std::string>{
                "vtxDist",
                "pandoraPrimary"
            },
            loader,
            std::vector<ana::SpillVar>{
                reco_true_muon_vtxDist,
                reco_true_muonPandoraPrimary
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees_muons.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(),
            std::vector<std::string>{
                "vtxDist",
                "pandoraPrimary"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                true_muon_vtxDist,
                true_muonPandoraPrimary
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees_protons.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(),
            std::vector<std::string>{
                "vtxDist",
                "pandoraPrimary"
            },
            loader,
            std::vector<ana::SpillMultiVar>{
                reco_true_proton_vtxDist,
                reco_true_protonPandoraPrimary
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees_protons.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(),
            std::vector<std::string>{
                "vtxDist",
                "pandoraPrimary"
            }, 
            loader,
            std::vector<ana::SpillMultiVar>{
                true_proton_vtxDist,
                true_protonPandoraPrimary
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees_vertex.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(),
            std::vector<std::string>{
                "vtxDist_z",
                "vtxDist_y",
                "vtxDist_x",
                "vtxDist"
            },
            loader,
            std::vector<ana::SpillVar>{
                reco_true_vtxPull_z,
                reco_true_vtxPull_y,
                reco_true_vtxPull_x,
                reco_true_vtxPull
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees_vertex.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(),
            std::vector<std::string>{
                "vtxDist_z",
                "vtxDist_y",
                "vtxDist_x",
                "vtxDist"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                true_vtxPull_z,
                true_vtxPull_y,
                true_vtxPull_x,
                true_vtxPull
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        loader.Go();

    }

    file_1muNp->mkdir("muons");
    for (auto const& tree: trees_muons) 
        tree->SaveTo(file_1muNp->GetDirectory("muons"));

    file_1muNp->mkdir("protons");
    for (auto const& tree: trees_protons) 
        tree->SaveTo(file_1muNp->GetDirectory("protons"));

    file_1muNp->mkdir("vertex");
    for (auto const& tree: trees_vertex) 
        tree->SaveTo(file_1muNp->GetDirectory("vertex"));

}