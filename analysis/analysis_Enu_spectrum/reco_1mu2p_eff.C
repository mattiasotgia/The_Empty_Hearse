
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
    cuts::reco::slice_1mu2p             &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_vtx_in_FV
);

const cut_type_t reco_true_1u2p = cut_type_t::BOTH_1mu2p;
const cut_type_t true_1u2p = cut_type_t::TRUE_1mu2p;
const cut_type_t reco_1u2p = cut_type_t::RECO;
const particle_data::int_type_t type = particle_data::int_type_t::true_visible_1mu2p;

// true + reco
const ana::SpillVar spill_reco_E_reco_true_cut  = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, reco_true_1u2p, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut  = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_true_1u2p, def_cut_truth);
const ana::SpillVar spill_reco_pT_reco_true_cut = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, reco_true_1u2p, def_cut_truth);

// true only
const ana::SpillVar spill_reco_E_true_cut       = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, true_1u2p, def_cut_truth);
const ana::SpillVar spill_true_E_true_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, true_1u2p, def_cut_truth);
const ana::SpillVar spill_reco_pT_true_cut      = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, true_1u2p, def_cut_truth);

// reco only
const ana::SpillVar spill_reco_E_reco_cut       = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, reco_1u2p, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_1u2p, def_cut_truth);
const ana::SpillVar spill_reco_pT_reco_cut      = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, reco_1u2p, def_cut_truth);

void reco_1mu2p_eff() {

    // Older sample run...
    // ana::SpectrumLoader cheated_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_stage1tocaf_cheated_flatcaf");
    // ana::SpectrumLoader nominal_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_stage1tocaf_nominal_flatcaf");

    // Newer sample run...
    ana::SpectrumLoader cheated_2D_Vtx_3D_Nu_Mva_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_down_ladder_cheated_2D_Vtx_3D_Nu_Mva");
    ana::SpectrumLoader cheated_2D_Vtx_3D_Nu_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_down_ladder_cheated_2D_Vtx_3D_Nu");
    ana::SpectrumLoader cheated_2D_Vtx_3D_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_down_ladder_cheated_2D_Vtx_3D");
    ana::SpectrumLoader cheated_2D_Vtx_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_down_ladder_cheated_2D_Vtx");
    ana::SpectrumLoader cheated_2D_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_down_ladder_cheated_2D");
    ana::SpectrumLoader nominal_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_down_ladder_nominal");


    std::map<std::string, ana::SpectrumLoader*> loaders_available = {
        {"cheated_2D_Vtx_3D_Nu_Mva",    &cheated_2D_Vtx_3D_Nu_Mva_loader},
        {"cheated_2D_Vtx_3D_Nu",        &cheated_2D_Vtx_3D_Nu_loader},
        {"cheated_2D_Vtx_3D",           &cheated_2D_Vtx_3D_loader},
        {"cheated_2D_Vtx",              &cheated_2D_Vtx_loader},
        {"cheated_2D",                  &cheated_2D_loader},
        {"nominal_reconstruction",      &nominal_loader}
    };

    // Running all :)
    std::vector<std::string> running_loaders = {
        "cheated_2D_Vtx_3D_Nu_Mva",
        "cheated_2D_Vtx_3D_Nu",
        "cheated_2D_Vtx_3D",
        "cheated_2D_Vtx",
        "cheated_2D",
        "nominal_reconstruction"
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

    std::unique_ptr<TFile> file_1mu2p(new TFile("efficiency_plot_1u2p.root", "RECREATE"));
    file_1mu2p->mkdir("efficiency_studies");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1mu2p->GetDirectory("efficiency_studies"));
};

#endif 
