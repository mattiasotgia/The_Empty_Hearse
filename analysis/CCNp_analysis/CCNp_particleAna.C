#pragma once 

#include "sbnana/CAFAna/Core/Tree.h"
#include "pdg_helper.h"
#include "slice_helper.h"

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

const ana::SpillVar reco_true_muon_chi2_mu          = SPILLVAR(-9999, vars::pdg::slice_muon_chi2_mu,        def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muon_chi2_proton      = SPILLVAR(-9999, vars::pdg::slice_muon_chi2_p,         def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muon_length           = SPILLVAR(-9999, vars::pdg::slice_muon_length,         def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muon_length_ratio     = SPILLVAR(-9999, vars::pdg::slice_muon_length_ratio,   def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muon_vtxDist          = SPILLVAR(-9999, vars::pdg::slice_muon_vtxDist,        def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muon_trackScore       = SPILLVAR(-9999, vars::pdg::slice_muon_trackScore,     def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muon_completeness     = SPILLVAR(-9999, vars::pdg::slice_muon_completeness,   def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muon_purity           = SPILLVAR(-9999, vars::pdg::slice_muon_purity,         def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_muonPandoraPrimary    = SPILLVAR(-9999, vars::pdg::slice_muonPandoraPrimary,  def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_slcEff                = SPILLVAR(-9999, vars::slice::slice_efficiency,        def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_slcPur                = SPILLVAR(-9999, vars::slice::slice_purity,            def_cut, reco_true_1uNp, def_cut_truth);

const ana::SpillVar true_muon_chi2_mu               = SPILLVAR(-9999, vars::pdg::slice_muon_chi2_mu,        def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muon_chi2_proton           = SPILLVAR(-9999, vars::pdg::slice_muon_chi2_p,         def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muon_length                = SPILLVAR(-9999, vars::pdg::slice_muon_length,         def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muon_length_ratio          = SPILLVAR(-9999, vars::pdg::slice_muon_length_ratio,   def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muon_vtxDist               = SPILLVAR(-9999, vars::pdg::slice_muon_vtxDist,        def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muon_trackScore            = SPILLVAR(-9999, vars::pdg::slice_muon_trackScore,     def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muon_completeness          = SPILLVAR(-9999, vars::pdg::slice_muon_completeness,   def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muon_purity                = SPILLVAR(-9999, vars::pdg::slice_muon_purity,         def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_muonPandoraPrimary         = SPILLVAR(-9999, vars::pdg::slice_muonPandoraPrimary,  def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_slcEff                     = SPILLVAR(-9999, vars::slice::slice_efficiency,        def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_slcPur                     = SPILLVAR(-9999, vars::slice::slice_purity,            def_cut, true_1uNp, def_cut_truth);

const ana::SpillMultiVar reco_true_proton_chi2_mu       = SPILLMULTIVAR({}, vars::pdg::slice_proton_chi2_mu,        def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_proton_chi2_proton   = SPILLMULTIVAR({}, vars::pdg::slice_proton_chi2_p,         def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_proton_length        = SPILLMULTIVAR({}, vars::pdg::slice_proton_length,         def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_proton_length_ratio  = SPILLMULTIVAR({}, vars::pdg::slice_proton_length_ratio,   def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_proton_vtxDist       = SPILLMULTIVAR({}, vars::pdg::slice_proton_vtxDist,        def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_proton_trackScore    = SPILLMULTIVAR({}, vars::pdg::slice_proton_trackScore,     def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_proton_completeness  = SPILLMULTIVAR({}, vars::pdg::slice_proton_completeness,   def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_proton_purity        = SPILLMULTIVAR({}, vars::pdg::slice_proton_purity,         def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_proton_depEnergy     = SPILLMULTIVAR({}, vars::pdg::slice_proton_depEnergy,      def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_protonPandoraPrimary = SPILLMULTIVAR({}, vars::pdg::slice_protonPandoraPrimary,  def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_slcEff_proton        = SPILLMULTIVAR({}, vars::slice_proton::slice_efficiency,   def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillMultiVar reco_true_slcPur_proton        = SPILLMULTIVAR({}, vars::slice_proton::slice_purity,       def_cut, reco_true_1uNp, def_cut_truth);

const ana::SpillMultiVar true_proton_chi2_mu            = SPILLMULTIVAR({}, vars::pdg::slice_proton_chi2_mu,        def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_proton_chi2_proton        = SPILLMULTIVAR({}, vars::pdg::slice_proton_chi2_p,         def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_proton_length             = SPILLMULTIVAR({}, vars::pdg::slice_proton_length,         def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_proton_length_ratio       = SPILLMULTIVAR({}, vars::pdg::slice_proton_length_ratio,   def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_proton_vtxDist            = SPILLMULTIVAR({}, vars::pdg::slice_proton_vtxDist,        def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_proton_trackScore         = SPILLMULTIVAR({}, vars::pdg::slice_proton_trackScore,     def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_proton_completeness       = SPILLMULTIVAR({}, vars::pdg::slice_proton_completeness,   def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_proton_purity             = SPILLMULTIVAR({}, vars::pdg::slice_proton_purity,         def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_proton_depEnergy          = SPILLMULTIVAR({}, vars::pdg::slice_proton_depEnergy,      def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_protonPandoraPrimary      = SPILLMULTIVAR({}, vars::pdg::slice_protonPandoraPrimary,  def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_slcEff_proton             = SPILLMULTIVAR({}, vars::slice_proton::slice_efficiency,   def_cut, true_1uNp, def_cut_truth);
const ana::SpillMultiVar true_slcPur_proton             = SPILLMULTIVAR({}, vars::slice_proton::slice_purity,       def_cut, true_1uNp, def_cut_truth);

void CCNp_particleAna () {

    std::unique_ptr<TFile> file_1muNp(new TFile("CCNp_particleAna.root", "RECREATE"));

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

    std::vector<std::unique_ptr<ana::Tree>> trees_muons, trees_protons;

    for (auto const& running_loader: running_loaders) {
        
        ana::SpectrumLoader loader (loaders_available.at(running_loader).c_str());

        trees_muons.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(),
            std::vector<std::string>{
                "event",
                "chi2_mu",
                "chi2_proton",
                "length",
                "length_ratio",
                "vtxDist",
                "trackScore",
                "completeness",
                "purity",
                "pandoraPrimary",
                "sliceEfficiency", 
                "slicePurity"
            },
            loader,
            std::vector<ana::SpillVar>{
                event,
                reco_true_muon_chi2_mu,
                reco_true_muon_chi2_proton,
                reco_true_muon_length,
                reco_true_muon_length_ratio,
                reco_true_muon_vtxDist,
                reco_true_muon_trackScore,
                reco_true_muon_completeness,
                reco_true_muon_purity,
                reco_true_muonPandoraPrimary,
                reco_true_slcEff,
                reco_true_slcPur
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees_muons.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(),
            std::vector<std::string>{
                "event",
                "chi2_mu",
                "chi2_proton",
                "length",
                "length_ratio",
                "vtxDist",
                "trackScore",
                "completeness",
                "purity",
                "pandoraPrimary",
                "sliceEfficiency", 
                "slicePurity"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                event,
                true_muon_chi2_mu,
                true_muon_chi2_proton,
                true_muon_length,
                true_muon_length_ratio,
                true_muon_vtxDist,
                true_muon_trackScore,
                true_muon_completeness,
                true_muon_purity,
                true_muonPandoraPrimary,
                true_slcEff,
                true_slcPur
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees_protons.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(),
            std::vector<std::string>{
                // "event",
                "chi2_mu",
                "chi2_proton",
                "length",
                "length_ratio",
                "vtxDist",
                "trackScore",
                "completeness",
                "purity",
                "depEnergy",
                "pandoraPrimary",
                "sliceEfficiency", 
                "slicePurity"
            },
            loader,
            std::vector<ana::SpillMultiVar>{
                // event,
                reco_true_proton_chi2_mu,
                reco_true_proton_chi2_proton,
                reco_true_proton_length,
                reco_true_proton_length_ratio,
                reco_true_proton_vtxDist,
                reco_true_proton_trackScore,
                reco_true_proton_completeness,
                reco_true_proton_purity,
                reco_true_proton_depEnergy,
                reco_true_protonPandoraPrimary,
                reco_true_slcEff_proton,
                reco_true_slcPur_proton
            }, 
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees_protons.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(),
            std::vector<std::string>{
                // "event",
                "chi2_mu",
                "chi2_proton",
                "length",
                "length_ratio",
                "vtxDist",
                "trackScore",
                "completeness",
                "purity",
                "depEnergy",
                "pandoraPrimary",
                "sliceEfficiency", 
                "slicePurity"
            }, 
            loader,
            std::vector<ana::SpillMultiVar>{
                // event,
                true_proton_chi2_mu,
                true_proton_chi2_proton,
                true_proton_length,
                true_proton_length_ratio,
                true_proton_vtxDist,
                true_proton_trackScore,
                true_proton_completeness,
                true_proton_purity,
                true_proton_depEnergy,
                true_protonPandoraPrimary,
                true_slcEff_proton,
                true_slcPur_proton
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

}