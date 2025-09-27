
#pragma once

#include "helper.h"
#include "selection.h"
#include "slice_helper.h"
#include "cheat_pid.h"

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

const ana::Cut cheat_cut = (
    cheatPid::cheatMuonCut && cheatPid::cheatAllCut_1uNp && 
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
const ana::SpillVar reco_true_E      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_slcEff = SPILLVAR(-9999, vars::slice::slice_efficiency,      def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_slcPur = SPILLVAR(-9999, vars::slice::slice_purity,          def_cut, reco_true_1uNp, def_cut_truth);

const ana::SpillVar reco_true_E_cheatedPid      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, cheat_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_slcEff_cheatedPid = SPILLVAR(-9999, vars::slice::slice_efficiency,      cheat_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar reco_true_slcPur_cheatedPid = SPILLVAR(-9999, vars::slice::slice_purity,          cheat_cut, reco_true_1uNp, def_cut_truth);

// true only
const ana::SpillVar true_E      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_slcEff = SPILLVAR(-9999, vars::slice::slice_efficiency,      def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_slcPur = SPILLVAR(-9999, vars::slice::slice_purity,          def_cut, true_1uNp, def_cut_truth);

const ana::SpillVar true_E_cheatedPid      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, cheat_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_slcEff_cheatedPid = SPILLVAR(-9999, vars::slice::slice_efficiency,      cheat_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_slcPur_cheatedPid = SPILLVAR(-9999, vars::slice::slice_purity,          cheat_cut, true_1uNp, def_cut_truth); 

// reco only
const ana::SpillVar reco_E      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar reco_slcEff = SPILLVAR(-9999, vars::slice::slice_efficiency,      def_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar reco_slcPur = SPILLVAR(-9999, vars::slice::slice_purity,          def_cut, reco_1uNp, def_cut_truth);

const ana::SpillVar reco_E_cheatedPid      = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, cheat_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar reco_slcEff_cheatedPid = SPILLVAR(-9999, vars::slice::slice_efficiency,      cheat_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar reco_slcPur_cheatedPid = SPILLVAR(-9999, vars::slice::slice_purity,          cheat_cut, reco_1uNp, def_cut_truth);

const ana::SpillVar sliceRecoCount ([](const caf::SRSpillProxy *spill) -> double {
    
    int num_protons_above50 = 0;
    int num_muons = 0;
    int num_pions = 0;
    int num_neutral_pions = 0;
    int num_gamma = 0;
    double length_muon = 0;
    double dep_E = 0;

    int G4ID_parent;
    int use_plane = 2;

    int nslc = 0;

    for (auto const& slice: spill->slc) {
        particle_data::int_type_t interaction = var_utils::classification_type(spill, &slice);
        if (
            interaction.num_protons_above50 > 0 &&
            interaction.num_muons == 1          &&
            interaction.num_pions == 0          &&
            interaction.num_neutral_pions == 0  &&
            interaction.num_gamma == 0          &&
            interaction.length_muon > 50.       &&
            interaction.all_contained           &&
            !interaction.unclassified           &&
            def_cut_truth(&slice) & def_cut(&slice)
        ) {
            nslc++;
        }
    }
    return nslc;
});

const ana::SpillVar sliceTrueCount ([](const caf::SRSpillProxy *spill) -> double {
    
    int num_protons_above50 = 0;
    int num_muons = 0;
    int num_pions = 0;
    int num_neutral_pions = 0;
    int num_gamma = 0;
    double length_muon = 0;
    double dep_E = 0;

    int G4ID_parent;
    int use_plane = 2;

    int nslc = 0;

    for (auto const& slice: spill->slc) {
        particle_data::int_type_t interaction = var_utils::classification_type(spill, &slice);
        if (
            interaction.num_protons_above50 > 0 &&
            interaction.num_muons == 1          &&
            interaction.num_pions == 0          &&
            interaction.num_neutral_pions == 0  &&
            interaction.num_gamma == 0          &&
            interaction.length_muon > 50.       &&
            interaction.all_contained           &&
            !interaction.unclassified           &&
            def_cut_truth(&slice)
        ) {
            nslc++;
        }
    }
    return nslc;
});

const ana::SpillVar sliceMCCount ([](const caf::SRSpillProxy* spill) -> double {
    return spill->nslc;
});

void CCNp_efficiencySliceIssue () {

    // I fucked up all :)
    std::map<std::string, std::string> loaders_available = {
        {"nominal",         "msotgia_v09_89_01_01p03_sliceIssue_nominal"},
        {"cheatedSlice",    "msotgia_v09_89_01_01p03_sliceIssue_cheatedSlice"},
        {"twoPlaneSlice",   "msotgia_v09_89_01_01p03_sliceIssue_twoPlaneSlice"}
    };

    // Running all :)
    std::vector<std::string> running_loaders = {
        "nominal",
        "cheatedSlice",
        "twoPlaneSlice"
    };

    std::vector<std::unique_ptr<ana::Tree>> trees;

    for (auto const& running_loader: running_loaders) {
        
        ana::SpectrumLoader loader (loaders_available.at(running_loader).c_str());

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(), 
            std::vector<std::string>{"E", "sliceEfficiency", "slicePurity", "recoCount", "trueCount", "mcCount"}, 
            loader,
            std::vector<ana::SpillVar>{reco_true_E, reco_true_slcEff, reco_true_slcPur, sliceRecoCount, sliceTrueCount, sliceMCCount},
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader + "_cheatedPid").c_str(), 
            std::vector<std::string>{"E", "sliceEfficiency", "slicePurity", "recoCount", "trueCount", "mcCount"}, 
            loader,
            std::vector<ana::SpillVar>{reco_true_E_cheatedPid, reco_true_slcEff_cheatedPid, reco_true_slcPur_cheatedPid, sliceRecoCount, sliceTrueCount, sliceMCCount},
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_" + running_loader).c_str(), 
            std::vector<std::string>{"E", "sliceEfficiency", "slicePurity", "recoCount", "trueCount", "mcCount"}, 
            loader,
            std::vector<ana::SpillVar>{reco_E, reco_slcEff, reco_slcPur, sliceRecoCount, sliceTrueCount, sliceMCCount},
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("reco_" + running_loader + "_cheatedPid").c_str(), 
            std::vector<std::string>{"E", "sliceEfficiency", "slicePurity", "recoCount", "trueCount", "mcCount"}, 
            loader,
            std::vector<ana::SpillVar>{reco_E_cheatedPid, reco_slcEff_cheatedPid, reco_slcPur_cheatedPid, sliceRecoCount, sliceTrueCount, sliceMCCount},
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(), 
            std::vector<std::string>{"E", "sliceEfficiency", "slicePurity", "recoCount", "trueCount", "mcCount"}, 
            loader,
            std::vector<ana::SpillVar>{true_E, true_slcEff, true_slcPur, sliceRecoCount, sliceTrueCount, sliceMCCount},
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader + "_cheatedPid").c_str(), 
            std::vector<std::string>{"E", "sliceEfficiency", "slicePurity", "recoCount", "trueCount", "mcCount"}, 
            loader,
            std::vector<ana::SpillVar>{true_E_cheatedPid, true_slcEff_cheatedPid, true_slcPur_cheatedPid, sliceRecoCount, sliceTrueCount, sliceMCCount},
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1muNp(new TFile("CCNp_efficiencySliceIssue.root", "RECREATE"));
    file_1muNp->mkdir("sliceIssue");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1muNp->GetDirectory("sliceIssue"));
};

