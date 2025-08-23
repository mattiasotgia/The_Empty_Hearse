#pragma once

#include "helper.h"
#include "selection.h"

#include "sbnana/CAFAna/Core/Tree.h"
#include <map>

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

#define SPILLVAR(_def, _var, _reco, _what, _true)                                                                               \
    var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(_def, _var, _reco, _what, _true)


const ana::SpillVar CCNp_true ([](const caf::SRSpillProxy *spill) -> double {
    
    int num_protons_above50 = 0;
    int num_muons = 0;
    int num_pions = 0;
    int num_neutral_pions = 0;
    int num_gamma = 0;
    double length_muon = 0;
    double dep_E = 0;

    int G4ID_parent;
    int use_plane = 2;

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
            return interaction.num_protons_above50;
        }
    }
    return num_protons_above50;
});

const ana::SpillVar CCNp_true_reco_mu ([](const caf::SRSpillProxy *spill) -> double {
    
    int num_protons_above50 = 0;

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
            def_cut_truth(&slice) && cuts::reco::slice_at_least_mu(&slice)
        ) {
            return interaction.num_protons_above50;
        }
    }
    return num_protons_above50;
});

const cut_type_t reco_true_1uNp = cut_type_t::BOTH_1muNp;
const cut_type_t true_1uNp = cut_type_t::TRUE_1muNp;
const cut_type_t reco_1uNp = cut_type_t::RECO;

const ana::SpillVar CCNp_reco = SPILLVAR(0, vars::reco::slice_Np,  def_cut, reco_true_1uNp, def_cut_truth);

void CCNp_count_2D() {

    ana::SpectrumLoader nominal_loader                     ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_nominal");
    ana::SpectrumLoader cheated_2d_loader                  ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d");
    ana::SpectrumLoader cheated_vtx_loader                 ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_vtx");
    ana::SpectrumLoader cheated_vtxSelection_loader        ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_vtxSelection");
    ana::SpectrumLoader cheated_3d_loader                  ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_3d");
    ana::SpectrumLoader cheated_nuH_loader                 ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_nu");
    ana::SpectrumLoader cheated_mva_loader                 ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_mva");
    ana::SpectrumLoader cheated_2d_vtx_loader              ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d_vtx");
    ana::SpectrumLoader cheated_2d_vtx_3d_loader           ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d_vtx_3d");
    ana::SpectrumLoader cheated_2d_vtx_3d_nu_loader        ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d_vtx_3d_nu");
    ana::SpectrumLoader cheated_2d_vtx_3d_nu_mva_loader    ("msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_2d_vtx_3d_nu_mva");

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
        // "cheated_2d",
        // "cheated_vtx",
        // "cheated_vtxSelection",
        // "cheated_3d",
        // "cheated_nuH",
        // "cheated_mva",
        // "cheated_2d_vtx",
        // "cheated_2d_vtx_3d",
        // "cheated_2d_vtx_3d_nu",
        // "cheated_2d_vtx_3d_nu_mva"
    };

    std::vector<std::unique_ptr<ana::Tree>> trees;

    for (auto const& running_loader: running_loaders) {
        
        ana::SpectrumLoader& loader = *loaders_available.at(running_loader);

        trees.emplace_back(std::make_unique<ana::Tree>(
            (running_loader).c_str(), 
            std::vector<std::string>{"event", "reco_Np", "true_Np", "true_Np_reco_mu"}, 
            loader,
            std::vector<ana::SpillVar>{event, CCNp_reco, CCNp_true, CCNp_true_reco_mu},
            cuts::reco::spill_CRTPMTNeutrino
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1muNp(new TFile("CCNp_study.root", "RECREATE"));
    file_1muNp->mkdir("Np");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1muNp->GetDirectory("Np"));

}