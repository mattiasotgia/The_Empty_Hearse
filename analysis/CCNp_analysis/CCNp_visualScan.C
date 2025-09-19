
#pragma once

#include "helper.h"
#include "selection.h"
#include "slice_helper.h"
#include "pdg_helper.h"

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

const ana::Var selected ([](const caf::SRSliceProxy* slice) -> double {
    return def_cut(slice);
});

const cut_type_t reco_true_1uNp = cut_type_t::BOTH_1muNp;
const cut_type_t true_1uNp = cut_type_t::TRUE_1muNp;
const cut_type_t reco_1uNp = cut_type_t::RECO;

// true only
const ana::SpillVar true_muonLength = SPILLVAR(-9999, vars::pdg::visualScanning::slice_muon_length, def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_protonMinLength = SPILLVAR(-9999, vars::pdg::visualScanning::slice_shortestProtonLength, def_cut, true_1uNp, def_cut_truth);

const ana::SpillVar true_vertexX = SPILLVAR(-9999, vars::truth::slice_vertexX, def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_vertexY = SPILLVAR(-9999, vars::truth::slice_vertexY, def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_vertexZ = SPILLVAR(-9999, vars::truth::slice_vertexZ, def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar true_selected = SPILLVAR(-9999, selected, def_cut, true_1uNp, def_cut_truth);

void CCNp_visualScan () {

    // I fucked up all :)
    std::map<std::string, std::string> loaders_available = {
        {"nominal",                       "/storage/gpfs_data/icarus/plain/user/cfarnese/Mattia_nuonly_full_nominal/*/caf*/mc*.flat.caf.root"},
        {"cheated_2d",                    "/storage/gpfs_data/icarus/plain/user/cfarnese/Mattia_cheating_v0989_2D/*/caf*/mc*.flat.caf.root"},
        {"cheated_2d_vtx",                "/storage/gpfs_data/icarus/plain/user/cfarnese/Mattia_cheating_v0989_2D_Vtx/*/caf*/mc*.flat.caf.root"},
        {"cheated_2d_vtx_3d",             "/storage/gpfs_data/icarus/plain/user/cfarnese/Mattia_cheating_v0989_2D_Vtx_3D/*/caf*/mc*.flat.caf.root"},
        {"cheated_2d_vtx_3d_mva",         "/storage/gpfs_data/icarus/plain/user/cfarnese/Mattia_cheating_v0989_2D_Vtx_3D_Mva/*/caf*/mc*.flat.caf.root"},
        {"cheated_slicing",               "msotgia_v09_89_01_01p03_down_slicing_ifdh_nonu_respun_cheated_slice"},
        {"cheated_2d_slicing",            "msotgia_v09_89_01_01p03_down_slicing_ifdh_nonu_respun_cheated_2d_slice"},
        {"cheated_2d_vt_slicing",         "msotgia_v09_89_01_01p03_down_slicing_ifdh_nonu_respun_cheated_2d_vtx_slice"},
        {"cheated_2d_vtx_3d_slicing",     "msotgia_v09_89_01_01p03_down_slicing_ifdh_nonu_respun_cheated_2d_vtx_3d_slice"},
        {"cheated_2d_vtx_3d_mva_slicing", "msotgia_v09_89_01_01p03_down_slicing_ifdh_nonu_respun_cheated_2d_vtx_3d_mva_slice"}
    };

    // Running all :)
    std::vector<std::string> running_loaders = {
        "nominal",
        "cheated_2d",
        "cheated_2d_vtx",
        "cheated_2d_vtx_3d",
        "cheated_2d_vtx_3d_mva"
    };

    std::vector<std::unique_ptr<ana::Tree>> trees;

    for (auto const& running_loader: running_loaders) {
        
        ana::SpectrumLoader loader (loaders_available.at(running_loader).c_str());

        trees.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(), 
            std::vector<std::string>{
                "muonLength",
                "protonMinLength",
                "vertexX",
                "vertexY",
                "vertexZ",
                "selected"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                true_muonLength,
                true_protonMinLength,
                true_vertexX,
                true_vertexY,
                true_vertexZ,
                true_selected
            },
            cuts::reco::spill_CRTPMTNeutrino,
            true
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1muNp(new TFile("CCNp_visualScan.root", "RECREATE"));
    file_1muNp->mkdir("scanning");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1muNp->GetDirectory("scanning"));
};

