
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

const ana::SpillCut valid_events ([](const caf::SRSpillProxy *spill) -> bool {

    std::vector<unsigned> golden_events = {
        16537
        // 473911
    };
    return std::find(golden_events.begin(), golden_events.end(), static_cast<int>(event(spill))) != golden_events.end();
});

std::vector<std::string> loaders = {
    // "msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_2d_vtx_3d_nu",
    // "msotgia_v09_89_01_01p03_down_singles_both_ifdh_cheated_2d_vtx_3d",
    "msotgia_v09_89_01_01p03_down_singles_both_ifdh_new_nominal"
};

void CCNp_test () {

    for (const auto& l: loaders) {
        ana::SpectrumLoader loader(l.c_str());
        ana::Tree("test", {"test_var"}, loader, {cheating::test_variables}, valid_events); 
        loader.Go();
    }
    return;
}