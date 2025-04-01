
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

#include "helper.h"
#include "selection.h"

#include "sbnana/CAFAna/Core/Tree.h"

#include "martero_util.h"

using level_t = logger::level;
using cut_type_t = var_utils::cut_type_t;

#define SPILLVAR(_def, _var, _reco, _what, _true) \
    var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(_def, _var, _reco, _what, _true)
#define SPILLMULTIVAR(_def, _var, _reco, _what, _true) \
    var_utils::make_spill_from_slice<ana::SpillMultiVar, ana::MultiVar, std::vector<double>>(_def, _var, _reco, _what, _true)

const ana::Cut def_cut = (
    // ana::kNoCut
    (cuts::reco::slice_1muNp || cuts::reco::slice_1mu1p)  &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    ana::kNoCut
    // cuts::truth::slice_vtx_in_FV
);

const cut_type_t local_cut_type = cut_type_t::BOTH_1muNp;

const ana::SpillVar count_events_mattia (
    const ana::Var &slice_var, 
    const ana::Cut &reco_cut = ana::kNoCut,  
    cut_type_t what_to_cut_on = cut_type_t::RECO, 
    const ana::Cut &truth_cut = ana::kNoCut,
    bool info = false, double treshold = 0.25
) {
    return ana::SpillVar([=](const caf::SRSpillProxy *spill) -> double {
        int selected_slices = 0;
        double slice_value = -9999;
        bool debug = false;

        for (auto const& slice: spill->slc) {

            if (var_utils::classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu1p) {
                selected_slices ++ ;
            }
        } // loop spill->slc
        return selected_slices;
    });
}

const ana::SpillVar count_slices = count_events_mattia (vars::truth::slice_neutrino_dE, def_cut, local_cut_type, def_cut_truth, true);
const ana::SpillVar spill_dE = SPILLVAR (-9999, vars::truth::slice_neutrino_dE, def_cut, local_cut_type, def_cut_truth);

void reco_1mu1p_checks() {


    ana::SpectrumLoader loader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    ana::Tree test("test", {"event", "mattia", "maria"}, loader_non_cheated, {event, spill_dE, k1mu1p_RecoMinusSlicetruthOverSlicetruth}, cuts::reco::spill_CRTPMTNeutrino);
    loader_non_cheated.Go();

    std::unique_ptr<TFile> file_reco_1muNp(new TFile("tests_mattia.root", "RECREATE"));
    test.SaveTo(file_reco_1muNp->mkdir("ana"));
};

#endif 
