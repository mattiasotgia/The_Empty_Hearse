
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

#include "helper.h"
#include "selection.h"

#include "sbnana/CAFAna/Core/Tree.h"

#define SPILLVAR(_def, _var, _reco, _what, _true)                                                                               \
    var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(_def, _var, _reco, _what, _true)
#define SPILLMULTIVAR(_def, _var, _reco, _what, _true)                                                                          \
    var_utils::make_spill_from_slice<ana::SpillMultiVar, ana::MultiVar, std::vector<double>>(_def, _var, _reco, _what, _true)

using level_t = logger::level;
using cut_type_t = var_utils::cut_type_t;

const ana::Cut def_cut = (
    // ana::kNoCut
    (cuts::reco::slice_1mu1p || cuts::reco::slice_1muNp) &&
    // cuts::reco::slice_1mu1p             &&
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
const particle_data::int_type_t type = particle_data::int_type_t::true_visible_1muNp;

// true + reco
const ana::SpillVar spill_reco_E_reco_true_cut  = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_true_cut  = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_true_1uNp, def_cut_truth);
const ana::SpillVar spill_reco_pT_reco_true_cut = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, reco_true_1uNp, def_cut_truth);

// true only
const ana::SpillVar spill_reco_E_true_cut       = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_true_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, true_1uNp, def_cut_truth);
const ana::SpillVar spill_reco_pT_true_cut      = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, true_1uNp, def_cut_truth);

// reco only
const ana::SpillVar spill_reco_E_reco_cut       = SPILLVAR(-9999, vars::reco::slice_neutrino_energy_1muNp,  def_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_true_E_reco_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_1uNp, def_cut_truth);
const ana::SpillVar spill_reco_pT_reco_cut      = SPILLVAR(-9999, vars::reco::slice_neutrino_pT_1muNp,      def_cut, reco_1uNp, def_cut_truth);

void reco_1muNp_eff() {

    ana::SpectrumLoader cheated_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_stage1tocaf_cheated_flatcaf");
    ana::SpectrumLoader nominal_loader("msotgia_v09_89_01_01p03_stage1_to_caf_reco_ana_stage1tocaf_nominal_flatcaf");

    std::unique_ptr<ana::Tree> cheated_reco_true_1uNp(new ana::Tree(
        "cheated_reco_true_1uNp", 
        {"event", "reco_E", "true_E", "reco_pT"},               ////< 3 labels 
        cheated_loader, 
        {event, spill_reco_E_reco_true_cut, spill_true_E_reco_true_cut, spill_reco_pT_reco_true_cut},     ////< 3 variables
        cuts::reco::spill_CRTPMTNeutrino
    )); 

    std::unique_ptr<ana::Tree> cheated_true_1uNp(new ana::Tree(
        "cheated_true_1uNp", 
        {"event", "reco_E", "true_E", "reco_pT"},               ////< 3 labels 
        cheated_loader, 
        {event, spill_reco_E_true_cut, spill_true_E_true_cut, spill_reco_pT_true_cut},     ////< 3 variables
        cuts::reco::spill_CRTPMTNeutrino
    )); 

    std::unique_ptr<ana::Tree> cheated_reco_1uNp(new ana::Tree(
        "cheated_reco_1uNp", 
        {"event", "reco_E", "true_E", "reco_pT"},               ////< 3 labels 
        cheated_loader, 
        {event, spill_reco_E_reco_cut, spill_true_E_reco_cut, spill_reco_pT_reco_cut},     ////< 3 variables
        cuts::reco::spill_CRTPMTNeutrino
    )); 
    
    std::unique_ptr<ana::Tree> nominal_reco_true_1uNp(new ana::Tree(
        "nominal_reco_true_1uNp", 
        {"event", "reco_E", "true_E", "reco_pT"},               ////< 3 labels 
        nominal_loader, 
        {event, spill_reco_E_reco_true_cut, spill_true_E_reco_true_cut, spill_reco_pT_reco_true_cut},     ////< 3 variables
	    cuts::reco::spill_CRTPMTNeutrino
    ));

    std::unique_ptr<ana::Tree> nominal_true_1uNp(new ana::Tree(
        "nominal_true_1uNp", 
        {"event", "reco_E", "true_E", "reco_pT"},               ////< 3 labels 
        nominal_loader, 
        {event, spill_reco_E_true_cut, spill_true_E_true_cut, spill_reco_pT_true_cut},     ////< 3 variables
	    cuts::reco::spill_CRTPMTNeutrino
    ));

    std::unique_ptr<ana::Tree> nominal_reco_1uNp(new ana::Tree(
        "nominal_reco_1uNp", 
        {"event", "reco_E", "true_E", "reco_pT"},               ////< 3 labels 
        nominal_loader, 
        {event, spill_reco_E_reco_cut, spill_true_E_reco_cut, spill_reco_pT_reco_cut},     ////< 3 variables
	    cuts::reco::spill_CRTPMTNeutrino
    ));

    cheated_loader.Go();
    nominal_loader.Go();

    std::unique_ptr<TFile> file_1mu1p(new TFile("efficiency_plot_1uNp.root", "RECREATE"));
    std::cout << "__ WRITING : cheated __ " << std::endl;
    cheated_reco_true_1uNp->SaveTo(file_1mu1p->mkdir("efficiency_plots"));
    cheated_true_1uNp->SaveTo(file_1mu1p->GetDirectory("efficiency_plots"));
    cheated_reco_1uNp->SaveTo(file_1mu1p->GetDirectory("efficiency_plots"));

    std::cout << "__ WRITING : nominal __ " << std::endl;
    nominal_reco_true_1uNp->SaveTo(file_1mu1p->GetDirectory("efficiency_plots"));
    nominal_true_1uNp->SaveTo(file_1mu1p->GetDirectory("efficiency_plots"));
    nominal_reco_1uNp->SaveTo(file_1mu1p->GetDirectory("efficiency_plots"));
};

#endif 
