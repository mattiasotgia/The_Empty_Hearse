
#pragma once

#include "helper.h"
#include "selection.h"
#include "slice_helper.h"

#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
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

const cut_type_t reco_true_1uNp = cut_type_t::BOTH_1muNp;
const cut_type_t true_1uNp = cut_type_t::TRUE_1muNp;
const cut_type_t reco_1uNp = cut_type_t::RECO;

// true + reco
const ana::SpillVar spill_true_E_reco_true_cut  = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_true_1uNp, def_cut_truth);

// true only
const ana::SpillVar spill_true_E_true_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, true_1uNp, def_cut_truth); 

// reco only
const ana::SpillVar spill_true_E_reco_cut       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy,       def_cut, reco_1uNp, def_cut_truth);

void CCNp_efficiencyHist() {

    std::map<std::string, std::string> loaders_available = {
        {"nominal",                  "/exp/icarus/data/users/cfarnese/ana_v096401/for_Mattia/nominal/*.root"},
        {"cheated_2d_vtx_3d_nu_mva", "/exp/icarus/data/users/cfarnese/ana_v096401/for_Mattia/full_cheating/*.root"}
    };

    // Running all :)
    std::vector<std::string> running_loaders = {
        "nominal",
        "cheated_2d_vtx_3d_nu_mva"
    };

    struct counts {
        double _reco_true, _reco, _true;
    };

    std::map<std::string, counts> EfficiencyPurityByLoader;

    for (auto const& running_loader: running_loaders) {
        
        ana::SpectrumLoader loader (loaders_available.at(running_loader).c_str());

        ana::Spectrum reco_true_spectrum("E_{#nu}", ana::Binning::Simple(1, 0.24, 2.4), loader, spill_true_E_reco_true_cut, cuts::reco::spill_CRTPMTNeutrino);
        ana::Spectrum true_spectrum("E_{#nu}", ana::Binning::Simple(1, 0.24, 2.4), loader, spill_true_E_true_cut, cuts::reco::spill_CRTPMTNeutrino);
        ana::Spectrum reco_spectrum("E_{#nu}", ana::Binning::Simple(1, 0.24, 2.4), loader, spill_true_E_reco_cut, cuts::reco::spill_CRTPMTNeutrino);

        loader.Go();

        EfficiencyPurityByLoader[running_loader] = {
            reco_true_spectrum.ToTH1(reco_true_spectrum.POT())->Integral(), 
            reco_spectrum.ToTH1(reco_spectrum.POT())->Integral(), 
            true_spectrum.ToTH1(true_spectrum.POT())->Integral()
        };
    }

    std::cout << "____ RESULTS ____ " << std::endl;

    for (auto const& running_loader: running_loaders) {
        std::cout << "Running on " << running_loader << ": reco_true = " << EfficiencyPurityByLoader[running_loader]._reco_true 
                  << ", reco = " << EfficiencyPurityByLoader[running_loader]._reco 
                  << ", true = " << EfficiencyPurityByLoader[running_loader]._true
                  << " (integrals from 0.24 to 2.4 GeV)" << std::endl;
        
        std::cout << "Running on " << running_loader << ": efficiency = " 
                  << EfficiencyPurityByLoader[running_loader]._reco_true/EfficiencyPurityByLoader[running_loader]._true
                  << ", purity = "
                  << EfficiencyPurityByLoader[running_loader]._reco_true/EfficiencyPurityByLoader[running_loader]._reco
                  << std::endl;
    }
};

