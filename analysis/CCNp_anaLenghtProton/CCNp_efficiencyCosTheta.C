
#pragma once

#include "helper.h"
#include "selection.h"
#include "pdg_helper.h"

#include "sbnana/CAFAna/Core/Tree.h"
#include <map>

using level_t = logger::level;
using cut_type_t = var_utils::cut_type_t;

const ana::Cut makeCut_Np (const double& protonLength = 0.) {
    return ana::Cut ([=](const caf::SRSliceProxy *slice) -> bool {
        int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
        if (ipfp_muon == -1) return false; // redundant, btw, but who cares...

        int num_protons = 0;
        int num_pions = 0;
        int num_showers = 0;
        bool atLeastOneShort = false;

        for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp) {
            if (int(ipfp) == ipfp_muon)
                continue;
            if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::proton)  num_protons++;
            if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::pion)    num_pions++;
            if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::shower)  num_showers++;
            if (
                var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::proton &&
                slice->reco.pfp[ipfp].trk.len <= protonLength
            ) atLeastOneShort = true;
        } // loop pfp
        return num_protons > 0 && num_pions == 0 && num_showers == 0 && !atLeastOneShort;
    });
} 


const ana::Cut def_cut = (
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// true ANALYSIS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace true_vars {
    const ana::SpillVar L0_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0),          reco_true_1uNp,    def_cut_truth, 0);
    const ana::SpillVar L0_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0),          true_1uNp,         def_cut_truth, 0); 
    const ana::SpillVar L0_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0),          reco_1uNp,         def_cut_truth, 0);
    const ana::SpillVar L1_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0.31622777), reco_true_1uNp,    def_cut_truth, 0.31622777);
    const ana::SpillVar L1_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0.31622777), true_1uNp,         def_cut_truth, 0.31622777); 
    const ana::SpillVar L1_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0.31622777), reco_1uNp,         def_cut_truth, 0.31622777);
    const ana::SpillVar L2_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0.59948425), reco_true_1uNp,    def_cut_truth, 0.59948425);
    const ana::SpillVar L2_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0.59948425), true_1uNp,         def_cut_truth, 0.59948425); 
    const ana::SpillVar L2_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(0.59948425), reco_1uNp,         def_cut_truth, 0.59948425);
    const ana::SpillVar L3_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(1.13646367), reco_true_1uNp,    def_cut_truth, 1.13646367);
    const ana::SpillVar L3_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(1.13646367), true_1uNp,         def_cut_truth, 1.13646367); 
    const ana::SpillVar L3_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(1.13646367), reco_1uNp,         def_cut_truth, 1.13646367);
    const ana::SpillVar L4_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(2.15443469), reco_true_1uNp,    def_cut_truth, 2.15443469);
    const ana::SpillVar L4_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(2.15443469), true_1uNp,         def_cut_truth, 2.15443469); 
    const ana::SpillVar L4_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(2.15443469), reco_1uNp,         def_cut_truth, 2.15443469);
    const ana::SpillVar L5_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(4.08423865), reco_true_1uNp,    def_cut_truth, 4.08423865);
    const ana::SpillVar L5_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(4.08423865), true_1uNp,         def_cut_truth, 4.08423865); 
    const ana::SpillVar L5_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(4.08423865), reco_1uNp,         def_cut_truth, 4.08423865);
    const ana::SpillVar L6_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(7.74263683), reco_true_1uNp,    def_cut_truth, 7.74263683);
    const ana::SpillVar L6_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(7.74263683), true_1uNp,         def_cut_truth, 7.74263683); 
    const ana::SpillVar L6_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(7.74263683), reco_1uNp,         def_cut_truth, 7.74263683);
    const ana::SpillVar L7_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(14.67799268), reco_true_1uNp,   def_cut_truth, 14.67799268);
    const ana::SpillVar L7_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(14.67799268), true_1uNp,        def_cut_truth, 14.67799268); 
    const ana::SpillVar L7_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(14.67799268), reco_1uNp,        def_cut_truth, 14.67799268);
    const ana::SpillVar L8_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(27.82559402), reco_true_1uNp,   def_cut_truth, 27.82559402);
    const ana::SpillVar L8_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(27.82559402), true_1uNp,        def_cut_truth, 27.82559402);
    const ana::SpillVar L8_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(-9999, vars::pdg::slice_muonProtonTrueAngle, def_cut && makeCut_Np(27.82559402), reco_1uNp,        def_cut_truth, 27.82559402);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////std::vector<//////>//{}//////////////////////////////////////
////////////////////////////////////////// RECO ANALYSIS ////////////////////////////////////////////////////////////std::vector<//////>//{}//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////std::vector<//////>//{}//////////////////////////////////////

namespace reco_vars {
    const ana::SpillVar L0_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0),          reco_true_1uNp,    def_cut_truth, 0);
    const ana::SpillVar L0_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0),          true_1uNp,         def_cut_truth, 0); 
    const ana::SpillVar L0_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0),          reco_1uNp,         def_cut_truth, 0);
    const ana::SpillVar L1_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0.31622777), reco_true_1uNp,    def_cut_truth, 0.31622777);
    const ana::SpillVar L1_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0.31622777), true_1uNp,         def_cut_truth, 0.31622777); 
    const ana::SpillVar L1_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0.31622777), reco_1uNp,         def_cut_truth, 0.31622777);
    const ana::SpillVar L2_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0.59948425), reco_true_1uNp,    def_cut_truth, 0.59948425);
    const ana::SpillVar L2_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0.59948425), true_1uNp,         def_cut_truth, 0.59948425); 
    const ana::SpillVar L2_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(0.59948425), reco_1uNp,         def_cut_truth, 0.59948425);
    const ana::SpillVar L3_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(1.13646367), reco_true_1uNp,    def_cut_truth, 1.13646367);
    const ana::SpillVar L3_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(1.13646367), true_1uNp,         def_cut_truth, 1.13646367); 
    const ana::SpillVar L3_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(1.13646367), reco_1uNp,         def_cut_truth, 1.13646367);
    const ana::SpillVar L4_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(2.15443469), reco_true_1uNp,    def_cut_truth, 2.15443469);
    const ana::SpillVar L4_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(2.15443469), true_1uNp,         def_cut_truth, 2.15443469); 
    const ana::SpillVar L4_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(2.15443469), reco_1uNp,         def_cut_truth, 2.15443469);
    const ana::SpillVar L5_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(4.08423865), reco_true_1uNp,    def_cut_truth, 4.08423865);
    const ana::SpillVar L5_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(4.08423865), true_1uNp,         def_cut_truth, 4.08423865); 
    const ana::SpillVar L5_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(4.08423865), reco_1uNp,         def_cut_truth, 4.08423865);
    const ana::SpillVar L6_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(7.74263683), reco_true_1uNp,    def_cut_truth, 7.74263683);
    const ana::SpillVar L6_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(7.74263683), true_1uNp,         def_cut_truth, 7.74263683); 
    const ana::SpillVar L6_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(7.74263683), reco_1uNp,         def_cut_truth, 7.74263683);
    const ana::SpillVar L7_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(14.67799268), reco_true_1uNp,   def_cut_truth, 14.67799268);
    const ana::SpillVar L7_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(14.67799268), true_1uNp,        def_cut_truth, 14.67799268); 
    const ana::SpillVar L7_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(14.67799268), reco_1uNp,        def_cut_truth, 14.67799268);
    const ana::SpillVar L8_reco_true_cut  = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(27.82559402), reco_true_1uNp,   def_cut_truth, 27.82559402);
    const ana::SpillVar L8_true_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(27.82559402), true_1uNp,        def_cut_truth, 27.82559402);
    const ana::SpillVar L8_reco_cut       = var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>({}, vars::pdg::slice_muonProtonRecoAngle, def_cut && makeCut_Np(27.82559402), reco_1uNp,        def_cut_truth, 27.82559402);
}



void CCNp_efficiencyCosTheta() {
    std::unique_ptr<TFile> file_1muNp(new TFile("CCNp_anaMuonProtonAngle.root", "RECREATE"));

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

    std::vector<std::unique_ptr<ana::Tree>> trees_recoVar, trees_trueVar;

    for (auto const& running_loader: running_loaders) {
        
        ana::SpectrumLoader loader (loaders_available.at(running_loader).c_str());

        trees_recoVar.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(), 
            std::vector<std::string>{"L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"}, 
            loader,
            std::vector<ana::SpillVar>{
                reco_vars::L0_reco_true_cut, 
                reco_vars::L1_reco_true_cut, 
                reco_vars::L2_reco_true_cut, 
                reco_vars::L3_reco_true_cut, 
                reco_vars::L4_reco_true_cut, 
                reco_vars::L5_reco_true_cut, 
                reco_vars::L6_reco_true_cut, 
                reco_vars::L7_reco_true_cut, 
                reco_vars::L8_reco_true_cut
            },
            cuts::reco::spill_CRTPMTNeutrino,
	    true
        ));
        trees_recoVar.emplace_back(std::make_unique<ana::Tree>(
            ("reco_" + running_loader).c_str(), 
            std::vector<std::string>{"L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"}, 
            loader,
            std::vector<ana::SpillVar>{
                reco_vars::L0_reco_cut, 
                reco_vars::L1_reco_cut, 
                reco_vars::L2_reco_cut, 
                reco_vars::L3_reco_cut, 
                reco_vars::L4_reco_cut, 
                reco_vars::L5_reco_cut, 
                reco_vars::L6_reco_cut, 
                reco_vars::L7_reco_cut, 
                reco_vars::L8_reco_cut
            },
            cuts::reco::spill_CRTPMTNeutrino,
	    true
        ));
        trees_recoVar.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(), 
            std::vector<std::string>{"L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"}, 
            loader,
            std::vector<ana::SpillVar>{
                reco_vars::L0_true_cut, 
                reco_vars::L1_true_cut, 
                reco_vars::L2_true_cut, 
                reco_vars::L3_true_cut, 
                reco_vars::L4_true_cut, 
                reco_vars::L5_true_cut, 
                reco_vars::L6_true_cut, 
                reco_vars::L7_true_cut, 
                reco_vars::L8_true_cut
            },
            cuts::reco::spill_CRTPMTNeutrino,
	    true
        ));

        trees_trueVar.emplace_back(std::make_unique<ana::Tree>(
            ("reco_true_" + running_loader).c_str(), 
            std::vector<std::string>{"L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"}, 
            loader,
            std::vector<ana::SpillVar>{
                true_vars::L0_reco_true_cut, 
                true_vars::L1_reco_true_cut, 
                true_vars::L2_reco_true_cut, 
                true_vars::L3_reco_true_cut, 
                true_vars::L4_reco_true_cut, 
                true_vars::L5_reco_true_cut, 
                true_vars::L6_reco_true_cut, 
                true_vars::L7_reco_true_cut, 
                true_vars::L8_reco_true_cut
            },
            cuts::reco::spill_CRTPMTNeutrino,
	    true
        ));
        trees_trueVar.emplace_back(std::make_unique<ana::Tree>(
            ("reco_" + running_loader).c_str(), 
            std::vector<std::string>{"L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"}, 
            loader,
            std::vector<ana::SpillVar>{
                true_vars::L0_reco_cut, 
                true_vars::L1_reco_cut, 
                true_vars::L2_reco_cut, 
                true_vars::L3_reco_cut, 
                true_vars::L4_reco_cut, 
                true_vars::L5_reco_cut, 
                true_vars::L6_reco_cut, 
                true_vars::L7_reco_cut, 
                true_vars::L8_reco_cut
            },
            cuts::reco::spill_CRTPMTNeutrino,
	    true
        ));
        trees_trueVar.emplace_back(std::make_unique<ana::Tree>(
            ("true_" + running_loader).c_str(), 
            std::vector<std::string>{"L0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"}, 
            loader,
            std::vector<ana::SpillVar>{
                true_vars::L0_true_cut, 
                true_vars::L1_true_cut, 
                true_vars::L2_true_cut, 
                true_vars::L3_true_cut, 
                true_vars::L4_true_cut, 
                true_vars::L5_true_cut, 
                true_vars::L6_true_cut, 
                true_vars::L7_true_cut, 
                true_vars::L8_true_cut
            },
            cuts::reco::spill_CRTPMTNeutrino,
	    true
        ));

        loader.Go();
    }

    file_1muNp->mkdir("efficiency_studies_recoVar");
    file_1muNp->mkdir("efficiency_studies_trueVar");

    for (auto const& tree: trees_recoVar) 
        tree->SaveTo(file_1muNp->GetDirectory("efficiency_studies_recoVar"));

    for (auto const& tree: trees_trueVar) 
        tree->SaveTo(file_1muNp->GetDirectory("efficiency_studies_trueVar"));

};

