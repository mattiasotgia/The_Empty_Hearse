#include "helper.h"
// #include "selection.h"
#include "proton_helper.h"

#include "sbnana/CAFAna/Core/Tree.h"
#include <map>

#define SPILLVAR(_def, _var, _reco, _what, _true)                                                                               \
    var_utils::make_spill_from_slice<ana::SpillVar, ana::Var, double>(_def, _var, _reco, _what, _true)
#define SPILLMULTIVAR(_def, _var, _reco, _what, _true)                                                                          \
    var_utils::make_spill_from_slice<ana::SpillMultiVar, ana::MultiVar, std::vector<double>>(_def, _var, _reco, _what, _true)

using level_t = logger::level;
using cut_type_t = var_utils::cut_type_t;

const ana::Cut def_cut = (
    // ana::kNoCut
    // (cuts::reco::slice_1mu1p || cuts::reco::slice_1muNp) &&
    cuts::reco::slice_1mu1p             &&
    cuts::reco::slice_at_least_mu       &&
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_vtx_in_FV
);

const ana::SpillVar spill_reco_length_reco_1u1p                 = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_reco_length_numuCC                    = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_reco_length_numuCC_2p                 = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_reco_length_numuCC_3p                 = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_reco_length_numuCC_Ngreater3p         = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_reco_length_dirt_1mu1pi               = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_reco_length_dirt_1mu2pi               = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_reco_length_dirt_1mu3pi               = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_reco_length_dirt_1muNpi               = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_reco_length_dirt_2p                   = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_reco_length_dirt_1pi1p                = SPILLVAR(-9999, vars::reco::slice_pid_proton_reco_length, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);

const ana::SpillVar spill_true_length_reco_1u1p                 = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_true_length_numuCC                    = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_true_length_numuCC_2p                 = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_true_length_numuCC_3p                 = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_true_length_numuCC_Ngreater3p         = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_true_length_dirt_1mu1pi               = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_true_length_dirt_1mu2pi               = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_true_length_dirt_1mu3pi               = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_true_length_dirt_1muNpi               = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_true_length_dirt_2p                   = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_true_length_dirt_1pi1p                = SPILLVAR(-9999, vars::reco::slice_pid_proton_true_length, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);

const ana::SpillVar spill_L_reco_true_ratio_reco_1u1p           = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_numuCC              = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_numuCC_2p           = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_numuCC_3p           = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_numuCC_Ngreater3p   = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_dirt_1mu1pi         = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_dirt_1mu2pi         = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_dirt_1mu3pi         = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_dirt_1muNpi         = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_dirt_2p             = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_dirt_1pi1p          = SPILLVAR(-9999, vars::reco::slice_pid_proton_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);

const ana::SpillVar spill_hit_completeness_reco_1u1p            = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_hit_completeness_numuCC               = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_hit_completeness_numuCC_2p            = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_hit_completeness_numuCC_3p            = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_hit_completeness_numuCC_Ngreater3p    = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_hit_completeness_dirt_1mu1pi          = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_hit_completeness_dirt_1mu2pi          = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_hit_completeness_dirt_1mu3pi          = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_hit_completeness_dirt_1muNpi          = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_hit_completeness_dirt_2p              = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_hit_completeness_dirt_1pi1p           = SPILLVAR(-9999, vars::reco::slice_proton_hit_completeness, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);

const ana::SpillVar spill_hit_purity_reco_1u1p                  = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_hit_purity_numuCC                     = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_hit_purity_numuCC_2p                  = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_hit_purity_numuCC_3p                  = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_hit_purity_numuCC_Ngreater3p          = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_hit_purity_dirt_1mu1pi                = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_hit_purity_dirt_1mu2pi                = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_hit_purity_dirt_1mu3pi                = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_hit_purity_dirt_1muNpi                = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_hit_purity_dirt_2p                    = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_hit_purity_dirt_1pi1p                 = SPILLVAR(-9999, vars::reco::slice_proton_hit_purity, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);

const ana::SpillVar spill_R_reco_1u1p                           = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_R_numuCC                              = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_R_numuCC_2p                           = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_R_numuCC_3p                           = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_R_numuCC_Ngreater3p                   = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_R_dirt_1mu1pi                         = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_R_dirt_1mu2pi                         = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_R_dirt_1mu3pi                         = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_R_dirt_1muNpi                         = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_R_dirt_2p                             = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_R_dirt_1pi1p                          = SPILLVAR(-9999, vars::reco::slice_leading_proton_R, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);

const ana::SpillVar spill_track_score_reco_1u1p                 = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_track_score_numuCC                    = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_track_score_numuCC_2p                 = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_track_score_numuCC_3p                 = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_track_score_numuCC_Ngreater3p         = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_track_score_dirt_1mu1pi               = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_track_score_dirt_1mu2pi               = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_track_score_dirt_1mu3pi               = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_track_score_dirt_1muNpi               = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_track_score_dirt_2p                   = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_track_score_dirt_1pi1p                = SPILLVAR(-9999, vars::reco::slice_proton_track_score, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);

const ana::SpillVar spill_reco_Np_reco_1u1p                     = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_reco_Np_numuCC                        = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_reco_Np_numuCC_2p                     = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_reco_Np_numuCC_3p                     = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_reco_Np_numuCC_Ngreater3p             = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_reco_Np_dirt_1mu1pi                   = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_reco_Np_dirt_1mu2pi                   = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_reco_Np_dirt_1mu3pi                   = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_reco_Np_dirt_1muNpi                   = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_reco_Np_dirt_2p                       = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_reco_Np_dirt_1pi1p                    = SPILLVAR(-9999, vars::reco::slice_Np, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);


const ana::SpillVar spill_true_Np_reco_1u1p                     = vars::truth::spill_Np(-9999, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_true_Np_numuCC                        = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_true_Np_numuCC_2p                     = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1mu2p,      def_cut_truth);
const ana::SpillVar spill_true_Np_numuCC_3p                     = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1mu3p,      def_cut_truth);
const ana::SpillVar spill_true_Np_numuCC_Ngreater3p             = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1muN3p,     def_cut_truth);
const ana::SpillVar spill_true_Np_dirt_1mu1pi                   = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1mu1pi,     def_cut_truth);
const ana::SpillVar spill_true_Np_dirt_1mu2pi                   = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1mu2pi,     def_cut_truth);
const ana::SpillVar spill_true_Np_dirt_1mu3pi                   = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1mu3pi,     def_cut_truth);
const ana::SpillVar spill_true_Np_dirt_1muNpi                   = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1muNpi,     def_cut_truth);
const ana::SpillVar spill_true_Np_dirt_2p                       = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_true_Np_dirt_1pi1p                    = vars::truth::spill_Np(-9999, def_cut, cut_type_t::BOTH_1pi1p,      def_cut_truth);

void reco_1mu1p_protons() {

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
            std::vector<std::string>{
                "event", 
                "reco_length_reco_1u1p",
                "reco_length_numuCC",
                "reco_length_numuCC_2p",
                "reco_length_numuCC_3p",
                "reco_length_numuCC_Ngreater3p",
                "reco_length_dirt_1mu1pi",
                "reco_length_dirt_1mu2pi",
                "reco_length_dirt_1mu3pi",
                "reco_length_dirt_1muNpi",
                "reco_length_dirt_2p",
                "reco_length_dirt_1pi1p",
                "true_length_reco_1u1p",
                "true_length_numuCC",
                "true_length_numuCC_2p",
                "true_length_numuCC_3p",
                "true_length_numuCC_Ngreater3p",
                "true_length_dirt_1mu1pi",
                "true_length_dirt_1mu2pi",
                "true_length_dirt_1mu3pi",
                "true_length_dirt_1muNpi",
                "true_length_dirt_2p",
                "true_length_dirt_1pi1p",
                "L_reco_true_ratio_reco_1u1p",
                "L_reco_true_ratio_numuCC",
                "L_reco_true_ratio_numuCC_2p",
                "L_reco_true_ratio_numuCC_3p",
                "L_reco_true_ratio_numuCC_Ngreater3p",
                "L_reco_true_ratio_dirt_1mu1pi",
                "L_reco_true_ratio_dirt_1mu2pi",
                "L_reco_true_ratio_dirt_1mu3pi",
                "L_reco_true_ratio_dirt_1muNpi",
                "L_reco_true_ratio_dirt_2p",
                "L_reco_true_ratio_dirt_1pi1p",
                "hit_completeness_reco_1u1p",
                "hit_completeness_numuCC",
                "hit_completeness_numuCC_2p",
                "hit_completeness_numuCC_3p",
                "hit_completeness_numuCC_Ngreater3p",
                "hit_completeness_dirt_1mu1pi",
                "hit_completeness_dirt_1mu2pi",
                "hit_completeness_dirt_1mu3pi",
                "hit_completeness_dirt_1muNpi",
                "hit_completeness_dirt_2p",
                "hit_completeness_dirt_1pi1p",
                "hit_purity_reco_1u1p",
                "hit_purity_numuCC",
                "hit_purity_numuCC_2p",
                "hit_purity_numuCC_3p",
                "hit_purity_numuCC_Ngreater3p",
                "hit_purity_dirt_1mu1pi",
                "hit_purity_dirt_1mu2pi",
                "hit_purity_dirt_1mu3pi",
                "hit_purity_dirt_1muNpi",
                "hit_purity_dirt_2p",
                "hit_purity_dirt_1pi1p",
                "R_reco_1u1p",
                "R_numuCC",
                "R_numuCC_2p",
                "R_numuCC_3p",
                "R_numuCC_Ngreater3p",
                "R_dirt_1mu1pi",
                "R_dirt_1mu2pi",
                "R_dirt_1mu3pi",
                "R_dirt_1muNpi",
                "R_dirt_2p",
                "R_dirt_1pi1p",
                "track_score_reco_1u1p",
                "track_score_numuCC",
                "track_score_numuCC_2p",
                "track_score_numuCC_3p",
                "track_score_numuCC_Ngreater3p",
                "track_score_dirt_1mu1pi",
                "track_score_dirt_1mu2pi",
                "track_score_dirt_1mu3pi",
                "track_score_dirt_1muNpi",
                "track_score_dirt_2p",
                "track_score_dirt_1pi1p",
                "reco_Np_reco_1u1p",
                "reco_Np_numuCC",
                "reco_Np_numuCC_2p",
                "reco_Np_numuCC_3p",
                "reco_Np_numuCC_Ngreater3p",
                "reco_Np_dirt_1mu1pi",
                "reco_Np_dirt_1mu2pi",
                "reco_Np_dirt_1mu3pi",
                "reco_Np_dirt_1muNpi",
                "reco_Np_dirt_2p",
                "reco_Np_dirt_1pi1p",
                "true_Np_reco_1u1p",
                "true_Np_numuCC",
                "true_Np_numuCC_2p",
                "true_Np_numuCC_3p",
                "true_Np_numuCC_Ngreater3p",
                "true_Np_dirt_1mu1pi",
                "true_Np_dirt_1mu2pi",
                "true_Np_dirt_1mu3pi",
                "true_Np_dirt_1muNpi",
                "true_Np_dirt_2p",
                "true_Np_dirt_1pi1p"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                event, 
                spill_reco_length_reco_1u1p,
                spill_reco_length_numuCC,
                spill_reco_length_numuCC_2p,
                spill_reco_length_numuCC_3p,
                spill_reco_length_numuCC_Ngreater3p,
                spill_reco_length_dirt_1mu1pi,
                spill_reco_length_dirt_1mu2pi,
                spill_reco_length_dirt_1mu3pi,
                spill_reco_length_dirt_1muNpi,
                spill_reco_length_dirt_2p,
                spill_reco_length_dirt_1pi1p,
                spill_true_length_reco_1u1p,
                spill_true_length_numuCC,
                spill_true_length_numuCC_2p,
                spill_true_length_numuCC_3p,
                spill_true_length_numuCC_Ngreater3p,
                spill_true_length_dirt_1mu1pi,
                spill_true_length_dirt_1mu2pi,
                spill_true_length_dirt_1mu3pi,
                spill_true_length_dirt_1muNpi,
                spill_true_length_dirt_2p,
                spill_true_length_dirt_1pi1p,
                spill_L_reco_true_ratio_reco_1u1p,
                spill_L_reco_true_ratio_numuCC,
                spill_L_reco_true_ratio_numuCC_2p,
                spill_L_reco_true_ratio_numuCC_3p,
                spill_L_reco_true_ratio_numuCC_Ngreater3p,
                spill_L_reco_true_ratio_dirt_1mu1pi,
                spill_L_reco_true_ratio_dirt_1mu2pi,
                spill_L_reco_true_ratio_dirt_1mu3pi,
                spill_L_reco_true_ratio_dirt_1muNpi,
                spill_L_reco_true_ratio_dirt_2p,
                spill_L_reco_true_ratio_dirt_1pi1p,
                spill_hit_completeness_reco_1u1p,
                spill_hit_completeness_numuCC,
                spill_hit_completeness_numuCC_2p,
                spill_hit_completeness_numuCC_3p,
                spill_hit_completeness_numuCC_Ngreater3p,
                spill_hit_completeness_dirt_1mu1pi,
                spill_hit_completeness_dirt_1mu2pi,
                spill_hit_completeness_dirt_1mu3pi,
                spill_hit_completeness_dirt_1muNpi,
                spill_hit_completeness_dirt_2p,
                spill_hit_completeness_dirt_1pi1p,
                spill_hit_purity_reco_1u1p,
                spill_hit_purity_numuCC,
                spill_hit_purity_numuCC_2p,
                spill_hit_purity_numuCC_3p,
                spill_hit_purity_numuCC_Ngreater3p,
                spill_hit_purity_dirt_1mu1pi,
                spill_hit_purity_dirt_1mu2pi,
                spill_hit_purity_dirt_1mu3pi,
                spill_hit_purity_dirt_1muNpi,
                spill_hit_purity_dirt_2p,
                spill_hit_purity_dirt_1pi1p,
                spill_R_reco_1u1p,
                spill_R_numuCC,
                spill_R_numuCC_2p,
                spill_R_numuCC_3p,
                spill_R_numuCC_Ngreater3p,
                spill_R_dirt_1mu1pi,
                spill_R_dirt_1mu2pi,
                spill_R_dirt_1mu3pi,
                spill_R_dirt_1muNpi,
                spill_R_dirt_2p,
                spill_R_dirt_1pi1p,
                spill_track_score_reco_1u1p,
                spill_track_score_numuCC,
                spill_track_score_numuCC_2p,
                spill_track_score_numuCC_3p,
                spill_track_score_numuCC_Ngreater3p,
                spill_track_score_dirt_1mu1pi,
                spill_track_score_dirt_1mu2pi,
                spill_track_score_dirt_1mu3pi,
                spill_track_score_dirt_1muNpi,
                spill_track_score_dirt_2p,
                spill_track_score_dirt_1pi1p,
                spill_reco_Np_reco_1u1p,
                spill_reco_Np_numuCC,
                spill_reco_Np_numuCC_2p,
                spill_reco_Np_numuCC_3p,
                spill_reco_Np_numuCC_Ngreater3p,
                spill_reco_Np_dirt_1mu1pi,
                spill_reco_Np_dirt_1mu2pi,
                spill_reco_Np_dirt_1mu3pi,
                spill_reco_Np_dirt_1muNpi,
                spill_reco_Np_dirt_2p,
                spill_reco_Np_dirt_1pi1p,
                spill_true_Np_reco_1u1p,
                spill_true_Np_numuCC,
                spill_true_Np_numuCC_2p,
                spill_true_Np_numuCC_3p,
                spill_true_Np_numuCC_Ngreater3p,
                spill_true_Np_dirt_1mu1pi,
                spill_true_Np_dirt_1mu2pi,
                spill_true_Np_dirt_1mu3pi,
                spill_true_Np_dirt_1muNpi,
                spill_true_Np_dirt_2p,
                spill_true_Np_dirt_1pi1p
            },
            cuts::reco::spill_CRTPMTNeutrino
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1mu1p(new TFile("proton_purity_completeness_plot_1u1p.root", "RECREATE"));
    file_1mu1p->mkdir("proton");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1mu1p->GetDirectory("proton"));
};
