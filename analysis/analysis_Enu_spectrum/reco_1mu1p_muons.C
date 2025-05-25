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

const ana::SpillVar spill_purity_reco_1u1p                  = SPILLVAR(-9999, vars::reco::slice_muon_hit_purity, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_purity_numuCC                     = SPILLVAR(-9999, vars::reco::slice_muon_hit_purity, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_purity_numuNC                     = SPILLVAR(-9999, vars::reco::slice_muon_hit_purity, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_numuNC);
const ana::SpillVar spill_purity_nueCC                      = SPILLVAR(-9999, vars::reco::slice_muon_hit_purity, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_nueCC);
const ana::SpillVar spill_purity_dirt_2mu                   = SPILLVAR(-9999, vars::reco::slice_muon_hit_purity, def_cut, cut_type_t::BOTH_2mu,        def_cut_truth);
const ana::SpillVar spill_purity_dirt_2p                    = SPILLVAR(-9999, vars::reco::slice_muon_hit_purity, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_purity_dirt_1muShortNp            = SPILLVAR(-9999, vars::reco::slice_muon_hit_purity, def_cut, cut_type_t::BOTH_1muShortNp, def_cut_truth);

const ana::SpillVar spill_completeness_reco_1u1p            = SPILLVAR(-9999, vars::reco::slice_muon_hit_completeness, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_completeness_numuCC               = SPILLVAR(-9999, vars::reco::slice_muon_hit_completeness, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_completeness_numuNC               = SPILLVAR(-9999, vars::reco::slice_muon_hit_completeness, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_numuNC);
const ana::SpillVar spill_completeness_nueCC                = SPILLVAR(-9999, vars::reco::slice_muon_hit_completeness, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_nueCC);
const ana::SpillVar spill_completeness_dirt_2mu             = SPILLVAR(-9999, vars::reco::slice_muon_hit_completeness, def_cut, cut_type_t::BOTH_2mu,        def_cut_truth);
const ana::SpillVar spill_completeness_dirt_2p              = SPILLVAR(-9999, vars::reco::slice_muon_hit_completeness, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_completeness_dirt_1muShortNp      = SPILLVAR(-9999, vars::reco::slice_muon_hit_completeness, def_cut, cut_type_t::BOTH_1muShortNp, def_cut_truth);

const ana::SpillVar spill_reco_length_reco_1u1p             = SPILLVAR(-9999, vars::reco::slice_pid_muon_reco_length, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_reco_length_numuCC                = SPILLVAR(-9999, vars::reco::slice_pid_muon_reco_length, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_reco_length_numuNC                = SPILLVAR(-9999, vars::reco::slice_pid_muon_reco_length, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_numuNC);
const ana::SpillVar spill_reco_length_nueCC                 = SPILLVAR(-9999, vars::reco::slice_pid_muon_reco_length, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_nueCC);
const ana::SpillVar spill_reco_length_dirt_2mu              = SPILLVAR(-9999, vars::reco::slice_pid_muon_reco_length, def_cut, cut_type_t::BOTH_2mu,        def_cut_truth);
const ana::SpillVar spill_reco_length_dirt_2p               = SPILLVAR(-9999, vars::reco::slice_pid_muon_reco_length, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_reco_length_dirt_1muShortNp       = SPILLVAR(-9999, vars::reco::slice_pid_muon_reco_length, def_cut, cut_type_t::BOTH_1muShortNp, def_cut_truth);

const ana::SpillVar spill_true_length_reco_1u1p             = SPILLVAR(-9999, vars::reco::slice_pid_muon_true_length, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_true_length_numuCC                = SPILLVAR(-9999, vars::reco::slice_pid_muon_true_length, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_true_length_numuNC                = SPILLVAR(-9999, vars::reco::slice_pid_muon_true_length, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_numuNC);
const ana::SpillVar spill_true_length_nueCC                 = SPILLVAR(-9999, vars::reco::slice_pid_muon_true_length, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_nueCC);
const ana::SpillVar spill_true_length_dirt_2mu              = SPILLVAR(-9999, vars::reco::slice_pid_muon_true_length, def_cut, cut_type_t::BOTH_2mu,        def_cut_truth);
const ana::SpillVar spill_true_length_dirt_2p               = SPILLVAR(-9999, vars::reco::slice_pid_muon_true_length, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_true_length_dirt_1muShortNp       = SPILLVAR(-9999, vars::reco::slice_pid_muon_true_length, def_cut, cut_type_t::BOTH_1muShortNp, def_cut_truth);

const ana::SpillVar spill_L_reco_true_ratio_reco_1u1p       = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_reco_true_ratio, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_numuCC          = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_numuNC          = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_reco_true_ratio, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_numuNC);
const ana::SpillVar spill_L_reco_true_ratio_nueCC           = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_reco_true_ratio, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_nueCC);
const ana::SpillVar spill_L_reco_true_ratio_dirt_2mu        = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_reco_true_ratio, def_cut, cut_type_t::BOTH_2mu,        def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_dirt_2p         = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_reco_true_ratio, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_L_reco_true_ratio_dirt_1muShortNp = SPILLVAR(-9999, vars::reco::slice_pid_muon_L_reco_true_ratio, def_cut, cut_type_t::BOTH_1muShortNp, def_cut_truth);

const ana::SpillVar spill_muon_track_score_reco_1u1p        = SPILLVAR(-9999, vars::reco::slice_muon_track_score, def_cut, cut_type_t::RECO,            def_cut_truth);
const ana::SpillVar spill_muon_track_score_numuCC           = SPILLVAR(-9999, vars::reco::slice_muon_track_score, def_cut, cut_type_t::BOTH_1mu1p,      def_cut_truth);
const ana::SpillVar spill_muon_track_score_numuNC           = SPILLVAR(-9999, vars::reco::slice_muon_track_score, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_numuNC);
const ana::SpillVar spill_muon_track_score_nueCC            = SPILLVAR(-9999, vars::reco::slice_muon_track_score, def_cut, cut_type_t::RECO,            def_cut_truth && cuts::truth::slice_nueCC);
const ana::SpillVar spill_muon_track_score_dirt_2mu         = SPILLVAR(-9999, vars::reco::slice_muon_track_score, def_cut, cut_type_t::BOTH_2mu,        def_cut_truth);
const ana::SpillVar spill_muon_track_score_dirt_2p          = SPILLVAR(-9999, vars::reco::slice_muon_track_score, def_cut, cut_type_t::BOTH_2p,         def_cut_truth);
const ana::SpillVar spill_muon_track_score_dirt_1muShortNp  = SPILLVAR(-9999, vars::reco::slice_muon_track_score, def_cut, cut_type_t::BOTH_1muShortNp, def_cut_truth);

void reco_1mu1p_muons() {

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
            (running_loader).c_str(), 
            std::vector<std::string>{
                "event", 
                "purity_reco_1u1p",
                "purity_numuCC",
                "purity_numuNC",
                "purity_nueCC",
                "purity_dirt_2mu",
                "purity_dirt_2p",
                "purity_dirt_1muShortNp",
                "completeness_reco_1u1p",
                "completeness_numuCC",
                "completeness_numuNC",
                "completeness_nueCC",
                "completeness_dirt_2mu",
                "completeness_dirt_2p",
                "completeness_dirt_1muShortNp",
                "reco_length_reco_1u1p",
                "reco_length_numuCC",
                "reco_length_numuNC",
                "reco_length_nueCC",
                "reco_length_dirt_2mu",
                "reco_length_dirt_2p",
                "reco_length_dirt_1muShortNp",
                "true_length_reco_1u1p",
                "true_length_numuCC",
                "true_length_numuNC",
                "true_length_nueCC",
                "true_length_dirt_2mu",
                "true_length_dirt_2p",
                "true_length_dirt_1muShortNp",
                "L_reco_true_ratio_reco_1u1p",
                "L_reco_true_ratio_numuCC",
                "L_reco_true_ratio_numuNC",
                "L_reco_true_ratio_nueCC",
                "L_reco_true_ratio_dirt_2mu",
                "L_reco_true_ratio_dirt_2p",
                "L_reco_true_ratio_dirt_1muShortNp",
                "spill_muon_track_score_reco_1u1p",
                "spill_muon_track_score_numuCC",
                "spill_muon_track_score_numuNC",
                "spill_muon_track_score_nueCC",
                "spill_muon_track_score_dirt_2mu",
                "spill_muon_track_score_dirt_2p",
                "spill_muon_track_score_dirt_1muShortNp"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                event, 
                spill_purity_reco_1u1p,
                spill_purity_numuCC,
                spill_purity_numuNC,
                spill_purity_nueCC,
                spill_purity_dirt_2mu,
                spill_purity_dirt_2p,
                spill_purity_dirt_1muShortNp,
                spill_completeness_reco_1u1p,
                spill_completeness_numuCC,
                spill_completeness_numuNC,
                spill_completeness_nueCC,
                spill_completeness_dirt_2mu,
                spill_completeness_dirt_2p,
                spill_completeness_dirt_1muShortNp,
                spill_reco_length_reco_1u1p,
                spill_reco_length_numuCC,
                spill_reco_length_numuNC,
                spill_reco_length_nueCC,
                spill_reco_length_dirt_2mu,
                spill_reco_length_dirt_2p,
                spill_reco_length_dirt_1muShortNp,
                spill_true_length_reco_1u1p,
                spill_true_length_numuCC,
                spill_true_length_numuNC,
                spill_true_length_nueCC,
                spill_true_length_dirt_2mu,
                spill_true_length_dirt_2p,
                spill_true_length_dirt_1muShortNp,
                spill_L_reco_true_ratio_reco_1u1p,
                spill_L_reco_true_ratio_numuCC,
                spill_L_reco_true_ratio_numuNC,
                spill_L_reco_true_ratio_nueCC,
                spill_L_reco_true_ratio_dirt_2mu,
                spill_L_reco_true_ratio_dirt_2p,
                spill_L_reco_true_ratio_dirt_1muShortNp,
                spill_muon_track_score_reco_1u1p,
                spill_muon_track_score_numuCC,
                spill_muon_track_score_numuNC,
                spill_muon_track_score_nueCC,
                spill_muon_track_score_dirt_2mu,
                spill_muon_track_score_dirt_2p,
                spill_muon_track_score_dirt_1muShortNp
            },
            cuts::reco::spill_CRTPMTNeutrino
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1mu1p(new TFile("muon_purity_completeness_plot_1u1p.root", "RECREATE"));
    file_1mu1p->mkdir("muon");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1mu1p->GetDirectory("muon"));
};
