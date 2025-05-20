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
    cuts::reco::slice_vtx_in_FV         &&
    cuts::reco::slice_barycenter        &&
    cuts::reco::slice_all_trk_contained
);

const ana::Cut def_cut_truth = (
    cuts::truth::slice_vtx_in_FV
);


const ana::SpillVar spill_true_E_numuCC             = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu1p, def_cut_truth);
const ana::SpillVar spill_true_E_numuCC_2p          = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu2p, def_cut_truth);
const ana::SpillVar spill_true_E_numuCC_3p          = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu3p, def_cut_truth);
const ana::SpillVar spill_true_E_numuCC_Ngreater3p  = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1muN3p, def_cut_truth);
const ana::SpillVar spill_true_E_numuNC             = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_numuNC);
const ana::SpillVar spill_true_E_nueCC              = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_nueCC);
const ana::SpillVar spill_true_E_numuCC_OOFV        = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu1p, !def_cut_truth);


const ana::SpillVar spill_true_E_QE                     = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_QE                 );
const ana::SpillVar spill_true_E_Res                    = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_Res                );
const ana::SpillVar spill_true_E_DIS                    = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_DIS                );
const ana::SpillVar spill_true_E_Coh                    = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_Coh                );
const ana::SpillVar spill_true_E_CohElastic             = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_CohElastic         );
const ana::SpillVar spill_true_E_ElectronScattering     = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_ElectronScattering );
const ana::SpillVar spill_true_E_IMDAnnihilation        = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_IMDAnnihilation    );
const ana::SpillVar spill_true_E_InverseBetaDecay       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_InverseBetaDecay   );
const ana::SpillVar spill_true_E_GlashowResonance       = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_GlashowResonance   );
const ana::SpillVar spill_true_E_AMNuGamma              = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_AMNuGamma          );
const ana::SpillVar spill_true_E_MEC                    = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_MEC                );
const ana::SpillVar spill_true_E_Diffractive            = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_Diffractive        );
const ana::SpillVar spill_true_E_EM                     = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_EM                 );
const ana::SpillVar spill_true_E_WeakMix                = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::RECO, def_cut_truth && cuts::truth::slice_WeakMix            );

const ana::SpillVar spill_true_E_dirt_1mu1pi        = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu1pi      , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_1mu2pi        = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu2pi      , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_1mu3pi        = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu3pi      , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_1muNpi        = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1muNpi      , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_2mu           = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_2mu         , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_2p            = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_2p          , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_1muShortNp    = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1muShortNp  , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_1pi1p         = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1pi1p       , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_1mu1p1pi0     = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu1p1pi0   , def_cut_truth);
const ana::SpillVar spill_true_E_dirt_1mu1pNpi0     = SPILLVAR(-9999, vars::truth::slice_neutrino_energy, def_cut, cut_type_t::BOTH_1mu1pNpi0   , def_cut_truth);


void reco_1mu1p_component_analysis() {

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
                "numuCC",
                "numuCC_2p",
                "numuCC_3p",
                "numuCC_Ngreater3p",
                "numuNC",
                "nueCC",
                "numuCC_OOFV",
                "genie_mode_QE",
                "genie_mode_Res",
                "genie_mode_DIS",
                "genie_mode_Coh",
                "genie_mode_CohElastic",
                "genie_mode_ElectronScattering",
                "genie_mode_IMDAnnihilation",
                "genie_mode_InverseBetaDecay",
                "genie_mode_GlashowResonance",
                "genie_mode_AMNuGamma",
                "genie_mode_MEC",
                "genie_mode_Diffractive",
                "genie_mode_EM",
                "genie_mode_WeakMix",
                "dirt_1mu1pi",
                "dirt_1mu2pi",
                "dirt_1mu3pi",
                "dirt_1muNpi",
                "dirt_2mu",
                "dirt_2p",
                "dirt_1pi1p",
                "dirt_1muShortNp",
                "dirt_1mu1p1pi0",
                "dirt_1mu1pNpi0"
            }, 
            loader,
            std::vector<ana::SpillVar>{
                event, 
                spill_true_E_numuCC,
                spill_true_E_numuCC_2p,
                spill_true_E_numuCC_3p,
                spill_true_E_numuCC_Ngreater3p,
                spill_true_E_numuNC,
                spill_true_E_nueCC,
                spill_true_E_numuCC_OOFV,
                spill_true_E_QE,
                spill_true_E_Res,
                spill_true_E_DIS,
                spill_true_E_Coh,
                spill_true_E_CohElastic,
                spill_true_E_ElectronScattering,
                spill_true_E_IMDAnnihilation,
                spill_true_E_InverseBetaDecay,
                spill_true_E_GlashowResonance,
                spill_true_E_AMNuGamma,
                spill_true_E_MEC,
                spill_true_E_Diffractive,
                spill_true_E_EM,
                spill_true_E_WeakMix,
                spill_true_E_dirt_1mu1pi,
                spill_true_E_dirt_1mu2pi,
                spill_true_E_dirt_1mu3pi,
                spill_true_E_dirt_1muNpi,
                spill_true_E_dirt_2mu,
                spill_true_E_dirt_2p,
                spill_true_E_dirt_1pi1p,
                spill_true_E_dirt_1muShortNp,
                spill_true_E_dirt_1mu1p1pi0,
                spill_true_E_dirt_1mu1pNpi0
            },
            cuts::reco::spill_CRTPMTNeutrino
        ));

        loader.Go();
    }

    std::unique_ptr<TFile> file_1mu1p(new TFile("pca_plot_1u1p.root", "RECREATE"));
    file_1mu1p->mkdir("component_analysis");
    for (auto const& tree: trees) 
        tree->SaveTo(file_1mu1p->GetDirectory("component_analysis"));
};
