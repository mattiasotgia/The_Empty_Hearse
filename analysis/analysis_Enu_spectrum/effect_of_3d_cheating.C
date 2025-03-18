#include "helper.h"
#include "selection.h"

const ana::SpillVar cheat_3d_completeness ([](const caf::SRSpillProxy *spill) -> double {
    if (spill->slc.size() == 0) return -1; 
    for (auto const& slice: spill->slc) {
        if (cuts::truth::slice_vtx_in_FV(&slice)) {

            // Select muon track 
            int ipfp_muon = var_utils::find_muon(slice, var_utils::dist_cut);

            if (ipfp_muon == -1) return -1;

            if (
                (slice.reco.pfp[ipfp_muon].trk.truth.bestmatch.hit_completeness < 0.7 && 
                slice.reco.pfp[ipfp_muon].trk.truth.bestmatch.hit_purity < 0.7 && 
                slice.reco.pfp[ipfp_muon].trk.len > 20) || spill->hdr.evt == 12271
            ) {
                std::cout << "event [" << spill->hdr.evt << "] hit completeness = " 
                          << slice.reco.pfp[ipfp_muon].trk.truth.bestmatch.hit_completeness << " and hit purity = " 
                          << slice.reco.pfp[ipfp_muon].trk.truth.bestmatch.hit_purity 
                          << " TRUTH is pdg = " << slice.reco.pfp[ipfp_muon].trk.truth.p.pdg << std::endl;
                          return -1;
            }
            
            return -1;
        }
    }
    return -1;
});

void effect_of_3d_cheating () {

    ana::SpectrumLoader loader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    ana::SpectrumLoader loader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_stage1tocaf_flatcafs");

    ana::Tree tests_no_cheat("tests", {"test_3d_cheating"}, loader_non_cheated, {cheat_3d_completeness}, ana::kNoSpillCut);
    loader_non_cheated.Go();

    ana::Tree tests_cheat("tests", {"test_3d_cheating"}, loader_cheated, {cheat_3d_completeness}, ana::kNoSpillCut);
    loader_cheated.Go();

}

// /pnfs/icarus/scratch/users/msotgia/bnb_nuonly_cheating/v09_89_01_01p03/production_for_analysis/out/2d/59/msotgia_v09_89_01_01p03_BNB_production_v09_89_01_01p03_stage0_74960770_23-8561e788-48ac-4b6c-ac0f-c2fab0d4582b.root



/*

Event 12271

*/