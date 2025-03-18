#include "helper.h"
#include "selection.h"

const ana::SpillVar vtx_dist ([](const caf::SRSpillProxy *spill) -> double {
    if (spill->slc.size() == 0) return -1; 
    for (auto const& slice: spill->slc) {
        if (cuts::truth::slice_vtx_in_FV(&slice)) {
            float dist = std::sqrt(
                std::pow(slice.vertex.x - slice.truth.position.x, 2) +
                std::pow(slice.vertex.y - slice.truth.position.y, 2) +
                std::pow(slice.vertex.z - slice.truth.position.z, 2)
            );
            if (dist > 50) std::cout << "event [" << spill->hdr.evt << "] vertex is distant diff = " << dist << std::endl;
            return dist;
        }
        return 0;
    }
    return 0;
});

void bad_vtx_event () {

    ana::SpectrumLoader loader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    ana::SpectrumLoader loader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_stage1tocaf_flatcafs");

    ana::Tree tests("tests", {"test_VTX_dist"}, loader_non_cheated, {vtx_dist}, ana::kNoSpillCut);
    loader_non_cheated.Go();

}

// /pnfs/icarus/scratch/users/msotgia/bnb_nuonly_cheating/v09_89_01_01p03/production_for_analysis/out/2d/59/msotgia_v09_89_01_01p03_BNB_production_v09_89_01_01p03_stage0_74960770_23-8561e788-48ac-4b6c-ac0f-c2fab0d4582b.root