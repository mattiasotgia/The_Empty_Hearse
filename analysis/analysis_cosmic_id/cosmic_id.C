#ifndef COSMIC_ID
#define COSMIC_ID

#include "selection_def.h"

void cosmic_id() {

    const auto binsNorm = ana::Binning::Simple(5, 0, 0.1);
    // ana::SpectrumLoader loaderNoCheated("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/event_1_processing/stage1_out_1muNp_nuCC.flat.caf.root");
    // ana::SpectrumLoader loaderNoCheated("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/event_1_processing_bare/stage1_out_1muNp_nuCC.flat.caf.root");
    // ana::SpectrumLoader loaderNoCheated("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/event_1_processing_bare/caf_nocut_test/stage1_out_1muNp_nuCC.flat.caf.root");
    // ana::SpectrumLoader loaderNoCheated("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/event_1_processing_bare/caf_nocut_test_false/stage1_out_1muNp_nuCC.flat.caf.root");

    ana::SpectrumLoader loader10evt_noSAM_fullcheat("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/50evt_tests/test_4_jobsub/stage1_full_cheat_essential_50.flat.caf.root");

    //    ana::SpectrumLoader loader10evt_noSAM_nocheat("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/50evt_tests/test_4_jobsub/stage1_no_cheat_essential_50.flat.caf.root");
    ana::SpectrumLoader loader10evt_noSAM_nocheat("/pnfs/icarus/scratch/users/msotgia/bnb_nuonly_cheating/v09_89_01_01p03/production_for_analysis/stage1tocaf/cheated/reco_ana/out/fb/59/msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_v09_89_01_01p03_stage1_cheated_18913322_0.flat.caf-745e4a70-db32-4719-8d1b-5ac4a394746f.root");

    ana::Spectrum dataNoCheated("", binsNorm, loader10evt_noSAM_nocheat, test, ana::kNoSpillCut);
    loader10evt_noSAM_nocheat.Go();

    ana::Spectrum dataFullCheated("", binsNorm, loader10evt_noSAM_fullcheat, test, ana::kNoSpillCut);
    loader10evt_noSAM_fullcheat.Go();
    return;
}

#endif
