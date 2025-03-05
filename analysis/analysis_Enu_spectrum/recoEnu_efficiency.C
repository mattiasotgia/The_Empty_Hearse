
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

// #include "selection_2025_numuCCana/mariaartero_selection.h"
#include "selection_2025_numuCCana/selection.h"
#include "plot.h"
#include "binning.h"

void recoEnu_efficiency() {

    // ana::SpectrumLoader loader150evt_noSAM("/pnfs/icarus/scratch/users/msotgia/bnb_nuonly_cheating/v09_89_01_01p03/production/out/**/**/msotgia*_stage0*root");

    ana::SpectrumLoader loader10evt_noSAM_fullcheat("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/50evt_tests/test_4_jobsub/stage1_full_cheat_essential_50.flat.caf.root");

    ana::SpectrumLoader loader10evt_noSAM_nocheat("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/50evt_tests/test_4_jobsub/stage1_no_cheat_essential_50.flat.caf.root");
    
    ana::Spectrum recoEnu_fullcheat("recoEnu", bins::nuE_15_normal, loader150evt_noSAM, selection::nu_energy_reco_1mu1p, ana::kNoSpillCut);
    loader10evt_noSAM_fullcheat.Go();

    ana::Spectrum recoEnu_nocheat("recoEnu", bins::nuE_15_normal, loader150evt_noSAM, selection::nu_energy_reco_1mu1p, ana::kNoSpillCut);
    loader10evt_noSAM_nocheat.Go();

};

#endif