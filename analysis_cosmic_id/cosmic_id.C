#ifndef COSMIC_ID
#define COSMIC_ID

#include "selection_def.h"

void cosmic_id() {

    const auto binsNorm = ana::Binning::Simple(5, 0, 0.1);
    // ana::SpectrumLoader loaderNoCheated("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/event_1_processing/stage1_out_1muNp_nuCC.flat.caf.root");
    // ana::SpectrumLoader loaderNoCheated("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/event_1_processing_bare/stage1_out_1muNp_nuCC.flat.caf.root");
    ana::SpectrumLoader loaderNoCheated("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/event_1_processing_bare/caf_nocut_test/stage1_out_1muNp_nuCC.flat.caf.root");
    // ana::SpectrumLoader loaderNoCheated("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/event_1_processing_bare/caf_nocut_test_false/stage1_out_1muNp_nuCC.flat.caf.root");

    

    ana::Spectrum dataNoCheated("", binsNorm, loaderNoCheated, test, ana::kNoSpillCut);
    loaderNoCheated.Go();

    return;
}

#endif