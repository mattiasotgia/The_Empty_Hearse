#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Core/Var.h"

#include <vector>     
#include <tuple>
#include <set>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

const ana::SpillVar run    = ana::SIMPLESPILLVAR(hdr.run);
const ana::SpillVar subrun = ana::SIMPLESPILLVAR(hdr.subrun);
const ana::SpillVar event  = ana::SIMPLESPILLVAR(hdr.evt);

ana::Binning simple = ana::Binning::Simple(30, 0, 1);

std::vector<double> bad_events{
    10201, 10202, 10203, 10204, 10205, 10206, 10207, 10208, 10209, 10210, 10211, 10212, 10213, 10214, 10215, 10216, 10217, 10218, 10219, 10220, 10221, 10222, 10223, 10224, 10225, 10226, 10227, 10228, 10229, 10230, 10231, 10232, 10233, 10234, 10235, 10236, 10237, 10238, 10239, 10240, 10241, 10242, 10243, 10244, 10245, 10246, 10247, 10248, 10249, 10250, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 12551, 12552, 12553, 12554, 12555, 12556, 12557, 12558, 12559, 12560, 12561, 12562, 12563, 12564, 12565, 12566, 12567, 12568, 12569, 12570, 12571, 12572, 12573, 12574, 12575, 12576, 12577, 12578, 12579, 12580, 12581, 12582, 12583, 12584, 12585, 12586, 12587, 12588, 12589, 12590, 12591, 12592, 12593, 12594, 12595, 12596, 12597, 12598, 12599, 12600, 14951, 14952, 14953, 14954, 14955, 14956, 14957, 14958, 14959, 14960, 14961, 14962, 14963, 14964, 14965, 14966, 14967, 14968, 14969, 14970, 14971, 14972, 14973, 14974, 14975, 14976, 14977, 14978, 14979, 14980, 14981, 14982, 14983, 14984, 14985, 14986, 14987, 14988, 14989, 14990, 14991, 14992, 14993, 14994, 14995, 14996, 14997, 14998, 14999, 15000
};

ana::SpillVar get_event_non_cheated([](const caf::SRSpillProxy *sr) -> double {

    std::ofstream writer;
    writer.open("events_non_cheated_writer.txt", std::ios_base::app);
    
    writer << run(sr) << " " << subrun(sr) << " " << event(sr) << std::endl;
    return 0;
});

ana::SpillVar get_event_cheated([](const caf::SRSpillProxy *sr) -> double {

    std::ofstream writer;
    writer.open("events_cheated_writer.txt", std::ios_base::app);
    
    writer << run(sr) << " " << subrun(sr) << " " << event(sr) << std::endl;
    return 0;
});

void count_event () {

    ana::SpectrumLoader non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    ana::SpectrumLoader cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_stage1tocaf_flatcafs");

    ana::Spectrum non_cheated_spectrum("", simple, non_cheated, get_event_non_cheated, ana::kNoSpillCut);
    ana::Spectrum cheated_spectrum("", simple, cheated, get_event_cheated, ana::kNoSpillCut);

    non_cheated.Go();
    cheated.Go();

    return;
}
