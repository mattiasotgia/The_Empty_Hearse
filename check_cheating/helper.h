#ifndef HELPER_H
#define HELPER_H

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"

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

namespace logger {
    enum level {
        debug, 
        info, 
        warning, 
        error
    };

    std::map<level, std::string> log_level_string = {
        {debug, "debug"}, 
        {info, "info"}, 
        {warning, "warning"}, 
        {error, "error"}
    };

    std::ostream& log(level L, std::string where = "cout") {
        
        static std::ofstream file;

        std::ostream* writer = &std::cout;

        if (where != "cout") {
            file.open(where.c_str(), std::ios_base::app);
            writer = &file;
        }
        return *writer << "[" << log_level_string[L] << "] ";
    };
} // namespace log

namespace utils {
    int is_true_track(const caf::Proxy<caf::SRPFP> &pfp, const caf::SRSpillProxy *SR = nullptr) {
        if (pfp.trk.truth.p.pdg != pfp.shw.truth.p.pdg) {
            if (SR != nullptr)
                logger::log(logger::level::error) << "Found different pdg in run:event " << run(SR) << ":" << event(SR) << std::endl;
            else
                logger::log(logger::level::error) << "Found different pdg: " << pfp.trk.truth.p.pdg << " (pfp.trk.truth.p.pdg) and " << pfp.shw.truth.p.pdg << " (pfp.shw.truth.p.pdg)" << std::endl;
            return -1;
        }

        int true_pdg = pfp.trk.truth.p.pdg;

        std::vector<int> track_pdgs{
            13, 211, 2212
        };

        return std::find(track_pdgs.begin(), track_pdgs.end(), std::abs(true_pdg)) != track_pdgs.end();
    };  
    
    template<class T>
    const bool sanity_check(T object) {

        // Check normal types
        if (std::is_same_v<T, int> || std::is_same_v<T, double> || std::is_same_v<T, float>) {
            if (std::isnan(object)) return false;
        }

        // Check vectors
        if (std::is_same_v<T, caf::SRVector3D>) {
            if (std::isnan(object.x) || std::isnan(object.y) || std::isnan(object.z)) return false;
        }

        return true;
    }
} // namespace utils


namespace bins {
    const ana::Binning simple = ana::Binning::Simple(15, 0, 1);
    const ana::Binning simple_log = ana::Binning::LogUniform(50, 0, 1);
} // namespace bins

namespace hard_code {

    std::vector<unsigned> bad_events {
        10201, 10202, 10203, 10204, 10205, 10206, 10207, 10208, 10209, 10210, 10211, 10212, 10213, 10214, 10215, 10216, 10217, 10218, 10219, 10220, 10221, 10222, 10223, 10224, 10225, 10226, 10227, 10228, 10229, 10230, 10231, 10232, 10233, 10234, 10235, 10236, 10237, 10238, 10239, 10240, 10241, 10242, 10243, 10244, 10245, 10246, 10247, 10248, 10249, 10250, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 12551, 12552, 12553, 12554, 12555, 12556, 12557, 12558, 12559, 12560, 12561, 12562, 12563, 12564, 12565, 12566, 12567, 12568, 12569, 12570, 12571, 12572, 12573, 12574, 12575, 12576, 12577, 12578, 12579, 12580, 12581, 12582, 12583, 12584, 12585, 12586, 12587, 12588, 12589, 12590, 12591, 12592, 12593, 12594, 12595, 12596, 12597, 12598, 12599, 12600, 14951, 14952, 14953, 14954, 14955, 14956, 14957, 14958, 14959, 14960, 14961, 14962, 14963, 14964, 14965, 14966, 14967, 14968, 14969, 14970, 14971, 14972, 14973, 14974, 14975, 14976, 14977, 14978, 14979, 14980, 14981, 14982, 14983, 14984, 14985, 14986, 14987, 14988, 14989, 14990, 14991, 14992, 14993, 14994, 14995, 14996, 14997, 14998, 14999, 15000
    };

    const bool is_bad_event(const caf::SRSpillProxy *sr) {
        return std::find(bad_events.begin(), bad_events.end(), static_cast<int>(event(sr))) != bad_events.end();
    };
} // namespace hard_code

namespace cheating {
    using level_t = logger::level;
 
    // 
    // Definition of a chi2 statistics to test the cheat level (of the reconstruction goodness)
    template<class T>
    const T chi2(T __a, T __b, std::string label = "undef") {
        double chi2 = 0; // default "as if truthmatched";
        if (__b == 0 ) {
            if (__a == __b) {
                logger::log(level_t::warning) << "Found truth value zero (checking label " 
                                              << label << "), but reco == truth, so this will be ignored"
                                              << std::endl;
                return chi2;
            } else {
                logger::log(level_t::warning) << "Found truth value zero (checking label " 
                                              << label << "), but reco != truth, so returning their difference squared..."
                                              << std::endl;
                return std::pow(std::abs(__a - __b), 2);
            }
        }
        chi2 = std::pow(std::abs(__a - __b)/(double)__b, 2);
        logger::log(level_t::info) << label << " chi2: " << chi2 << "(__a = " << __a << "; __b = " << __b << ")" << std::endl;
        return chi2;
    };


    // The idea is to build a strong estimate (chi2-like) to test the goodness of 
    // this cheating... This is integrated over all the variables
    // At the moment the only possibility is ( V_reco - V_truth ) / V_truth
    // Minimal selection applied: iscc + pdg == 14 + event is not in the bad_list...
    const ana::MultiVar score ([](const caf::SRSliceProxy *slice) -> std::vector<double> {

        std::vector<double> score;

        score.push_back(chi2<double>(slice->vertex.x, slice->truth.position.x, "vertex_position"));
        score.push_back(chi2<double>(slice->vertex.y, slice->truth.position.y, "vertex_position"));
        score.push_back(chi2<double>(slice->vertex.z, slice->truth.position.z, "vertex_position"));

        // recob::PFParticle and other stuff
        for (auto const &pfp: slice->reco.pfp) {
            if (utils::is_true_track(pfp) == -1) continue;

            bool is_track = utils::is_true_track(pfp);

            // Track score as additional score value...
            // score.push_back(chi2<double>(pfp.trackScore, is_track ? 1 : 0, "track_score"));
            
            if (is_track) {
                // score.push_back(chi2<double>(pfp.trk.truth.bestmatch.hit_completeness, 1, "hit_completeness"));
                // score.push_back(chi2<double>(pfp.trk.truth.bestmatch.hit_purity, 1, "hit_purity"));

                // score.push_back(chi2<double>(pfp.trk.start.x, pfp.trk.truth.p.start.x, "track_start_point"));
                // score.push_back(chi2<double>(pfp.trk.start.y, pfp.trk.truth.p.start.y, "track_start_point"));
                // score.push_back(chi2<double>(pfp.trk.start.z, pfp.trk.truth.p.start.z, "track_start_point"));

                // score.push_back(chi2<double>(pfp.trk.end.x, pfp.trk.truth.p.end.x, "track_end_point"));
                // score.push_back(chi2<double>(pfp.trk.end.y, pfp.trk.truth.p.end.y, "track_end_point"));
                // score.push_back(chi2<double>(pfp.trk.end.z, pfp.trk.truth.p.end.z, "track_end_point")); 

                // score.push_back(chi2<double>(pfp.trk.len, pfp.trk.truth.p.length, "track_lenght"));

                // score.push_back(chi2<double>(
                //     slice->vertex.x-pfp.trk.start.x, 
                //     slice->truth.position.x-pfp.trk.truth.p.start.x, "vertex_start_position"
                // ));
                // score.push_back(chi2<double>(
                //     slice->vertex.y-pfp.trk.start.y, 
                //     slice->truth.position.y-pfp.trk.truth.p.start.y, "vertex_start_position"
                // ));
                // score.push_back(chi2<double>(
                //     slice->vertex.z-pfp.trk.start.z, 
                //     slice->truth.position.z-pfp.trk.truth.p.start.z, "vertex_start_position"
                // ));
            } else {
                score.push_back(chi2<double>(pfp.shw.truth.bestmatch.hit_completeness, 1));
                score.push_back(chi2<double>(pfp.shw.truth.bestmatch.hit_purity, 1));

                score.push_back(chi2<double>(pfp.shw.start.x, pfp.shw.truth.p.start.x, "shower_start_point"));
                score.push_back(chi2<double>(pfp.shw.start.y, pfp.shw.truth.p.start.y, "shower_start_point"));
                score.push_back(chi2<double>(pfp.shw.start.z, pfp.shw.truth.p.start.z, "shower_start_point"));

                score.push_back(chi2<double>(pfp.shw.end.x, pfp.shw.truth.p.end.x, "shower_end_point"));
                score.push_back(chi2<double>(pfp.shw.end.y, pfp.shw.truth.p.end.y, "shower_end_point"));
                score.push_back(chi2<double>(pfp.shw.end.z, pfp.shw.truth.p.end.z, "shower_end_point")); 
                
                score.push_back(chi2<double>(pfp.shw.len, pfp.shw.truth.p.length, "shower_lenght"));

                score.push_back(chi2<double>(
                    slice->vertex.x-pfp.shw.start.x,
                    slice->truth.position.x-pfp.shw.truth.p.start.x, "vertex_start_position"
                ));
                score.push_back(chi2<double>(
                    slice->vertex.y-pfp.shw.start.y,
                    slice->truth.position.y-pfp.shw.truth.p.start.y, "vertex_start_position"
                ));
                score.push_back(chi2<double>(
                    slice->vertex.z-pfp.shw.start.z,
                    slice->truth.position.z-pfp.shw.truth.p.start.z, "vertex_start_position"
                ));
            }

        } // loop pfp
        return score;
    }); // score ana::MultiVar

    std::map<bool, std::string> boolean_print = {
        {true, "true"},
        {false, "false"}
    };

    const ana::SpillMultiVar test_variables ([](const caf::SRSpillProxy *spill) -> std::vector<double> {
        std::vector<double> tested_variables;

        std::cout << "This is a spill, run:event = " << run(spill) << ":" << event(spill) << std::endl;

        std::cout << "Looking at truth interaction..." << std::endl;
        std::cout << "--> Found nnu = " << spill->mc.nnu << " true nu(s)" << std::endl;

        int inu = 0;
        for (auto const& nu: spill->mc.nu) {

            std::cout << "Nu no. " << inu << ": pdg = "  << nu.pdg << ", initpdg = " << nu.initpdg << ", index = " << nu.index << std::endl;
            
            std::cout << "          isnc = " << boolean_print[nu.isnc] << ", "
                      << "iscc = " << boolean_print[nu.iscc] << ", "
                      << "isvtxcont = " << boolean_print[nu.isvtxcont] << ", "
                      << "is_numucc_primary = " << boolean_print[nu.is_numucc_primary] << ", "
                      <<std::endl; 
            
            // std::cout << "          vtx = ("
            //           << nu.vtx.x << ", " 
            //           << nu.vtx.y << ", " 
            //           << nu.vtx.z << "), " << std::endl;

            std::cout << "          position = ("
                      << nu.position.x << ", " 
                      << nu.position.y << ", " 
                      << nu.position.z << "), " << std::endl;

            std::cout << "          npiplus = " << nu.npiplus << ", "
                      << "npiminus = " << nu.npiminus << ", "
                      << "npizero = " << nu.npizero << ", "
                      << "nproton = " << nu.nproton << ", "
                      << "nneutron = " << nu.nneutron << " (BEFORE final state interaction!)"
                      << std::endl;
            
            std::cout << "This nu has nprim = " << nu.nprim << " primary particle(s) " << std::endl;
            std::cout << "The primaries pf this interaction have pdg in {" << std::endl;

            for (auto const& prim: nu.prim) 
                std::cout << "\t" << prim.pdg << ", \tlenght = " << prim.length << " cm" << std::endl;

            std::cout << "}" << std::endl;

            inu ++ ;
        } // loop over nu(s)


        std::cout << "Looking at slices..." << std::endl;
        std::cout << "--> Found nslc = " << spill->nslc << " slice(s) " << std::endl;
        
        int islc = 0;
        for (auto const& slice: spill->slc) {

            std::cout << "Slice no. " << islc << " with ID = " << slice.self << std::endl;
            std::cout << "This is a " << ((slice.truth.iscc && slice.truth.pdg == 14) ? "GOOD" : "BAD") << " slice" << std::endl;

            std::cout << "Reco vertex in (" 
                      << slice.vertex.x << ", " 
                      << slice.vertex.y << ", " 
                      << slice.vertex.z << ") and true vertex in ("
                      << slice.truth.position.x << ", " 
                      << slice.truth.position.y << ", " 
                      << slice.truth.position.z << ")" 
                      << std::endl;
            std::cout << "This slice has npfp = " << slice.reco.npfp << " pfp(s)" << std::endl;
            std::cout << "Looping though primaries = {" << std::flush;

            for (auto const& iprimary: slice.primary) std::cout << iprimary << ", " << std::flush;

            std::cout << "} " << std::endl;

            int ipfp = 0;
            for (auto const& pfp: slice.reco.pfp) {

                
                std::cout << "[ipfp = " << ipfp << "] --> " 
                          << "ID = " << pfp.id << ", "
                          << "parent = " << pfp.parent << ", "
                          << "parent_is_primary = " << boolean_print[pfp.parent_is_primary] << ", "
                          << "trackScore = " << pfp.trackScore << ", "
                          << "slcID = " << pfp.slcID << std::endl;
                
                if (pfp.parent == -1) continue; // Neutrino PFP :) all other variables are empty


                ///////////////////////////////////////////
                ////////////////// TRACK //////////////////
                ///////////////////////////////////////////

                            //[ipfp = 0]
                std::cout << "           --> TRK properties: " << std::endl;
                std::cout << "                             .trk:         " 
                          << "len = " << pfp.trk.len << " cm, "
                          << "start = (" << pfp.trk.start.x << ", " << pfp.trk.start.y << ", " <<pfp.trk.start.z << "), "
                          << "end = (" << pfp.trk.end.x << ", " << pfp.trk.end.y << ", " <<pfp.trk.end.z << "), "
                          << std::endl
                          << "                                           "
                          << "dir = (" << pfp.trk.dir.x << ", " << pfp.trk.dir.y << ", " <<pfp.trk.dir.z << "), "
                          << "dir_end = (" << pfp.trk.dir_end.x << ", " << pfp.trk.dir_end.y << ", " <<pfp.trk.dir_end.z << "), "
                          << "costh = " << pfp.trk.costh << ", "
                          << "phi = " << pfp.trk.phi << ", "
                          << std::endl;

                std::cout << "                             .trk.truth.p: "
                          << "len = " << pfp.trk.truth.p.length << " cm, "
                          << "start = (" << pfp.trk.truth.p.start.x << ", " << pfp.trk.truth.p.start.y << ", " <<pfp.trk.truth.p.start.z << "), "
                          << "end = (" << pfp.trk.truth.p.end.x << ", " << pfp.trk.truth.p.end.y << ", " <<pfp.trk.truth.p.end.z << "), "
                          << "gen = (" << pfp.trk.truth.p.gen.x << ", " << pfp.trk.truth.p.gen.y << ", " <<pfp.trk.truth.p.gen.z << "), "
                          << std::endl
                          << "                                           "
                          << "pdg = " << pfp.trk.truth.p.pdg << ", "
                          << "cont_tpc = " << boolean_print[pfp.trk.truth.p.cont_tpc] << ", "
                          << "crosses_tpc = " << boolean_print[pfp.trk.truth.p.crosses_tpc] << ", "
                          << "contained = " << boolean_print[pfp.trk.truth.p.contained] << ", "
                          << "cryostat = " << pfp.trk.truth.p.cryostat
                          << std::endl;
                std::cout << "                             .trk.truth.bestmatch: "
                          << "hit_completeness = " << pfp.trk.truth.bestmatch.hit_completeness << ", "
                          << "hit_purity = " << pfp.trk.truth.bestmatch.hit_purity
                          << std::endl;

                ///////////////////////////////////////////
                ///////////////// SHOWER //////////////////
                ///////////////////////////////////////////
                
                std::cout << "           --> SHW properties: " << std::endl;
                std::cout << "                             .shw:         " 
                          << "len = " << pfp.shw.len << " cm, "
                          << "start = (" << pfp.shw.start.x << ", " << pfp.shw.start.y << ", " <<pfp.shw.start.z << "), "
                          << "end = (" << pfp.shw.end.x << ", " << pfp.shw.end.y << ", " <<pfp.shw.end.z << "), "
                          << std::endl
                          << "                                           "
                          << "dir = (" << pfp.shw.dir.x << ", " << pfp.shw.dir.y << ", " <<pfp.shw.dir.z << "), "
                          << "conversion_gap = " << pfp.shw.conversion_gap << " cm, "
                          << "density = " << pfp.shw.density << " MeV/cm, "
                          << "open_angle = " << pfp.shw.open_angle << " rad"
                          << std::endl;
                std::cout << "                             .shw.truth.p: "
                          << "len = " << pfp.shw.truth.p.length << " cm, "
                          << "start = (" << pfp.shw.truth.p.start.x << ", " << pfp.shw.truth.p.start.y << ", " <<pfp.shw.truth.p.start.z << "), "
                          << "end = (" << pfp.shw.truth.p.end.x << ", " << pfp.shw.truth.p.end.y << ", " <<pfp.shw.truth.p.end.z << "), "
                          << "gen = (" << pfp.shw.truth.p.gen.x << ", " << pfp.shw.truth.p.gen.y << ", " <<pfp.shw.truth.p.gen.z << "), "
                          << std::endl
                          << "                                           "
                          << "pdg = " << pfp.shw.truth.p.pdg << ", "
                          << "cont_tpc = " << boolean_print[pfp.shw.truth.p.cont_tpc] << ", "
                          << "crosses_tpc = " << boolean_print[pfp.shw.truth.p.crosses_tpc] << ", "
                          << "contained = " << boolean_print[pfp.shw.truth.p.contained] << ", "
                          << "cryostat = " << pfp.shw.truth.p.cryostat
                          << std::endl;
                std::cout << "                             .shw.truth.bestmatch: "
                          << "hit_completeness = " << pfp.shw.truth.bestmatch.hit_completeness << ", "
                          << "hit_purity = " << pfp.shw.truth.bestmatch.hit_purity
                          << std::endl;

                ipfp ++ ;
            } // loop over pfp(s)

            islc ++ ;
        } // loop over slice(s)

        return tested_variables;
    });

    const ana::Cut cut_numuCC ([](const caf::SRSliceProxy *slice) -> bool {
        logger::log(level_t::info) << "This is a " << ((slice->truth.iscc && slice->truth.pdg == 14) ? "GOOD" : "BAD") << " slice" << std::endl;
        return slice->truth.iscc && slice->truth.pdg == 14;
    }); 

    const ana::SpillCut cut_bad_events([](const caf::SRSpillProxy *spill) -> bool {
        // logger::log(level_t::info) << "This is a " << (hard_code::is_bad_event(spill) ? "BAD" : "GOOD") << " run:event" << std::endl;
        return !hard_code::is_bad_event(spill);
    });

} // namespace cheating

#endif // HELPER_H