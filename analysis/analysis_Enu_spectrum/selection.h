#ifndef SELECTION_H
#define SELECTION_H

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "helper.h"

#include <iostream>
#include <string>
#include <vector>

namespace cuts {

    bool in_FV (double x, double y, double z) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z))
            return false;
        // need to add a check to avoid having the vtx slice in one cryo and the end track in the other cryo

        // fiducial volume for the dangling cable
        if (x > 210.0 && y > 60.0 && z > 290.0 && z < 390.0)
            return false;
        // Cathode padding
        // if(std::abs(x)>208.79 && std::abs(x)<211.79) return false;

        return (((x < -61.94 - 25 && x > -358.49 + 25) || // cheack on x direction cryo E/0
                 (x > 61.94 + 25 && x < 358.49 - 25))  && // cheack on x direction cryo W/1
                ((y > -181.86 + 25 && y < 134.96 - 25) && // cheack on y direction
                 (z > -894.95 + 30 && z < 894.95 - 50))); // cheack on z direction
    } // bool in_FV

    bool in_contained (double x, double y, double z, double dist = 5) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z))
            return false;

        // 5 cm containment for Y in both cryo
        return (((x < -61.94 - dist && x > -358.49 + dist) || // cheack on x direction cryo E/0
                 (x > 61.94 + dist && x < 358.49 - dist))  && // cheack on x direction cryo W/1
                ((y > -181.86 + dist && y < 134.96 - dist) && // cheack on y direction
                (z > -894.95 + dist && z < 894.95 - dist)));  // cheack on z direction
    } // bool is_contained  

    bool in_active (double x, double y, double z) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z))
            return false;

        return (((x < -61.94 && x > -358.49) || // cheack on x direction cryo E/0
                 (x > 61.94 && x < 358.49))  && // cheack on x direction cryo W/1
                ((y > -181.86 && y < 134.96) && // cheack on y direction
                 (z > -894.95 && z < 894.95))); // cheack on z direction
    } // bool in_active


    bool in_detector (double x, double y, double z) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z))
            return false;
        return (((x < -61.94 + 5 && x > -358.49 - 5) || // cheack on x direction cryo E/0
                 (x > 61.94 - 5 && x < 358.49 + 5))  && // cheack on x direction cryo W/1
                ((y > -181.86 - 5 && y < 134.96 + 5) && // cheack on y direction
                 (z > -894.95 - 5 && z < 894.95 + 5))); // cheack on z direction
    } // bool in_detector

    namespace truth {
        const ana::Cut slice_numuCC = IsNumuCC;

        const ana::SpillCut spill_has_numuCC ([](const caf::SRSpillProxy *spill) -> bool {
            for (auto const& nu: spill->mc.nu) {
                if (nu.pdg == std::numeric_limits<int>::min()) continue;
                if (nu.iscc && nu.pdg == 14) return true;
            }
            return false;
        });

        const ana::SpillCut spill_all_trk_contained_MC ([](const caf::SRSpillProxy *spill) -> bool {
            // Cut the spill if the spill->mc.nu event (iscc && pdg == 14)) has anything 
            // or daughters going out of the contained volume
            
            for (auto const& nu: spill->mc.nu) {

                if (!nu.iscc || nu.pdg != 14) continue;

                for (auto const& prim: nu.prim) {
                    if (prim.G4ID < 0)     continue;
                    if (prim.cryostat < 0) continue;
                    if (
                        std::abs(prim.pdg) != 13   && 
                        std::abs(prim.pdg) != 2212 && 
                        std::abs(prim.pdg) != 211  && 
                        std::abs(prim.pdg) != 11
                    ) continue; // why are we selecting only charged primaries?
                        
                    if (!in_contained(prim.end.x, prim.end.y, prim.end.z)) return false;

                    if (prim.daughters.size() < 0) continue;
                    int parent_g4id;

                    for (auto const& true_particle: spill->true_particles) {
                        parent_g4id = true_particle.parent;

                        if (parent_g4id != prim.G4ID) continue;
                        if (
                            std::abs(true_particle.pdg) != 13   && 
                            std::abs(true_particle.pdg) != 2212 && 
                            std::abs(true_particle.pdg) != 211  && 
                            std::abs(true_particle.pdg) != 11
                        ) continue; // why are we selecting only charged primaries?

                        /* L’idea e’ mettere soglia su quello che deposita energia…. 
                         * Quindi tutto quello che sia carico e che di solito abbiamo
                         * Per dire in qualche punto ci deve essere anche una soglia 
                         * sui fotoni singoli, per cui se hanno troppa energia depositata 
                         * possono essere identificati (non avevamo messo soglia sui pi0 
                         * pero si sui due fotoni del pi0)
                        */

                        if (!in_contained(true_particle.end.x, true_particle.end.y, true_particle.end.z)) return false;
                    } // loop spill->true_particles
                } // loop nu.prim
            } // loop spill->mc.nu
        }); // const ana::SpillCut spill_all_trk_contained_MC

        const ana::SpillCut spill_all_trk_contained_truth ([](const caf::SRSpillProxy *spill) -> bool {
            // Cut the spill if the spill->slc slice (truth.iscc && truth.pdg == 14)) has anything 
            // or daughters going out of the contained volume

            for (auto const& slc: spill->slc) {
                if (!slc.truth.iscc || slc.truth.pdg != 14) continue;

                for (auto const& prim: slc.truth.prim) {
                    if (prim.G4ID < 0)     continue;
                    if (prim.cryostat < 0) continue;
                    if (
                        std::abs(prim.pdg) != 13   &&   // mu
                        std::abs(prim.pdg) != 2212 &&   // p
                        std::abs(prim.pdg) != 211  &&   // pi
                        std::abs(prim.pdg) != 11        // e 
                    ) continue; // why are we selecting only charged primaries?
                        
                    if (!in_contained(prim.end.x, prim.end.y, prim.end.z)) return false;

                    if (prim.daughters.size() < 0) continue;
                    int parent_g4id;

                    for (auto const& true_particle: spill->true_particles) {
                        parent_g4id = true_particle.parent;

                        if (parent_g4id != prim.G4ID) continue;
                        if (
                            std::abs(true_particle.pdg) != 13   &&  // mu
                            std::abs(true_particle.pdg) != 2212 &&  // p
                            std::abs(true_particle.pdg) != 211  &&  // pi
                            std::abs(true_particle.pdg) != 11       // e 
                        ) continue; // why are we selecting only charged primaries?
                        
                        /* L’idea e’ mettere soglia su quello che deposita energia…. 
                         * Quindi tutto quello che sia carico e che di solito abbiamo
                         * Per dire in qualche punto ci deve essere anche una soglia 
                         * sui fotoni singoli, per cui se hanno troppa energia depositata 
                         * possono essere identificati (non avevamo messo soglia sui pi0 
                         * pero si sui due fotoni del pi0)
                        */

                        if (!in_contained(true_particle.end.x, true_particle.end.y, true_particle.end.z)) return false;
                    } // loop spill->true_particles
                } // loop slc.truth.prim
            } // loop spill->slc
        }); // const ana::SpillCut spill_all_trk_contained_truth

    } // namespace truth

    namespace reco {
        const ana::Cut slice_all_trk_contained ([](const caf::SRSliceProxy *slice) -> bool {
            for (auto const& pfp: slice->reco.pfp) {
                if (std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x) || std::isnan(pfp.trk.len)) 
                    continue;

                // Check meaningful points
                // if (!isInDetector(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z)) continue;
                // if (!isInDetector(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z)) continue;
                // if (!(islc.reco.pfp[ipfp].parent_is_primary )) continue; // Skip secondaries
                // if (islc.reco.pfp[ipfp].trackScore<0.4) continue;        // Want to check only tracks??

                if ((pfp.trk.start.x * islc.vertex.x) < 0) return false; // PFP crossing cryostats :(
                if (!in_contained(pfp.trk.end.x, pfp.trk.end.y, pfp.trk.end.z, 5.)) 
                    return false;
            }
        }); // const ana::Cut slice_all_contained
    } // namespace reco
} //namespace cuts

namespace vars {
    namespace truth {

    } // namespace truth

    namespace reco {

    } // namespace reco
}

#endif