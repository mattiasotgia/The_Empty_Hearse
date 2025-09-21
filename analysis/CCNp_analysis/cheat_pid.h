#ifndef CLASSIFICATION_TYPES_CHEAT_PID
#define CLASSIFICATION_TYPES_CHEAT_PID

#include "selection.h"


namespace cheatPid {
    int cheatMuon (const caf::SRSliceProxy &slice, double dist_mucut) {

        /* Find the pfp number of the muon 
        */

        // Select muon as longest track
        double max_length = -1.0;
        int ipfp_mu = -1;
        TVector3 reco_vtx;
        TVector3 reco_start;
        
        reco_vtx.SetXYZ(slice.vertex.x, slice.vertex.y, slice.vertex.z);
        
        for (std::size_t ipfp = 0; ipfp < slice.reco.npfp; ++ipfp) {
            if (std::isnan(slice.reco.pfp[ipfp].trk.start.x) || std::isnan(slice.reco.pfp[ipfp].trk.len))
                continue;
        
            // if(!isInDetector(slice.reco.pfp[ipfp].trk.start.x,slice.reco.pfp[ipfp].trk.start.y,slice.reco.pfp[ipfp].trk.start.z))   continue;
            // if(!isInDetector(slice.reco.pfp[ipfp].trk.end.x,slice.reco.pfp[ipfp].trk.end.y,slice.reco.pfp[ipfp].trk.end.z))         continue;

            reco_start.SetXYZ(slice.reco.pfp[ipfp].trk.start.x, slice.reco.pfp[ipfp].trk.start.y, slice.reco.pfp[ipfp].trk.start.z);
            
            /* For µ the selection requires a trackScore > 0.5 */
            if (slice.reco.pfp[ipfp].trackScore < 0.5)
                continue;

            if (
                slice.reco.pfp[ipfp].trk.len > max_length   && 
                slice.reco.pfp[ipfp].trk.len > 50           && // The selection of the longest muon being 50+ cm should not be here? // separate Cut below
                (reco_vtx - reco_start).Mag() < dist_mucut  && 
                cuts::in_contained (slice.reco.pfp[ipfp].trk.end.x, slice.reco.pfp[ipfp].trk.end.y, slice.reco.pfp[ipfp].trk.end.z, 5.) && // separate Cut below
                slice.reco.pfp[ipfp].trk.end.x * slice.vertex.x > 0 && 
                slice.reco.pfp[ipfp].parent_is_primary && 
                std::abs(slice.reco.pfp[ipfp].trk.truth.p.pdg) == 13
            ) {
                max_length = slice.reco.pfp[ipfp].trk.len;
                ipfp_mu = ipfp;
            }
        } // loop of pfp to find muon
        return ipfp_mu;
    } // int find_muon

    particle_data::particle_t cheatPid (const caf::SRSliceProxy &slice, int ipfp, int dist_cut) {
        /* This utility select the correct particle of the event slice
         * return 1 PROTONS
         * return 2 PIONS
         * return 3 SHOWER
         * return 9 other -> nan, not primary, too far, below energy threshold...
        */

        TVector3 rec_vtx;
        rec_vtx.SetXYZ(slice.vertex.x, slice.vertex.y, slice.vertex.z);

        if (!(slice.reco.pfp[ipfp].parent_is_primary))
            return particle_data::particle_t::not_primary;

        // Sanity checks :)
        if (
            std::isnan(slice.reco.pfp[ipfp].trk.start.x)    || 
            std::isnan(slice.reco.pfp[ipfp].trk.end.x)      || 
            std::isnan(slice.reco.pfp[ipfp].trk.len)
        ) return particle_data::particle_t::undefined;
        
        // There is always a muon, for a 1mu1p we need 2 tracks - 1 muon = 1 only proton with threshold
        // if (int(ipfp)==ipfp_mu) continue;
        /* Consider only primary tracks which are 20cm 
         * close to the vertex, either vtx-start or vtx-end
        */

        TVector3 start (slice.reco.pfp[ipfp].trk.start.x, slice.reco.pfp[ipfp].trk.start.y, slice.reco.pfp[ipfp].trk.start.z);
        TVector3 end   (slice.reco.pfp[ipfp].trk.end.x,   slice.reco.pfp[ipfp].trk.end.y,   slice.reco.pfp[ipfp].trk.end.z);
        
        /* Some reminder of ternary operator
         * condition ? result_if_true : result_if_false
        */

        double min_dist = ((start - rec_vtx).Mag() < (end - rec_vtx).Mag()) ? (start - rec_vtx).Mag() : (end - rec_vtx).Mag();
        
        // if (min_dist > 50.0) continue;
        // if (min_dist > 10.0) return 9;
        if (min_dist > 50.0)
            return particle_data::particle_t::too_far;

        if (slice.reco.pfp[ipfp].trackScore >= 0.5) {
            if (
                std::isnan(slice.reco.pfp[ipfp].trk.start.x)    || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.x)      || 
                std::isnan(slice.reco.pfp[ipfp].trk.len)
            ) return particle_data::particle_t::undefined;

            if (
                std::isnan(slice.reco.pfp[ipfp].trk.start.y) || 
                std::isnan(slice.reco.pfp[ipfp].trk.start.z) || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.y)   || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.z)
            ) return particle_data::particle_t::undefined;

            // Skip low energy tagged pions
            TVector3 start_mom_V3;
        
            if (std::abs(slice.reco.pfp[ipfp].trk.truth.p.pdg) == 211) // pions
                start_mom_V3.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.z
                );
            
            if (
                std::abs(slice.reco.pfp[ipfp].trk.truth.p.pdg) == 211 && 
                (rec_vtx - start).Mag() < dist_cut  && 
                std::sqrt ( 
                    std::pow(particle_data::masses::pion, 2) + std::pow(start_mom_V3.Mag() * particle_data::GeV, 2) 
                ) - particle_data::masses::pion >= 25.0 && 
                slice.reco.pfp[ipfp].parent_is_primary &&
                (slice.reco.pfp[ipfp].trk.truth.p.genE * particle_data::GeV - particle_data::masses::pion) >= 25.0
            )  return particle_data::particle_t::pion;

            // Skip low energy protons
            if (std::abs(slice.reco.pfp[ipfp].trk.truth.p.pdg) == 2212)
                start_mom_V3.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.z
                );

            if (
                std::abs(slice.reco.pfp[ipfp].trk.truth.p.pdg) == 2212 && 
                (rec_vtx - start).Mag() < dist_cut  && 
                std::sqrt ( 
                    std::pow(particle_data::masses::proton, 2) + std::pow(start_mom_V3.Mag() * particle_data::GeV, 2) 
                ) - particle_data::masses::proton >= 50.0 && 
                slice.reco.pfp[ipfp].parent_is_primary && 
                (slice.reco.pfp[ipfp].trk.truth.p.genE * particle_data::GeV - particle_data::masses::proton) >= 50.0
            ) return particle_data::particle_t::proton;

        }

        if (slice.reco.pfp[ipfp].trackScore < 0.5) {
            if (slice.reco.pfp[ipfp].trackScore >= 0.4 && std::abs(slice.reco.pfp[ipfp].trk.truth.p.pdg) == 2212) {
                TVector3 start_mom_V3_2;
                start_mom_V3_2.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.z
                );

                if (
                    std::sqrt( 
                        std::pow(particle_data::masses::proton, 2) + std::pow(start_mom_V3_2.Mag() * particle_data::GeV, 2)
                    ) - particle_data::masses::proton >= 50.0 && 
                    (rec_vtx - start).Mag() < dist_cut && 
                    slice.reco.pfp[ipfp].parent_is_primary && 
                    (slice.reco.pfp[ipfp].trk.truth.p.genE * particle_data::GeV - particle_data::masses::proton) >= 50.0
                ) return particle_data::particle_t::proton;
            }

            if (!(slice.reco.pfp[ipfp].trackScore >= 0.4 && std::abs(slice.reco.pfp[ipfp].trk.truth.p.pdg) != 2212)) {

                // int use_plane2 = slice.reco.pfp[ipfp].trk.calo[2].nhit>slice.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
                int use_plane2 = 2;

                if (std::isnan(slice.reco.pfp[ipfp].shw.plane[use_plane2].energy))
                    return particle_data::particle_t::undefined;
                if (
                    slice.reco.pfp[ipfp].shw.plane[use_plane2].energy * particle_data::GeV < 25.0 && 
                    (slice.reco.pfp[ipfp].trk.truth.p.genE * particle_data::GeV) < 25.0
                )
                    return particle_data::particle_t::low_energy;
                if (
                    slice.reco.pfp[ipfp].shw.plane[use_plane2].energy * particle_data::GeV > 25.0 && 
                    slice.reco.pfp[ipfp].parent_is_primary && 
                    (slice.reco.pfp[ipfp].trk.truth.p.genE * particle_data::GeV) > 25.0
                )
                    return particle_data::particle_t::shower;
            }
        }
        return particle_data::particle_t::undefined;
    } // int id_pfp


    const ana::Cut cheatAllCut_1uNp ([](const caf::SRSliceProxy *slice) -> bool {
            int ipfp_muon = cheatMuon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return false; // redundant, btw, but who cares...

            int num_protons = 0;
            int num_pions = 0;
            int num_showers = 0;

            for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp) {
                if (int(ipfp) == ipfp_muon)
                    continue;

                if (cheatPid(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::proton)  num_protons++;
                if (cheatPid(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::pion)    num_pions++;
                if (cheatPid(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::shower)  num_showers++;
            } // loop pfp

            // std::cout << "FOUND RECO NUM proton #. " << num_protons << std::endl;
            return num_protons >= 1 && num_pions == 0 && num_showers == 0;
        });

        const ana::Cut cheatMuonCut ([](const caf::SRSliceProxy *slice) -> bool {
            /* We should at least found one muon, otherwise something bad happens
            */

            // Using dist_cut = 10 [cm]
            int ipfp_muon = cheatMuon(*slice, var_utils::dist_cut); 
            if (ipfp_muon == -1) return false; // no muon found ahaha
            return true;
        }); // const ana::Cut slice_at_least_mu

} // namespace cheatPid

#endif