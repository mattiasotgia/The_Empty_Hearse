#ifndef SELECTION_H
#define SELECTION_H

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TVector3.h"
#include "TProfile.h"

#include "helper.h"

#include <iostream>
#include <string>
#include <vector>

namespace particle_data {
    namespace masses {
        double pion = 139.57039;        // PDG value 2024 [MeV]
        double proton = 938.27208816;   // PDG value 2024 [MeV]
    }
} // namespace particle_data

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

    bool in_contained (double x, double y, double z, double dist = 5.) {
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

        const ana::Cut slice_vtx_in_FV ([](const caf::SRSliceProxy *slice) -> bool {
            return in_FV (slice->truth.position.x, slice->truth.position.y, slice->truth.position.z);
        }); // const ana::Cut slice_vtx_in_FV

        const ana::Cut slice_vtx_in_contained ([](const caf::SRSliceProxy *slice) -> bool {
            return in_contained (slice->truth.position.x, slice->truth.position.y, slice->truth.position.z, 5.);
        }); // const ana::Cut slice_vtx_in_contained

        const ana::Cut slice_vtx_in_active ([](const caf::SRSliceProxy *slice) -> bool {
            return in_active (slice->truth.position.x, slice->truth.position.y, slice->truth.position.z);
        }); // const ana::Cut slice_vtx_in_active

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
                // if (!isInDetector(slice.reco.pfp[ipfp].trk.start.x,slice.reco.pfp[ipfp].trk.start.y,slice.reco.pfp[ipfp].trk.start.z)) continue;
                // if (!isInDetector(slice.reco.pfp[ipfp].trk.end.x,slice.reco.pfp[ipfp].trk.end.y,slice.reco.pfp[ipfp].trk.end.z)) continue;
                // if (!(slice.reco.pfp[ipfp].parent_is_primary )) continue; // Skip secondaries
                // if (slice.reco.pfp[ipfp].trackScore<0.4) continue;        // Want to check only tracks??

                if ((pfp.trk.start.x * slice.vertex.x) < 0) return false; // PFP crossing cryostats :(
                if (!in_contained(pfp.trk.end.x, pfp.trk.end.y, pfp.trk.end.z, 5.)) 
                    return false;
            }
        }); // const ana::Cut slice_all_contained
        
        const ana::Cut slice_vtx_in_FV ([](const caf::SRSliceProxy *slice) -> bool {
            return in_FV (slice->vertex.x, slice->vertex.y, slice->vertex.z);
        }); // const ana::Cut slice_vtx_in_FV

        const ana::Cut slice_vtx_in_contained ([](const caf::SRSliceProxy *slice) -> bool {
            return in_contained (slice->vertex.x, slice->vertex.y, slice->vertex.z, 5.);
        }); // const ana::Cut slice_vtx_in_contained

        const ana::Cut slice_vtx_in_active ([](const caf::SRSliceProxy *slice) -> bool {
            return in_active (slice->vertex.x, slice->vertex.y, slice->vertex.z);
        }); // const ana::Cut slice_vtx_in_active
    } // namespace reco
} //namespace cuts

namespace var_utils {

    std::string dEdx_temp = 
        "/exp/icarus/app/users/msotgia/analysis/sbnana_v09_93_01_thesis_analysis/analysis/dEdxrestemplates.root"
    TFile* file = TFile::Open(dEdx_temp.c_str());

    auto dedx_range_pro = (TProfile *)file->Get("dedx_range_pro");
    auto dedx_range_ka  = (TProfile *)file->Get("dedx_range_ka");
    auto dedx_range_pi  = (TProfile *)file->Get("dedx_range_pi");
    auto dedx_range_mu  = (TProfile *)file->Get("dedx_range_mu");

    struct chi2
    {
        double muon;
        double proton;
        double kaon;
        double pi;
    }; 
    

    chi2 chi2_ALG (std::vector<double> &dEdx, std::vector<double> &RR, double rr_min, double rr_max) {
        /* The output is chi2s
         * Those are used to classify final state particles 
        */

        double threshold = 0.5;
        double max_rr = rr_max; // Max value for the residual range (RR)
        double min_rr = rr_min; // Min value for the residual range (RR)

        std::vector<float> trkdedx;
        std::vector<float> trkres;
        std::vector<double> vpida;

        for (std::size_t i(0); i < dEdx.size(); ++i) {
            if (i == 0 || i == dEdx.size() - 1)
                continue;
            if (RR[i] < max_rr && RR[i] > rr_min) {
                trkdedx.push_back(dEdx[i]);
                trkres.push_back(RR[i]);
            }
        }

        int npt = 0;
        double chi2pro = 0;
        double chi2ka = 0;
        double chi2pi = 0;
        double chi2mu = 0;
        double avgdedx = 0;
        double PIDA = 0;

        int used_trkres = 0;
        for (unsigned i = 0; i < trkdedx.size(); ++i) { 
            // hits
            // ignore the first and the last point
            // if (i == 0 || i == trkdedx.size() - 1) continue;
            avgdedx += trkdedx[i];
            if (trkres[i] < 26) {
                PIDA += trkdedx[i] * pow(trkres[i], 0.42);
                vpida.push_back(trkdedx[i] * pow(trkres[i], 0.42));
                used_trkres++;
            }
            if (trkdedx[i] > 100 || trkdedx[i] < threshold)
                continue; // protect against large pulse height
            
            int bin = dedx_range_pro->FindBin(trkres[i]);
            if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
                
                double bincpro = dedx_range_pro->GetBinContent(bin);
                if (bincpro < 1e-6)  // for 0 bin content, using neighboring bins
                    bincpro = 
                        (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
                
                double bincka = dedx_range_ka->GetBinContent(bin);
                if (bincka < 1e-6)
                    bincka =
                        (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;

                double bincpi = dedx_range_pi->GetBinContent(bin);
                if (bincpi < 1e-6)
                    bincpi =
                        (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
                
                double bincmu = dedx_range_mu->GetBinContent(bin);
                if (bincmu < 1e-6)
                    bincmu =
                        (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
                
                double binepro = dedx_range_pro->GetBinError(bin);
                if (binepro < 1e-6)
                    binepro =
                        (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;

                double bineka = dedx_range_ka->GetBinError(bin);
                if (bineka < 1e-6)
                    bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;

                double binepi = dedx_range_pi->GetBinError(bin);
                if (binepi < 1e-6)
                    binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;

                double binemu = dedx_range_mu->GetBinError(bin);
                if (binemu < 1e-6)
                    binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;

                // double errke = 0.05*trkdedx[i];   //5% KE resolution

                double errdedx = 0.04231 + 0.0001783 * trkdedx[i] * trkdedx[i]; // resolution on dE/dx
                
                errdedx *= trkdedx[i];
                
                chi2pro += pow((trkdedx[i] - bincpro) / std::sqrt(pow(binepro, 2) + pow(errdedx, 2)), 2);
                chi2ka += pow((trkdedx[i] - bincka) / std::sqrt(pow(bineka, 2) + pow(errdedx, 2)), 2);
                chi2pi += pow((trkdedx[i] - bincpi) / std::sqrt(pow(binepi, 2) + pow(errdedx, 2)), 2);
                chi2mu += pow((trkdedx[i] - bincmu) / std::sqrt(pow(binemu, 2) + pow(errdedx, 2)), 2);
                // std::cout<<i<<" "<<trkdedx[i]<<" "<<trkres[i]<<" "<<bincpro<<std::endl;
                ++npt;
            }
        } // hits

        return {chi2mu / npt, chi2pro / npt, chi2ka / npt, chi2pi / npt};
    } // chi2 chi2_ALG

    int find_muon (const caf::SRSliceProxy &slice, int dist_mucut) {

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

            // int use_plane = slice.reco.pfp[ipfp].trk.calo[2].nhit>slice.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
            int use_plane = 2; // Always using collection?
            
            // compute new chi2
            std::vector<double> dedx;
            std::vector<double> rr;
            for (std::size_t ihit = 0; ihit < slice.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit) {
                dedx.push_back(slice.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
                rr.push_back(slice.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
            } // calo points

            // input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
            // output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt}
            chi2 chi2_values = chi2_ALG(dedx, rr, 0.0, 25.0);

            if (
                slice.reco.pfp[ipfp].trk.len > max_length   && 
                slice.reco.pfp[ipfp].trk.len > 50           &&
                (reco_vtx - reco_start).Mag() < dist_mucut  && 
                chi2_values.muon < 30                       && // chi2 cuts on proton and muon 
                chi2_values.proton > 60                     && // chi2 cuts on proton and muon 
                in_contained (slice.reco.pfp[ipfp].trk.end.x, slice.reco.pfp[ipfp].trk.end.y, slice.reco.pfp[ipfp].trk.end.z, 5.) &&
                slice.reco.pfp[ipfp].trk.end.x * slice.vertex.x > 0 && 
                slice.reco.pfp[ipfp].parent_is_primary
            ) {
                max_length = slice.reco.pfp[ipfp].trk.len;
                ipfp_mu = ipfp;
            }
        } // loop of pfp to find muon
        return ipfp_mu;
    } // int find_muon

    int id_pfp (const caf::SRSliceProxy &slice, int ipfp, int dist_cut) {
        /* This utility select the correct particle of the event slice
         * return 1 PROTONS
         * return 2 PIONS
         * return 3 SHOWER
         * return 9 other -> nan, not primary, too far, below energy threshold...
        */

        TVector3 rec_vtx;
        rec_vtx.SetXYZ(slice.vertex.x, slice.vertex.y, slice.vertex.z);

        if (!(slice.reco.pfp[ipfp].parent_is_primary))
            return 9;

        // Sanity checks :)
        if (
            std::isnan(slice.reco.pfp[ipfp].trk.start.x)    || 
            std::isnan(slice.reco.pfp[ipfp].trk.end.x)      || 
            std::isnan(slice.reco.pfp[ipfp].trk.len)
        ) return 9;
        
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
            return 9;

        // int use_plane = slice.reco.pfp[ipfp].trk.calo[2].nhit>slice.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
        int use_plane = 2;
        // compute new chi2
        std::vector<double> dedx;
        std::vector<double> rr;
        for (std::size_t ihit = 0; ihit < slice.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit) {
            dedx.push_back(slice.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
            rr.push_back(slice.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
        } // calo points

        // input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
        // output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt}
        chi2 chi2_values = chi2_ALG(dedx, rr, 0.0, 25.0);

        if (slice.reco.pfp[ipfp].trackScore >= 0.5) {
            if (
                std::isnan(slice.reco.pfp[ipfp].trk.start.x)    || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.x)      || 
                std::isnan(slice.reco.pfp[ipfp].trk.len)
            ) return 9;

            if (
                std::isnan(slice.reco.pfp[ipfp].trk.start.y) || 
                std::isnan(slice.reco.pfp[ipfp].trk.start.z) || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.y)   || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.z)
            ) return 9;

            // Skip low energy tagged pions
            TVector3 start_mom_V3;
        
            if (chi2_values.proton >= 100) 
                start_mom_V3.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.z
                );
            
            if (
                chi2_values.proton >= 100           && 
                (rec_vtx - start).Mag() < dist_cut  && 
                sqrt ( pow(particle_data::masses::pion, 2) + pow(start_mom_V3.Mag() * 1000, 2) ) - particle_data::masses::pion >= 25.0 && 
                slice.reco.pfp[ipfp].parent_is_primary
            )  return 2;

            // Skip low energy protons
            if (chi2_values.proton < 100)
                start_mom_V3.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.z
                );

            if (
                chi2_values.proton < 100            && 
                (rec_vtx - start).Mag() < dist_cut  && 
                sqrt ( pow(particle_data::masses::proton, 2) + pow(start_mom_V3.Mag() * 1000, 2) ) - particle_data::masses::proton >= 50.0 && 
                slice.reco.pfp[ipfp].parent_is_primary
            ) return 1;

        }

        if (slice.reco.pfp[ipfp].trackScore < 0.5) {
            if (slice.reco.pfp[ipfp].trackScore >= 0.4 && chi2_values.proton < 100) {
                TVector3 start_mom_V3_2;
                start_mom_V3_2.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.z
                );

                if (
                    sqrt(pow(particle_data::masses::proton 2) + pow(start_mom_V3_2.Mag() * 1000, 2)) - particle_data::masses::proton >= 50.0 && 
                    (rec_vtx - start).Mag() < dist_cut && 
                    slice.reco.pfp[ipfp].parent_is_primary
                ) return 1;
            }

            if (!(slice.reco.pfp[ipfp].trackScore >= 0.4 && chi2_values.proton < 100)) {

                // int use_plane2 = slice.reco.pfp[ipfp].trk.calo[2].nhit>slice.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
                int use_plane2 = 2;

                if (std::isnan(slice.reco.pfp[ipfp].shw.plane[use_plane2].energy))
                    return 9;
                if (slice.reco.pfp[ipfp].shw.plane[use_plane2].energy * 1000 < 25.0)
                    return 9;
                if (slice.reco.pfp[ipfp].shw.plane[use_plane2].energy * 1000 > 25.0 && slice.reco.pfp[ipfp].parent_is_primary)
                    return 3;
            }
        }
        return 9;
    } // int id_pfp
} // namespace var_utils

namespace vars {
    namespace truth {
        const ana::Var slice_neutrino_energy ([](const caf::SRSliceProxy *slice) -> double {
            /* This slice is the selected slice (so reco 1µNp or truth 1µNp)
             * Here only the computation should be performed, since the cuts will be applied
             * at the Tree/Spectrum stage
            */
            return slice->truth.E // true neutrino energy in GeV
        });
    } // namespace truth

    namespace reco {
        const ana::Var slice_neutrino_energy_1muNp ([](const caf::SRSliceProxy *slice) -> double {
            /* This slice is the selected slice (so reco 1µNp or truth 1µNp)
             * Here only the computation should be performed, since the cuts will be applied
             * at the Tree/Spectrum stage
            */
        }); // const ana::Var slice_neutrino_energy_reco_1muNp
    } // namespace reco
} // namespace vars

#endif