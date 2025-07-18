
#ifndef PDG_HELPER_H
#define PDG_HELPER_H

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TMath.h"

#include "helper.h"
#include "selection.h"

#include <iostream>
#include <string>
#include <vector>

namespace vars {
    namespace pdg {

        var_utils::chi2 get_chi2 (const caf::SRSliceProxy* slice, const std::size_t& iPfp, const int& usePlane = 2) {
            // compute new chi2
            std::vector<double> dedx;
            std::vector<double> rr;
            for (std::size_t iHit = 0; iHit < slice->reco.pfp[iPfp].trk.calo[usePlane].points.size(); ++iHit) {
                dedx.push_back(slice->reco.pfp[iPfp].trk.calo[usePlane].points[iHit].dedx);
                rr.push_back(slice->reco.pfp[iPfp].trk.calo[usePlane].points[iHit].rr);
            } // calo points

            // input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
            // output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt}
            return var_utils::chi2_ALG(dedx, rr, 0.0, 25.0);
        }

        bool zero_vtXDist (const caf::SRSliceProxy* slice, const std::size_t& iPfp) {
            // std::cout << "Vertex in (" 
            //           << slice->truth.position.x << ", " 
            //           << slice->truth.position.y << ", " 
            //           << slice->truth.position.z << "), startpoint in (" 
            //           << slice->reco.pfp[iPfp].trk.truth.p.start.x << ", " 
            //           << slice->reco.pfp[iPfp].trk.truth.p.start.y << ", " 
            //           << slice->reco.pfp[iPfp].trk.truth.p.start.z << "), delta (" 
            //           << slice->truth.position.x -  slice->reco.pfp[iPfp].trk.truth.p.start.x << ", " 
            //           << slice->truth.position.y -  slice->reco.pfp[iPfp].trk.truth.p.start.y << ", "  
            //           << slice->truth.position.z -  slice->reco.pfp[iPfp].trk.truth.p.start.z << ")" 
            //           << std::endl; 
            return 0 == (
                std::pow(slice->truth.position.x - slice->reco.pfp[iPfp].trk.truth.p.start.x, 2) + 
                std::pow(slice->truth.position.y - slice->reco.pfp[iPfp].trk.truth.p.start.y, 2) + 
                std::pow(slice->truth.position.z - slice->reco.pfp[iPfp].trk.truth.p.start.z, 2)
            );
        }
        
        const ana::Var slice_muon_chi2_mu ([](const caf::SRSliceProxy* slice) -> double {
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) {
                    return get_chi2(slice, iPfp).muon;
                }
                iPfp++;
            }
            return -9999;
        });

        const ana::Var slice_muon_chi2_p ([](const caf::SRSliceProxy* slice) -> double {
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) {
                    return get_chi2(slice, iPfp).proton;
                }
                iPfp++;
            }
            return -9999;
        });

        const ana::Var slice_muon_length ([](const caf::SRSliceProxy* slice) -> double {
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) {
                    // std::cout << "This is 13 (µ) and primary, pfp.trk.truth.p.parent = " 
                    //           << pfp.trk.truth.p.parent << std::endl; ////////////////////////////////////// COMMENT THIS OUT JUST AFTER CHECK
                    return pfp.trk.len;
                }
                iPfp++;
            }
            return -9999;
        });

        const ana::Var slice_muon_length_ratio ([](const caf::SRSliceProxy* slice) -> double {
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) {
                    return pfp.trk.len / pfp.trk.truth.p.length;
                }
                iPfp++;
            }
            return -9999;
        });

        const ana::Var slice_muonPandoraPrimary ([](const caf::SRSliceProxy* slice) -> double {
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) {
                    return pfp.parent_is_primary;
                }
                iPfp++;
            }
            return -9999;
        });

        const ana::Var slice_muon_vtxDist ([](const caf::SRSliceProxy* slice) -> double {
            TVector3 reco_vtx;
            TVector3 reco_start;

            reco_vtx.SetXYZ(slice->vertex.x, slice->vertex.y, slice->vertex.z);
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) {
                    reco_start.SetXYZ(pfp.trk.start.x, pfp.trk.start.y, pfp.trk.start.z);
                }
                iPfp++;
            }

            return (reco_vtx - reco_start).Mag();
        });

        const ana::Var slice_muon_trackScore ([](const caf::SRSliceProxy* slice) -> double {
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) 
                    return pfp.trackScore;
                iPfp++;
            }
            return -9999;
        });

        const ana::Var slice_muon_completeness ([](const caf::SRSliceProxy* slice) -> double {
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) 
                    return pfp.trk.truth.bestmatch.hit_completeness;
                iPfp++;
            }
            return -9999;
        });

        const ana::Var slice_muon_purity ([](const caf::SRSliceProxy* slice) -> double {
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 13 && zero_vtXDist(slice, iPfp)) 
                    return pfp.trk.truth.bestmatch.hit_purity;
                iPfp++;
            }
            return -9999;
        });

    const ana::MultiVar slice_proton_chi2_mu ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> chis;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) {
                    chis.push_back(get_chi2(slice, iPfp).muon);
                }
                iPfp++;
            }
            return chis;
        });

        const ana::MultiVar slice_proton_chi2_p ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> chis;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) {
                    chis.push_back(get_chi2(slice, iPfp).proton);
                }
                iPfp++;
            }
            return chis;
        });

        const ana::MultiVar slice_proton_length ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> lengths;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) {
                    lengths.push_back(pfp.trk.len);
                }
                iPfp++;
            }
            return lengths;
        });

        const ana::MultiVar slice_proton_length_ratio ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> lengths;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) {
                    lengths.push_back(pfp.trk.len / pfp.trk.truth.p.length);
                }
                iPfp++;
            }
            return lengths;
        });


        const ana::MultiVar slice_proton_vtxDist ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            TVector3 reco_vtx;
            TVector3 reco_start;
            std::vector<double> vtxDist;

            reco_vtx.SetXYZ(slice->vertex.x, slice->vertex.y, slice->vertex.z);
            std::size_t iPfp = 0;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) {
                    reco_start.SetXYZ(pfp.trk.start.x, pfp.trk.start.y, pfp.trk.start.z);
                    vtxDist.push_back((reco_vtx - reco_start).Mag());
                }
                iPfp++;
            }

            return vtxDist;
        });

        const ana::MultiVar slice_proton_trackScore ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> trackScore;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) 
                    trackScore.push_back(pfp.trackScore);
                iPfp++;
            }
            return trackScore;
        });

        const ana::MultiVar slice_proton_completeness ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> hits;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) 
                    hits.push_back(pfp.trk.truth.bestmatch.hit_completeness);
                iPfp++;
            }
            return hits;
        });

        const ana::MultiVar slice_proton_purity ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> hits;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) 
                    hits.push_back(pfp.trk.truth.bestmatch.hit_purity);
                iPfp++;
            }
            return hits;
        });

        const ana::MultiVar slice_proton_depEnergy ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> depEnergy;
            TVector3 startMomentum;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) {
                    startMomentum.SetXYZ(
                        pfp.trk.rangeP.p_proton * pfp.trk.dir.x, 
                        pfp.trk.rangeP.p_proton * pfp.trk.dir.y, 
                        pfp.trk.rangeP.p_proton * pfp.trk.dir.z
                    );
                    depEnergy.push_back(
                        std::sqrt ( 
                            std::pow(particle_data::masses::proton, 2) + std::pow(startMomentum.Mag() * particle_data::GeV, 2) 
                        ) - particle_data::masses::proton
                    );
                } 
                iPfp++;
            }
            return depEnergy;
        });

        const ana::MultiVar slice_protonPandoraPrimary ([](const caf::SRSliceProxy* slice) -> std::vector<double> {
            std::size_t iPfp = 0;
            std::vector<double> primTag;
            for (auto const& pfp: slice->reco.pfp) {
                if (std::abs(pfp.trk.truth.p.pdg) == 2212 && zero_vtXDist(slice, iPfp)) {
                    primTag.push_back(pfp.parent_is_primary);
                }
                iPfp++;
            }
            return primTag;
        });

    }
} 

#endif