#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h" //after v09_44 release
//#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/Binnings.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NumuCuts.h"
#include "TVector3.h"


#include <algorithm>
#include <vector>

#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "stdio.h"
#include "TProfile.h"
using namespace ana;

const double mass_muon = 105.6583755; //105.658;
const double mass_proton = 938.27208816; // 938.3
const double mass_pion = 139.57039; //139.570

//ofstream MyFile22("signal_mc.txt");
//ofstream MyFile("Debug.txt");
//ofstream MyFile("Selected_9435.txt");
//ofstream MyFile("test_spacepoints_11982.txt");
//ofstream MyFile("Optim_pscore_prescaled.txt");
//ofstream MyFile2("counting_shortsignal.txt");

//ofstream MyFile("Output_tpcind2transparent.txt");
//ofstream MyFile("Test_MC_reco.txt");
//ofstream MyFile2("output_pT_class_true.txt");

//ofstream MyFile2("new_signal_mc_CH.txt");
//ofstream MyFile3("signal_mc_ML.txt");


//ofstream MyFile("test_taskforce_30_1muNp_MC.txt");
//ofstream MyFile("test_taskforce_30_1muNp_DATA_20cm.txt");
//ofstream MyFile("test_taskforce_30_1muNp_DATA_10cm_EW_WE.txt");
//ofstream MyFile2("test_taskforce_30_loose_DATA_10cm.txt");


TFile* file = TFile::Open(
    "/exp/icarus/app/users/msotgia/analysis/sbnana_v09_93_01_thesis_analysis/analysis/dEdxrestemplates.root"
);

//TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/marterop/sbnana_v09_78_06/mc_test/dEdxrestemplates.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");



//Data
//double min_time =  -0.4;
//double max_time = 1.5;

//MC
//double min_time =  0.0;
//double max_time = 1.6;



//Daniele's function to compute trigger efficiency

int npoints=6;
double e[6]={0,100,200,300,400,500};
double eff[6]={0,.56,.80,.98,1.,1.};

double TriggerEff(double Edep){
  // Edep is the deposited energy in MeV
  if(Edep<=0) return eff[0];
  if(Edep>e[npoints-1]) return eff[npoints-1];
  for(int i=0;i<npoints-1;i++) {
     if(Edep>=e[i]&&Edep<=e[i+1]) {
       return eff[i]+(eff[i+1]-eff[i])/(e[i+1]-e[i])*(Edep-e[i]);
     }
  }
  return 0; // should never happen
}

std::vector<double> chi2_ALG(std::vector<double> &dEdx,std::vector<double> &RR, double rr_min, double rr_max)
{
    //The output is chi2s



    double threshold = 0.5;
    double max_rr = rr_max;
    double min_rr = rr_min;

    std::vector<float> trkdedx;
    std::vector<float> trkres;
    std::vector<double> vpida;
    

    for(std::size_t i(0); i<dEdx.size(); ++i){
      if(i==0 || i==dEdx.size()-1)continue;
        if(RR[i]<max_rr && RR[i]>rr_min ){trkdedx.push_back(dEdx[i]);trkres.push_back(RR[i]);}       
    }


    int npt = 0;
    double chi2pro = 0;
    double chi2ka = 0;
    double chi2pi = 0;
    double chi2mu = 0;
    double avgdedx = 0;
    double PIDA = 0;

    int used_trkres = 0;
    for (unsigned i = 0; i < trkdedx.size(); ++i) { //hits
      //ignore the first and the last point
      //if (i == 0 || i == trkdedx.size() - 1) continue;
      avgdedx += trkdedx[i];
      if (trkres[i] < 26) {
        PIDA += trkdedx[i] * pow(trkres[i], 0.42);
        vpida.push_back(trkdedx[i] * pow(trkres[i], 0.42));
        used_trkres++;
      }
      if (trkdedx[i] > 100 || trkdedx[i]<threshold) continue; //protect against large pulse height
      int bin = dedx_range_pro->FindBin(trkres[i]);
      if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
        double bincpro = dedx_range_pro->GetBinContent(bin);
        if (bincpro < 1e-6) { //for 0 bin content, using neighboring bins
          bincpro =
            (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
        }
        double bincka = dedx_range_ka->GetBinContent(bin);
        if (bincka < 1e-6) {
          bincka =
            (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;
        }
        double bincpi = dedx_range_pi->GetBinContent(bin);
        if (bincpi < 1e-6) {
          bincpi =
            (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
        }
        double bincmu = dedx_range_mu->GetBinContent(bin);
        if (bincmu < 1e-6) {
          bincmu =
            (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
        }
        double binepro = dedx_range_pro->GetBinError(bin);
        if (binepro < 1e-6) {
          binepro =
            (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;
        }
        double bineka = dedx_range_ka->GetBinError(bin);
        if (bineka < 1e-6) {
          bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;
        }
        double binepi = dedx_range_pi->GetBinError(bin);
        if (binepi < 1e-6) {
          binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;
        }
        double binemu = dedx_range_mu->GetBinError(bin);
        if (binemu < 1e-6) {
          binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;
        }
        //double errke = 0.05*trkdedx[i];   //5% KE resolution
        double errdedx = 0.04231 + 0.0001783 * trkdedx[i] * trkdedx[i]; //resolution on dE/dx
        errdedx *= trkdedx[i];
        chi2pro += pow((trkdedx[i] - bincpro) / std::sqrt(pow(binepro, 2) + pow(errdedx, 2)), 2);
        chi2ka += pow((trkdedx[i] - bincka) / std::sqrt(pow(bineka, 2) + pow(errdedx, 2)), 2);
        chi2pi += pow((trkdedx[i] - bincpi) / std::sqrt(pow(binepi, 2) + pow(errdedx, 2)), 2);
        chi2mu += pow((trkdedx[i] - bincmu) / std::sqrt(pow(binemu, 2) + pow(errdedx, 2)), 2);
        //std::cout<<i<<" "<<trkdedx[i]<<" "<<trkres[i]<<" "<<bincpro<<std::endl;
        ++npt;
      }
    } //hits
        std::vector<double> chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt};

    return chi2s;
}



bool isInFV (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
  //need to add a check to avoid having the vtx slice in one cryo and the end track in the other cryo

    //fiducial volume for the dangling cable 
  if(x>210.0 && y > 60.0 && z> 290.0 && z< 390.0 ) return false;
  //Cathode padding
  //if(std::abs(x)>208.79 && std::abs(x)<211.79) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
			( x >  61.94 + 25 && x <  358.49 - 25 )) &&
		  ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
		  ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
}

bool isInContained (double x, double y, double z, double dist)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;


//5 cm containment for Y in both cryo
  return (( ( x < -61.94 - dist && x > -358.49 + dist ) ||
			( x >  61.94 + dist && x <  358.49 - dist )) &&
		  ( ( y > -181.86 + dist && y < 134.96 - dist ) &&
		  ( z > -894.95 + dist && z < 894.95 - dist ) ));

/*
//10 cm containment for Y in both cryo
  return (( ( x < -61.94 - dist && x > -358.49 + dist ) ||
			( x >  61.94 + dist && x <  358.49 - dist )) &&
		  ( ( y > -181.86 + 10 && y < 134.96 - 10 ) &&
		  ( z > -894.95 + dist && z < 894.95 - dist ) ));
*/

}


bool isInActive (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 && x > -358.49 ) ||
			( x >  61.94 && x <  358.49)) &&
		  ( ( y > -181.86 && y < 134.96)  &&
		    ( z > -894.95 && z < 894.95) ));
}

bool isInContained (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
  return (( ( x < -61.94 - 5 && x > -358.49 + 5 ) ||
			( x >  61.94 + 5 && x <  358.49 - 5 )) &&
		  ( ( y > -181.86 + 5 && y < 134.96 - 5 ) &&
		  ( z > -894.95 + 5 && z < 894.95 - 5 ) ));

}

bool isInDetector (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
  return (( ( x < -61.94 + 5 && x > -358.49 - 5 ) ||
			( x >  61.94 - 5 && x <  358.49 + 5 )) &&
		  ( ( y > -181.86 - 5 && y < 134.96 + 5 ) &&
		  ( z > -894.95 - 5 && z < 894.95 + 5 ) ));

}


//functions for automatic selection 

bool all_contained ( const caf::Proxy<caf::SRSlice>& islc ) { 

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;
    //check meaningful points 
    //if(!isInDetector(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z))continue;
    //if(!isInDetector(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z))continue;
    //if (!(islc.reco.pfp[ipfp].parent_is_primary ))continue; //skip secondaries
    //if(islc.reco.pfp[ipfp].trackScore<0.4)continue; //Want to check only tracks??
    if((islc.reco.pfp[ipfp].trk.start.x*islc.vertex.x)<0){return false;} //not contained if they cross cryostats
    //if not contained return false
    if(!isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0)){return false;}
      
    }   
    
return true;
}


bool all_contained_mc ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRTrueInteraction>& nu ) { 
    //Check only those pfp that are visible 
                
                for ( auto const& ipart : nu.prim ){
                    if ( ipart.G4ID < 0 )  continue;
                    if ( ipart.cryostat < 0 )  continue;
                    int iG4ID_parent;
                    //check if charged primaries are contained: 
                    if(abs(ipart.pdg)==13 || abs(ipart.pdg)==2212 || abs(ipart.pdg)==211 || abs(ipart.pdg)==11){
                        if(isInContained(ipart.end.x,ipart.end.y,ipart.end.z)==false){return false;}
                    }                   
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){
                        if(abs(itrue.pdg)==13 || abs(itrue.pdg)==2212 || abs(itrue.pdg)==211 || abs(itrue.pdg)==11){
                        if(isInContained(itrue.end.x,itrue.end.y,itrue.end.z)==false){return false;}
                    }                             
                        }
                        }
                    
                        } 
 
                }  
return true;
} 

bool all_contained_truth ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc ) { 
    //Check only those pfp that are visible 
                
                for ( auto const& ipart : islc.truth.prim ){
                    if ( ipart.G4ID < 0 )  continue;
                    if ( ipart.cryostat < 0 )  continue;
                    int iG4ID_parent;
                    double dep_E=0;   
                    //check if charged primaries are contained: 
                    if(abs(ipart.pdg)==13 || abs(ipart.pdg)==2212 || abs(ipart.pdg)==211 || abs(ipart.pdg)==11){
                        if(isInContained(ipart.end.x,ipart.end.y,ipart.end.z)==false){return false;}
                    }                   
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){
                        if(abs(itrue.pdg)==13 || abs(itrue.pdg)==2212 || abs(itrue.pdg)==211 || abs(itrue.pdg)==11){
                        if(isInContained(itrue.end.x,itrue.end.y,itrue.end.z)==false){return false;}
                    }                             
                        }
                        }
                    
                        } 
 
                }  
return true;
} 

int find_muon ( const caf::Proxy<caf::SRSlice>& islc, int dist_mucut) { 

    //Select muon as longest track
    double max_length=-1.0;
    int ipfp_mu = -1;
    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 RecoStart;
        for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
        if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;
    //if(!isInDetector(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z))continue;
    //if(!isInDetector(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z))continue;
        RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
        if(islc.reco.pfp[ipfp].trackScore<0.5)continue;
        //if(islc.reco.pfp[ipfp].trackScore<0.4)continue;
   
  //int use_plane = islc.reco.pfp[ipfp].trk.calo[2].nhit>islc.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
    int use_plane = 2;
    //compute new chi2
    std::vector<double> output;
    std::vector<double> dedx;
    std::vector<double> rr;
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit ){
                        dedx.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
                        rr.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
                            } // calo points
                  //input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
//output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt} 
    output = chi2_ALG(dedx,rr,0.0,25.0);


        if(islc.reco.pfp[ipfp].trk.len>max_length && ((RecoVtx-RecoStart).Mag()<dist_mucut) && islc.reco.pfp[ipfp].trk.len>50 && 
        output[0]<30 && output[1]>60 && 
        isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) && 
        (islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0 && islc.reco.pfp[ipfp].parent_is_primary){
        max_length=islc.reco.pfp[ipfp].trk.len;
        ipfp_mu=ipfp;
            }
        }//loop of pfp to find muon
return ipfp_mu;
}

int id_pfp ( const caf::Proxy<caf::SRSlice>& islc, int ipfp, int dist_cut ) { 
    //return 1 PROTONS
    //return 2 PIONS
    //return 3 SHOWER
    //return 9 other -> nan, not primary, too far, below energy threshold... 


    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
//    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
    //skip secondaries
    if (!(islc.reco.pfp[ipfp].parent_is_primary ))return 9;
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;
    //if(int(ipfp)==ipfp_mu)continue;     //There is always a muon, for a 1mu1p we need 2 tracks - 1 muon = 1 only proton with threshold
    //consider only primary tracks which are 20cm close to the vertex, either vtx-start or vtx-end
    TVector3 start(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
    TVector3 end(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z);
    //condition ? result_if_true : result_if_false
    double min_dist = ((start-RecoVtx).Mag()< (end-RecoVtx).Mag() ? (start-RecoVtx).Mag() : (end-RecoVtx).Mag());
    //if(min_dist>50.0)continue;
    //if(min_dist>10.0) return 9;
    if(min_dist>50.0) return 9;
     

    //int use_plane = islc.reco.pfp[ipfp].trk.calo[2].nhit>islc.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
    int use_plane = 2;
    //compute new chi2
    std::vector<double> output;
    std::vector<double> dedx;
    std::vector<double> rr;
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit ){
                        dedx.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
                        rr.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
                            } // calo points
                  //input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
//output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt} 
    output = chi2_ALG(dedx,rr,0.0,25.0);
    if(islc.reco.pfp[ipfp].trackScore>=0.5){
    if (std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;
    if (std::isnan(islc.reco.pfp[ipfp].trk.start.y) || std::isnan(islc.reco.pfp[ipfp].trk.start.z)|| std::isnan(islc.reco.pfp[ipfp].trk.end.y) || std::isnan(islc.reco.pfp[ipfp].trk.end.z) )return 9;
    

    //skip low energy tagged pions
    TVector3 Start_mom_v;
    if(output[1]>=100 ){ Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);}
    if(output[1]>=100 && ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(mass_pion,2)+pow(Start_mom_v.Mag()*1000,2))-mass_pion)>=25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 2;}}
    //skip low energy protons
    if(output[1]<100 ){ Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);}
    if(output[1]<100 && ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(mass_proton,2)+pow(Start_mom_v.Mag()*1000,2))-mass_proton)>=50.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 1;}}
            
    }
    if(islc.reco.pfp[ipfp].trackScore<0.5){
        if(islc.reco.pfp[ipfp].trackScore>=0.4 && islc.reco.pfp[ipfp].trackScore<0.5 && output[1]<100 ){
            TVector3 Start_mom_v2;
            Start_mom_v2.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
            if((sqrt(pow(mass_proton,2)+pow(Start_mom_v2.Mag()*1000,2))-mass_proton)>=50.0 && ((RecoVtx-start).Mag()<dist_cut) && islc.reco.pfp[ipfp].parent_is_primary){return 1;}
        } 
    if(!(islc.reco.pfp[ipfp].trackScore>=0.4 && islc.reco.pfp[ipfp].trackScore<0.5 && output[1]<100 )){
    //int use_plane2 = islc.reco.pfp[ipfp].trk.calo[2].nhit>islc.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;  
    int use_plane2 = 2;
    if(std::isnan(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy))return 9;
    if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000<25.0)return 9;
    if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000>25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 3;}}}
        }
    
    return 9;
}

//reco energy with fermi momentum = binding energy 

double Neutrino_energy_reco ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) {
    
    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;
    float E_mu,E_p;    
    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV
    E_mu=sqrt(p_mu_tot*p_mu_tot+(mass_muon*mass_muon)/(1000*1000));

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    double p_p_tot = sqrt(p_p_x*p_p_x+p_p_y*p_p_y+p_p_z*p_p_z);
    E_p=sqrt(p_p_tot*p_p_tot+(938.272*938.272)/(1000*1000));

    //Total momentum
    p_tot_x=p_p_x+p_mu_x;
    p_tot_y=p_p_y+p_mu_y;
    p_tot_z=p_p_z+p_mu_z;

    double R = 37215.516+1000*(p_mu_z+p_p_z-E_p-E_mu);  //MeV
    double E_b=21.8; //MeV Binding energy
    double p_T=1000*sqrt(p_tot_x*p_tot_x+p_tot_y*p_tot_y);
    double p_L=R/2-p_T*p_T/(2*R)-(37215.516-939.565+E_b)*(37215.516-939.565+E_b)/(2*R);
    double E_nu = 1000*(p_mu_z+p_p_z)-p_L;

return E_nu;
}

double Neutrino_energy_reco_Np ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int dist_emucut) {
    

        float p_mu_x = -1; float p_mu_y=-1; float p_mu_z=-1;
        float p_p_x=-1; float p_p_y =-1; float p_p_z=-1;
        float p_tot_x =-1; float p_tot_y=-1; float p_tot_z=-1;
        double E_mu =0 ; double E_p=0; 
        int ipfp_pro = -1; TVector3 Start_mom;

        p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
        p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
        p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
        double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV
        E_mu=1000*sqrt(p_mu_tot*p_mu_tot+(mass_muon*mass_muon)/(1000*1000));
        
        for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
            if(int(ipfp)==ipfp_mu)continue;
            if(id_pfp(islc, ipfp,dist_emucut)==1){
            TVector3 Start_mom_v2;
            Start_mom_v2.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
            E_p += (sqrt(pow(mass_proton,2)+pow(Start_mom_v2.Mag()*1000,2))-mass_proton);
            ipfp_pro=ipfp;
            }
                    }

       double E_nu = (E_mu + E_p)/1000;




return E_nu;
}

double Transverse_mom_reco ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) {
    
    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;
    float E_mu,E_p;    
    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV
    E_mu=sqrt(p_mu_tot*p_mu_tot+(mass_muon*mass_muon)/(1000*1000));

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    double p_p_tot = sqrt(p_p_x*p_p_x+p_p_y*p_p_y+p_p_z*p_p_z);
    E_p=sqrt(p_p_tot*p_p_tot+(938.272*938.272)/(1000*1000));

    //Total momentum
    p_tot_x=p_p_x+p_mu_x;
    p_tot_y=p_p_y+p_mu_y;
    p_tot_z=p_p_z+p_mu_z;

    double p_T=1000*sqrt(p_tot_x*p_tot_x+p_tot_y*p_tot_y);


return p_T;
}

double Transverse_mom_reco_Np ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu) {

        
        float p_p_x=0; float p_p_y =0; float p_p_z=0;
        float p_tot_x =0; float p_tot_y=0; float p_tot_z=0;
        int ipfp_pro = -1; 

        float p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
        float p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
        float p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;

        for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
            if(int(ipfp)==ipfp_mu)continue;
            if(id_pfp(islc, ipfp,10)==1){
                p_p_x +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x;
                p_p_y +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y;
                p_p_z +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z;
                ipfp_pro=ipfp;

            }
                    }

    if(ipfp_mu!=-1 && ipfp_pro!=-1){
            p_tot_x=p_p_x+p_mu_x;
            p_tot_y=p_p_y+p_mu_y;
            p_tot_z=p_p_z+p_mu_z;
                      
    }
        
       double p_T = (sqrt(p_tot_x*p_tot_x+p_tot_y*p_tot_y));


return p_T;
}


double T3D_angle_mup ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) { 

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    TVector3 mu_vector(p_mu_x,p_mu_y,p_mu_z);

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    TVector3 pro_vector(p_p_x,p_p_y,p_p_z);

     
    return cos(mu_vector.Angle(pro_vector)); 
    }

double T3D_angle_mup_true ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) { 

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=islc.reco.pfp[ipfp_mu].trk.truth.p.startp.x; //GeV
    p_mu_y=islc.reco.pfp[ipfp_mu].trk.truth.p.startp.y; //GeV
    p_mu_z=islc.reco.pfp[ipfp_mu].trk.truth.p.startp.z; //GeV
    TVector3 mu_vector(p_mu_x,p_mu_y,p_mu_z);

    //Proton momentum
    p_p_x=islc.reco.pfp[ipfp_pro].trk.truth.p.startp.x; //GeV
    p_p_y=islc.reco.pfp[ipfp_pro].trk.truth.p.startp.y; //GeV
    p_p_z=islc.reco.pfp[ipfp_pro].trk.truth.p.startp.z; //GeV
    TVector3 pro_vector(p_p_x,p_p_y,p_p_z);

    
    return cos(mu_vector.Angle(pro_vector)); 
    }

double Transverse_angle ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) { 

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    double p_p_tot = sqrt(p_p_x*p_p_x+p_p_y*p_p_y+p_p_z*p_p_z);

    //Total momentum
    p_tot_x=p_p_x+p_mu_x;
    p_tot_y=p_p_y+p_mu_y;
    p_tot_z=p_p_z+p_mu_z;

    //TRANSVERSE PLANE - ANGLE! 
    double norm_mu=sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y);
    double norm_pro=sqrt(p_p_x*p_p_x+p_p_y*p_p_y);    
       
     
    return (p_mu_x*p_p_x+p_mu_y*p_p_y)/(norm_mu*norm_pro); 
    }





//bool automatic_selection ( const caf::Proxy<caf::SRSlice>& islc ) {

bool automatic_selection ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc ,int dist_smucut  ){
         
        int ipfp_mu = -1;
        int ipfp_pro = -1;

        {
        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/)){
        {            

        {
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

        if(islc.barycenterFM.deltaZ_Trigger < 100 && islc.barycenterFM.deltaZ_Trigger>0)    
    {
            if(all_contained(islc)){
                ipfp_mu=find_muon(islc,dist_smucut);
                if(ipfp_mu!=-1){
                    int num_protons = 0;
                    int num_pions   = 0;
                    int num_showers = 0; 

                    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
                        if(int(ipfp)==ipfp_mu)continue;
                        //if(!isInDetector(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z))continue;
                        //if(!isInDetector(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z))continue;
                        if(id_pfp(islc, ipfp,dist_smucut)==1){num_protons+=1;}
                        if(id_pfp(islc, ipfp,dist_smucut)==2){num_pions+=1;}
                        if(id_pfp(islc, ipfp,dist_smucut)==3){num_showers+=1;}
                    }
                    if(num_protons==1 && num_pions==0 && num_showers==0){
                        return true;
                            }//1mu1p 
                        }//muon with conditions found
                    
                    }//all tracks of slice contained  
        }//new Barycenter match
    }//fiducial condition
 
    }//signal
    }//only neutrinos in active!
    }//check no nan in true info
    }//islc.truth.index>= 0

   
return false;
}

/*
int automatic_selection_mu_index ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_smucut , int cut_baryc) {
     
         
        int ipfp_mu = -1;
        int ipfp_pro = -1;

        {
        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z))){
        {            

        {
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

        if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0 )    
    {
            if(all_contained(islc)){
                ipfp_mu=find_muon(islc,dist_smucut);
                if(ipfp_mu!=-1){return ipfp_mu;}//all tracks of slice contained 
                    
                    }//muon with conditions found
        }//new Barycenter match
    }//fiducial condition
 
    }//signal
    }//only neutrinos in active!
    }//check no nan in true info
    }//islc.truth.index>= 0

   
return -1;
}
*/

         
int automatic_selection_mu_index ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_smucut , int cut_baryc) {

        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/)){
        {            

        {
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

        if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0)   
    {
            if(all_contained(islc))
            {
                ipfp_mu=find_muon(islc,dist_smucut);
                if(ipfp_mu!=-1){
                    int num_protons =0;
                    int num_pions =0;
                    int num_showers =0; 

                    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
                        if(int(ipfp)==ipfp_mu)continue;
                        //if(!isInDetector(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z))continue;
                        //if(!isInDetector(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z))continue;                        
                        if(id_pfp(islc, ipfp,dist_smucut)==1){num_protons+=1;}
                        if(id_pfp(islc, ipfp,dist_smucut)==2){num_pions+=1;}
                        if(id_pfp(islc, ipfp,dist_smucut)==3){num_showers+=1;}
                    }
                    if(num_protons>0 && num_pions==0 && num_showers==0){
                        return ipfp_mu;
                            }//1muNp 
                        }//muon with conditions found
                    
                    }//all tracks of slice contained  
        }//new Barycenter match
    }//fiducial condition
 
    }//signal
    }//only neutrinos in active!
    }//check no nan in true info

   
return false;
}

int automatic_selection_p_index ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc , int dist_smucut) {
         

        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/)){
        {            

        {
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

        if(islc.barycenterFM.deltaZ_Trigger < 100 && islc.barycenterFM.deltaZ_Trigger>0)    
    {
            if(all_contained(islc)){
                ipfp_mu=find_muon(islc,dist_smucut);
                if(ipfp_mu!=-1){
                    int num_protons = 0;
                    int num_pions   = 0;
                    int num_showers = 0; 

                    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
                        if(int(ipfp)==ipfp_mu)continue;
                        //if(!isInDetector(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z))continue;
                        //if(!isInDetector(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z))continue;
                        if(id_pfp(islc, ipfp,dist_smucut)==1){num_protons+=1; ipfp_pro=ipfp;}
                        if(id_pfp(islc, ipfp,dist_smucut)==2){num_pions+=1;}
                        if(id_pfp(islc, ipfp,dist_smucut)==3){num_showers+=1;}
                    }
                    if(num_protons==1 && num_pions==0 && num_showers==0){
                        return ipfp_pro;
                            }//1mu1p 
                        }//all tracks of slice contained 
                    
                    }//muon with conditions found
        }//new Barycenter match
    }//fiducial condition
 
    }//signal
    }//only neutrinos in active!
    }//check no nan in true info

   
return -1;
}

bool automatic_selection_1muNp ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_cut, int cut_baryc ){

        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/)){
        {            

        {
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

        if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0)   
    {
            if(all_contained(islc))
            {
                ipfp_mu=find_muon(islc,dist_cut);
                if(ipfp_mu!=-1){
                    int num_protons =0;
                    int num_pions =0;
                    int num_showers =0; 

                    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
                        if(int(ipfp)==ipfp_mu)continue;
                        //if(!isInDetector(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z))continue;
                        //if(!isInDetector(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z))continue;
                        if(id_pfp(islc, ipfp,dist_cut)==1){num_protons+=1;}
                        if(id_pfp(islc, ipfp,dist_cut)==2){num_pions+=1;}
                        if(id_pfp(islc, ipfp,dist_cut)==3){num_showers+=1;}
                    }
                    if(num_protons>0 && num_pions==0 && num_showers==0){
                        return true;
                            }//1mu1p 
                        }//muon with conditions found
                    
                    }//all tracks of slice contained  
        }//new Barycenter match
    }//fiducial condition
 
    }//signal
    }//only neutrinos in active!
    }//check no nan in true info

   
return false;
}

bool automatic_selection_1mu1p ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_cut, int cut_baryc ){

    if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/)){
    {            

    {
    int ipfp_mu = -1;
    int ipfp_pro = -1;
if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

    if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0)   
{
        if(all_contained(islc))
        {
            ipfp_mu=find_muon(islc,dist_cut);
            if(ipfp_mu!=-1){
                int num_protons =0;
                int num_pions =0;
                int num_showers =0; 

                for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
                    if(int(ipfp)==ipfp_mu)continue;
                    //if(!isInDetector(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z))continue;
                    //if(!isInDetector(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z))continue;
                    if(id_pfp(islc, ipfp,dist_cut)==1){num_protons+=1;}
                    if(id_pfp(islc, ipfp,dist_cut)==2){num_pions+=1;}
                    if(id_pfp(islc, ipfp,dist_cut)==3){num_showers+=1;}
                }
                if(num_protons==1 && num_pions==0 && num_showers==0){
                    return true;
                        }//1mu1p 
                    }//muon with conditions found
                
                }//all tracks of slice contained  
    }//new Barycenter match
}//fiducial condition

}//signal
}//only neutrinos in active!
}//check no nan in true info


return false;
}


bool automatic_selection_nuCC ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_cut, int cut_baryc ){

        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/)){
        {            

        {
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

        if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0)   
    {
            if(all_contained(islc)){
                ipfp_mu=find_muon(islc,dist_cut);
                if(ipfp_mu!=-1){ return true;}//muon with conditions found
                    
                    }//all tracks of slice contained  
        }//new Barycenter match
    }//fiducial condition
 
    }//signal
    }//only neutrinos in active!
    }//check no nan in true info

   
return false;
}

int automatic_selection_mu_index_nuCC ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_smucut , int cut_baryc) {

        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/)){
        {            

        {
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

        if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0)   
    {
            if(all_contained(islc)){
                ipfp_mu=find_muon(islc,dist_smucut);
                if(ipfp_mu!=-1){return ipfp_mu;}//muon with conditions found
                    
                    }//all tracks of slice contained  
        }//new Barycenter match
    }//fiducial condition
 
    }//signal
    }//only neutrinos in active!
    }//check no nan in true info

   
return false;
}

int classification_type ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc ) { 

    // return 2 -> Visible and true 1mu1p
    // return 5 -> Visible and true 1muNp
    // return 3 -> Other 

    //CHECK WHICH TYPE OF INTERACTION 
    if(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)) 
        return 99;
    
    if(
        std::isnan(islc.truth.position.x) || 
        std::isnan(islc.truth.position.y) || 
        std::isnan(islc.truth.position.z)
    ) 
        return 99;
                
    if ( 
        abs(islc.truth.pdg) == 14 && 
        islc.truth.iscc && 
        !std::isnan(islc.truth.position.x) && 
        !std::isnan(islc.truth.position.y) && 
        !std::isnan(islc.truth.position.z) &&
		isInActive(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z) ){

            if(isInFV(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z)){

                int num_protons_above50 = 0;
                int num_muons = 0; 
                double length_muon = 0;                 
                double dep_E=0;
                
                for ( auto const& ipart : islc.truth.prim ){
                    if ( ipart.G4ID < 0 )  continue;

                    if(abs(ipart.pdg) == 13){
                        num_muons+=1; 
                        length_muon =ipart.length;
                    
                    }  // muons
                    
                    if(abs(ipart.pdg)==211){
                        return 3;
                    }  //Charged pions
                    int iG4ID_parent;  
                    //int use_plane = ipart.plane[ipart.cryostat][2].nhit>ipart.plane[ipart.cryostat][1].nhit ? 2:1;
                    int use_plane = 2;
                    
                    if(abs(ipart.pdg)==111){            //Neutral pions - reject if any of its gamma > 25 MeV
                        if(ipart.daughters.size()>0){
                            for ( auto const& itrue2 : sr->true_particles ){
                                iG4ID_parent=itrue2.parent;
                                //sum depE daughters 
                                if(iG4ID_parent==ipart.G4ID && abs(itrue2.pdg) == 22){
                                    if(itrue2.plane[ipart.cryostat][use_plane].visE*1000>25){
                                        return 3;
                                    };
                                }
                            }
                        }
                    }
                    if(abs(ipart.pdg) == 22){                    
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }}}
                        dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                        } 
                    if(abs(ipart.pdg)== 22 && dep_E>25.0){return 3;}   //primary photons   

                    if(abs(ipart.pdg)== 2212){                    
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }}}
                        dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                        }
                    if(abs(ipart.pdg)== 2212 && dep_E>50.0){num_protons_above50+=1;} //protons
                    dep_E=0;  
                    }
                if(num_muons == 1 && num_protons_above50 == 1){
                    if(length_muon > 50. && all_contained_truth(sr, islc)){return 2; }
                }
                if(num_muons == 1 && num_protons_above50 > 1){
                    if(length_muon > 50. && all_contained_truth(sr, islc)){return 5; }
                }
             }//fiducial 
             }//active
         
return 3;
}

int classification_type_MC ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRTrueInteraction>& nu ) { 

    // return 2 -> Visible and true 1mu1p
    // return 5 -> Visible and true 1muNp
    // return 3 -> Other 

            //CHECK WHICH TYPE OF INTERACTION 
if(std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z)) return 99;
if ( abs(nu.pdg) == 14 && nu.iscc && !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) &&
		 isInActive(nu.position.x,nu.position.y,nu.position.z) ){          
             if(isInFV(nu.position.x,nu.position.y,nu.position.z)){

                 int num_protons_above50 = 0;
                 int num_muons = 0; 
                 double length_muon2 = 0;                 
                double dep_E=0;
                for ( auto const& ipart : nu.prim ){
                    if ( ipart.G4ID < 0 )  continue;

                    if(abs(ipart.pdg) == 13){num_muons+=1; length_muon2 =ipart.length;}  // muons
                    if(abs(ipart.pdg)==211){return 3;}  //Charged pions
                    int iG4ID_parent;  
                    //int use_plane = ipart.plane[ipart.cryostat][2].nhit>ipart.plane[ipart.cryostat][1].nhit ? 2:1;
                    int use_plane = 2;
                    
                    if(abs(ipart.pdg)==111){            //Neutral pions - reject if any of its gamma > 25 MeV
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue2 : sr->true_particles ){
                        iG4ID_parent=itrue2.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID && abs(itrue2.pdg) == 22){if(itrue2.plane[ipart.cryostat][use_plane].visE*1000>25){return 3;};
                        }}}                        
                        }
                    if(abs(ipart.pdg) == 22){                    
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }}}
                        dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                        } 
                    if(abs(ipart.pdg)== 22 && dep_E>25.0){return 3;}   //primary photons   

                    if(abs(ipart.pdg)== 2212){                    
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }}}
                        dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                        }
                    if(abs(ipart.pdg)== 2212 && dep_E>50.0){num_protons_above50+=1;} //protons
                    dep_E=0;  
                    }
                if(num_muons == 1 && num_protons_above50 == 1){
                    if(length_muon2 > 50. && all_contained_mc(sr, nu)){return 2; }
                }
                if(num_muons == 1 && num_protons_above50 > 1){
                    if(length_muon2 > 50. && all_contained_mc(sr, nu)){return 5; }
                }
             }//fiducial 
             }//active
         
return 3;
}


const SpillVar kNu_energy_true_weight_cosmics([](const caf::SRSpillProxy* sr)-> double
{
    double weight = 1;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if(automatic_selection_1muNp(sr, islc,10,100 )){
        ipfp_mu = automatic_selection_mu_index(sr,islc,10,100);
    if(ipfp_mu!=-1){
    if(!(std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z) || std::isnan(islc.truth.E) || std::isnan(islc.reco.pfp[ipfp_mu].trk.truth.p.plane[islc.reco.pfp[ipfp_mu].trk.truth.p.cryostat][2].visE*1000))){
        weight = 1-0.36*std::sin((1.27*7.3*islc.truth.baseline)/(1000*islc.truth.E))*std::sin((1.27*7.3*islc.truth.baseline)/(1000*islc.truth.E));
        }
     if((std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z) || std::isnan(islc.truth.E) || std::isnan(islc.reco.pfp[ipfp_mu].trk.truth.p.plane[islc.reco.pfp[ipfp_mu].trk.truth.p.cryostat][2].visE*1000))){
        weight = 1.0;
     }
       
    }
    
    }//selected by the automatic selection
    
    }//loop over slices

 	return weight;
});

double flashtime = 0;

const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* spill){
  for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    double min_time =-1; double max_time =-1;
    if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    if(match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){
        flashtime = match.flashGateTime;
        /*
        if((spill->hdr.run == 10066 && spill->hdr.evt == 17465) ||  (spill->hdr.run == 9383 && spill->hdr.evt == 25629) || (spill->hdr.run == 9367 && spill->hdr.evt == 31831) ||(spill->hdr.run == 9558 && spill->hdr.evt == 21984) ||(spill->hdr.run == 9560 && spill->hdr.evt == 1502) ||(spill->hdr.run == 9796 && spill->hdr.evt == 43883) ||(spill->hdr.run == 9731 && spill->hdr.evt == 99026) ||(spill->hdr.run == 9837 && spill->hdr.evt == 3527)||(spill->hdr.run == 9844 && spill->hdr.evt == 16860) ||(spill->hdr.run == 9840 && spill->hdr.evt == 32953)||(spill->hdr.run == 9954 && spill->hdr.evt == 5378) ){
        cout << "flashGateTime " <<  match.flashGateTime << endl;
        }*/
          return true;} 

    //if(match.flashGateTime > -0.3 && match.flashGateTime < 1.3 && match.flashClassification == 0){  return true;} 



  }
  return false;
});

const SpillCut kGoodRuns([](const caf::SRSpillProxy* spill){

double good_runs[368] = {9301,9302,9303,9307,9308,9309,9310,9311,9312,9313,9314,9315,9316,9317,9318,9327,9328,9329,9330,9331,9332,9333,9334,9335,9336,9337,9338,9339,9340,9341,9342,9343,9344,9345,9346,9347,9353,9354,9355,9356,9357,9358,9359,9360,9361,9362,9363,9364,9365,9366,9367,9368,9369,9370,9371,9372,9373,9374,9375,9376,9377,9378,9379,9380,9383,9384,9385,9386,9387,9388,9389,9390,9391,9392,9393,9394,9409,9412,9417,9418,9419,9420,9421,9422,9423,9424,9425,9426,9427,9430,9431,9432,9435,9436,9437,9438,9439,9441,9444,9445,9448,9450,9457,9458,9460,9464,9472,9473,9474,9475,9477,9478,9481,9482,9495,9499,9504,9515,9516,9517,9518,9531,9533,9534,9558,9560,9562,9563,9564,9565,9568,9569,9570,9582,9583,9584,9585,9586,9587,9588,9589,9590,9593,9594,9595,9597,9598,9599,9602,9610,9624,9625,9626,9627,9631,9634,9642,9646,9647,9648,9649,9658,9672,9673,9674,9675,9687,9688,9689,9690,9691,9692,9693,9694,9695,9696,9697,9698,9699,9700,9703,9704,9705,9714,9715,9716,9717,9720,9721,9722,9723,9724,9725,9726,9727,9728,9729,9730,9731,9732,9733,9734,9735,9743,9744,9745,9746,9747,9748,9749,9750,9751,9752,9753,9754,9755,9756,9757,9758,9760,9761,9762,9763,9764,9765,9766,9779,9780,9781,9782,9783,9785,9788,9791,9792,9793,9794,9795,9796,9797,9807,9814,9834,9835,9836,9837,9838,9840,9841,9844,9847,9848,9849,9850,9851,9853,9854,9855,9858,9859,9860,9861,9862,9863,9867,9868,9869,9870,9871,9892,9893,9894,9896,9897,9911,9913,9914,9915,9917,9918,9919,9920,9921,9922,9923,9924,9925,9926,9929,9931,9932,9933,9934,9935,9940,9941,9942,9943,9944,9945,9946,9947,9948,9949,9950,9951,9953,9954,9955,9956,9959,9960,9961,9962,9968,9970,9971,9972,9973,9974,9975,9977,9978,9979,9980,9981,9982,9983,9984,9985,9986,9987,9988,9989,9991,10009,10010,10015,10030,10040,10048,10054,10055,10056,10057,10058,10059,10060,10061,10062,10063,10064,10065,10066,10067,10068,10084,10085,10088,10089,10090,10091,10092,10093,10095,10096,10097,10098};

  for(int i=0; i<368; i++){
      if(good_runs[i]==int(spill->hdr.run)){return true;}
  }
  return false;
});


