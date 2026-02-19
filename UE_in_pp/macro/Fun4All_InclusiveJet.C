#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <g4centrality/PHG4CentralityReco.h>
#include <caloreco/RawClusterBuilderTopo.h>

//#include "HIJetReco.C"
//#include "HIJetRecoOrig.C"
#include "HIJetReco_nosub.C"
#include <InclusiveJet.h>

#include <g4mbd/MbdDigitization.h>
#include <mbd/MbdReco.h>
#include <globalvertex/GlobalVertexReco.h>

#include <stdio.h>
#include <string.h>

#include <Calo_Calib.C>
//#include "Calo_Calib_120125.C"

#include "mdctreemaker/MDCTreeMaker.h"

#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <jetbase/ClusterJetInput.h>

#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/DetermineTowerRho.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <jetbackground/TowerRho.h>

#include <globalvertex/GlobalVertex.h>

#include "/sphenix/user/egm2153/calib_study/sPHENIX-software/sPHENIX-macros/common/Sys_Calo_Local.C"
//#include <Sys_Calo.C>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libInclusiveJet.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libg4mbd.so)
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libmdctreemaker.so)

#endif
void Fun4All_InclusiveJet(int nEvents = 100, int seg = 0, int isSim = 0, const char *filelistcalo = "dst_calo_cluster.list", const char *type = "run22_jet10", const char *filelistjet = "dst_jet.list")
{
  
  Fun4AllServer *se = Fun4AllServer::instance();
  int verbosity = 0;

  std::ifstream calolist;
  std::string caloline;
  std::ifstream jetlist;
  std::string jetline;
  std::string mcfilelistcalo;
  std::string filelisttruth;
  std::string filelistglobal;
  std::string filelisttruthparticle;
  std::string filelisttracks;

  int runnumber = 0;
  int segment = 0;

  if (!isSim) {
    calolist.open(filelistcalo);
    for (int i = 0; i < seg + 1; i++) 
    {
      getline(calolist, caloline);
    }
    // uncomment to add additional data files: start
    //jetlist.open(filelistjet);
    //for (int i = 0; i < seg + 1; i++) 
    //{
    //  getline(jetlist, jetline);
    //}
    // uncomment to add additional data files:end
    std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(caloline);
    runnumber = runseg.first;
    segment = runseg.second;
    // uncomment to add additional data files: start
    //std::pair<int, int> jetrunseg = Fun4AllUtils::GetRunSegment(jetline);
    //if ((runseg.first != jetrunseg.first) || (runseg.second != jetrunseg.second)) {
    //  std::cout << "input files don't match. exiting now" << std::endl;
    //  return;
    //}
    // uncomment to add additional data files:end
  } else {
    std::ostringstream sseg;
    sseg << std::setw(4) << std::setfill('0') << seg;
    std::cout << sseg.str() << std::endl;
    if (!strcmp(type,"run22_jet10")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet10_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet10_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet10_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet10_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run22_jet20")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet20_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet20_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet20_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet20_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run22_jet30")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet30_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet30_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet30_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet30_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run22_jet50")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet50_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet50_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet50_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_jet50_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run22_mb")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_mb_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_mb_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_mb_dst_global_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_mb_dst_truthinfo_" + sseg.str() + ".list";
      filelisttracks = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_mb_dst_tracks_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run22_herwig_jet10")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_herwig_jet10_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_herwig_jet10_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_herwig_jet10_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_herwig_jet10_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run22_herwig_jet30")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_herwig_jet30_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_herwig_jet30_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_herwig_jet30_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run22_herwig_jet30_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run28_jet10")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet10_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet10_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet10_dst_mbd_epd_" + sseg.str() + ".list"; 
      //filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet10_dst_global_" + sseg.str() + ".list"; // used for tracking info
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet10_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run28_jet15")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet15_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet15_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet15_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet15_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run28_jet20")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet20_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet20_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet20_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet20_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run28_jet30")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet30_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet30_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet30_dst_mbd_epd_" + sseg.str() + ".list"; 
      //filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet30_dst_global_" + sseg.str() + ".list"; // used for tracking info
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet30_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run28_jet50")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet50_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet50_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet50_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet50_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run28_jet70")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet70_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet70_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet70_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run28_jet70_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run21_herwig_jet10")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_herwig_jet10_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_herwig_jet10_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_herwig_jet10_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_herwig_jet10_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run21_herwig_jet30")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_herwig_jet30_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_herwig_jet30_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_herwig_jet30_dst_mbd_epd_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_herwig_jet30_g4hits_" + sseg.str() + ".list";
    } else if (!strcmp(type,"run21_mb")) {
      mcfilelistcalo = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_mb_dst_calo_cluster_" + sseg.str() + ".list";
      filelisttruth = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_mb_dst_truth_jet_" + sseg.str() + ".list";
      filelistglobal = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_mb_dst_global_" + sseg.str() + ".list"; 
      filelisttruthparticle = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mc_list_files/run21_mb_dst_truthinfo_" + sseg.str() + ".list";
    } else {
      std::cout << "unknown simulation dst type" << std::endl;
    }
  }

  std::string outfilename;
  if (!isSim) {
    if (!strcmp(type,"1.5mrad_calofitting")) {
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/output_1.5mrad_ana468_calofitting_jettrigeff_" + to_string(runnumber) + "_" + to_string(segment) + ".root";
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/output_1.5mrad_ana468_calofitting_mb_events_" + to_string(runnumber) + "_" + to_string(segment) + ".root";
    } else if (!strcmp(type,"1.5mrad")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/output1.5mrad_ana468_" + to_string(runnumber) + "_" + to_string(segment) + ".root";
      } else if (!strcmp(type,"0mrad")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/output_0mrad_ana468_" + to_string(runnumber) + "_" + to_string(segment) + ".root";
    } else if (!strcmp(type,"jetcalo")) {
      //outfilename = "test_output_ana509_" + to_string(runnumber) + "_" + to_string(segment) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/output_ana509_v2_" + to_string(runnumber) + "_" + to_string(segment) + ".root"; // edited
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/output_ana509_jet_trig_" + to_string(runnumber) + "_" + to_string(segment) + ".root";
    } else if (!strcmp(type, "mbd_eff")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/output_mbd_eff_" + to_string(runnumber) + "_" + to_string(segment) + ".root";
    } else {
      std::cout << "unknown data dst type" << std::endl;
    }
  } else {
    if (!strcmp(type,"run22_jet10")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_jet10_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_jet10_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type, "run22_jet20")) {
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_jet20_3sigma_output_" + to_string(seg) + ".root";
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_jet20_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    }  else if (!strcmp(type,"run22_jet30")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_jet30_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_jet30_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type, "run22_jet50")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_jet50_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_jet50_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run22_mb")) { 
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_mb_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_mb_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run22_herwig_jet10")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_herwig_jet10_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_herwig_jet10_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run22_herwig_jet30")) {
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_herwig_jet30_3sigma_output_" + to_string(seg) + ".root";
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_herwig_jet30_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run28_jet10")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet10_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet10_total_events_zvtx_lt_60_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run28_jet15")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet15_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet15_total_events_zvtx_lt_60_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type, "run28_jet20")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet20_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet20_total_events_zvtx_lt_60_output_" + to_string(seg) + ".root";
    }  else if (!strcmp(type,"run28_jet30")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet30_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet30_total_events_zvtx_lt_60_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type, "run28_jet50")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet50_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet50_total_events_zvtx_lt_60_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type, "run28_jet70")) {
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet70_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet70_total_events_zvtx_lt_60_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run21_mb")) { 
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_mb_3sigma_output_" + to_string(seg) + ".root";
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_mb_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run21_herwig_jet10")) {
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_herwig_jet10_3sigma_output_" + to_string(seg) + ".root";
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_herwig_jet10_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
      //outfilename = "test_sim_run21_herwig_jet10_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run21_herwig_jet30")) {
      //outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_herwig_jet30_3sigma_output_" + to_string(seg) + ".root";
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_herwig_jet30_total_events_zvtx_lt_30_output_" + to_string(seg) + ".root";
    } else if (!strcmp(type,"run21_mb")) { 
      outfilename = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_mb_output_" + to_string(seg) + ".root";
    }
     else {
      std::cout << "unknown simulation dst type" << std::endl;
    }
  }
  
  se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();
  if (!isSim) {
    rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
    rc->set_uint64Flag("TIMESTAMP", runnumber);
  } else {
    rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
    rc->set_uint64Flag("TIMESTAMP", 28);
  }
    CDBInterface::instance()->Verbosity(1);
    
  MbdReco *mbdreco = new MbdReco();
  se->registerSubsystem(mbdreco);

  GlobalVertexReco *gvertex = new GlobalVertexReco();
  //gvertex->Verbosity(2);
  se->registerSubsystem(gvertex);

  Process_Calo_Calib();
  se->Print();

  //Register_Tower_sys();
  
  RawClusterBuilderTopo* ClusterBuilder = new RawClusterBuilderTopo("HcalRawClusterBuilderTopo");
  ClusterBuilder->Verbosity(verbosity);
  //ClusterBuilder->setInputTowerNodePrefix("TOWERINFO_CALIB_SYST1");
  ClusterBuilder->set_nodename("TOPOCLUSTER_ALLCALO");
  ClusterBuilder->set_enable_HCal(true);
  ClusterBuilder->set_enable_EMCal(true);
  ClusterBuilder->set_noise(0.0053, 0.0351, 0.0684); // 3sigma of pedestal noise
  ClusterBuilder->set_significance(4.0, 2.0, 1.0);
  ClusterBuilder->allow_corner_neighbor(true);
  ClusterBuilder->set_do_split(true);
  ClusterBuilder->set_minE_local_max(1.0, 2.0, 0.5);
  ClusterBuilder->set_R_shower(0.025);
  ClusterBuilder->set_use_only_good_towers(true);
  ClusterBuilder->set_absE(true);
  se->registerSubsystem(ClusterBuilder);

  RawClusterBuilderTopo* ClusterBuilder2 = new RawClusterBuilderTopo("HcalRawClusterBuilderTopo2");
  ClusterBuilder2->Verbosity(verbosity);
  ClusterBuilder2->set_nodename("TOPOCLUSTER_ALLCALO_2SIGMA");
  ClusterBuilder2->set_enable_HCal(true);
  ClusterBuilder2->set_enable_EMCal(true);
  ClusterBuilder2->set_noise(0.0035, 0.0234, 0.0456); // 2sigma of pedestal noise
  ClusterBuilder2->set_significance(4.0, 2.0, 1.0);
  ClusterBuilder2->allow_corner_neighbor(true);
  ClusterBuilder2->set_do_split(true);
  ClusterBuilder2->set_minE_local_max(1.0, 2.0, 0.5);
  ClusterBuilder2->set_R_shower(0.025);
  ClusterBuilder2->set_use_only_good_towers(true);
  ClusterBuilder2->set_absE(true);
  se->registerSubsystem(ClusterBuilder2);
  

  Enable::VERBOSITY = verbosity;
  HIJetReco();

  JetCalib *jetCalib04 = new JetCalib("JetCalib04");
  jetCalib04->set_InputNode("AntiKt_Tower_r04");
  jetCalib04->set_OutputNode("AntiKt_Tower_r04_calib");
  jetCalib04->set_JetRadius(0.4);
  jetCalib04->set_ZvrtxNode("GlobalVertexMap");
  jetCalib04->set_ApplyZvrtxDependentCalib(true);
  jetCalib04->set_ApplyEtaDependentCalib(true);
  se->registerSubsystem(jetCalib04);
    
  InclusiveJet *myJetVal = new InclusiveJet("AntiKt_Tower_r04", "AntiKt_Truth_r04", outfilename.c_str());
  myJetVal->doJetTriggerCut(false);
  myJetVal->setLeadPtCut(7.0);
  myJetVal->doJetLeadPtCut(true);
  myJetVal->setPtRange(2, 100);
  myJetVal->setEtaRange(-2.5, 2.5);
  if (isSim) myJetVal->doTruth(1);
  myJetVal->doSeeds(0);
  myJetVal->doTowers(0);
  if (isSim) myJetVal->doTruthParticles(0);
  myJetVal->doTracks(0);
  myJetVal->doTopoclusters(0);
  myJetVal->doEmcalClusters(0);
  myJetVal->doMBDeff(false);
  se->registerSubsystem(myJetVal);
  
  //MDCTreeMaker *tt = new MDCTreeMaker("MDCTreeMaker", outfilename, isSim, 1, 0);
  //tt->set_useMBD(true); 
  //tt->set_useEMCal(false);
  //tt->set_useHCal(false);
  //se->registerSubsystem( tt ); 

  if (!isSim) {
    // uncomment to add additional data files:start
    //Fun4AllInputManager *in3 = new Fun4AllDstInputManager("DSTjet");
    //in3->AddListFile(jetline,1);
    //se->registerInputManager(in3);
    // uncomment to add additional data files:end

    Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTcalo");
    in2->AddListFile(caloline);
    se->registerInputManager(in2);

    Fun4AllInputManager *intrue2 = new Fun4AllRunNodeInputManager("DST_GEO");
    std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
    intrue2->AddFile(geoLocation);
    se->registerInputManager(intrue2);

  } else {
    Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTcalo");
    in2->AddListFile(mcfilelistcalo);
    se->registerInputManager(in2);

    Fun4AllInputManager *intrue = new Fun4AllDstInputManager("DSTtruth");
    intrue->AddListFile(filelisttruth);
    if (strcmp(type, "run21_herwig_jet10") && strcmp(type, "run21_herwig_jet30")) se->registerInputManager(intrue);

    Fun4AllInputManager *in3 = new Fun4AllDstInputManager("DSTglobal");
    in3->AddListFile(filelistglobal);
    se->registerInputManager(in3);

    Fun4AllInputManager *in4 = new Fun4AllDstInputManager("DSTtruthparticle");
    in4->AddListFile(filelisttruthparticle);
    se->registerInputManager(in4);

    if(!strcmp(type,"run22_mb")) {
      Fun4AllInputManager *in5 = new Fun4AllDstInputManager("DSTtracks");
      in5->AddListFile(filelisttracks);
      se->registerInputManager(in5);
    }

    Fun4AllInputManager *intrue2 = new Fun4AllRunNodeInputManager("DST_GEO");
    std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
    intrue2->AddFile(geoLocation);
    se->registerInputManager(intrue2);

  }
  
  se->run(nEvents);
  se->End();
  CDBInterface::instance()->Print();  // print used DB files
  gSystem->Exit(0);
  return;

}
