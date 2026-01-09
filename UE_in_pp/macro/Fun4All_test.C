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

#include <g4mbd/MbdDigitization.h>
#include <mbd/MbdReco.h>
#include <globalvertex/GlobalVertexReco.h>

#include <stdio.h>
#include <string.h>

#include "Calo_Calib.C"

#include <globalvertex/GlobalVertex.h>

#include "HIJetReco_nosub.C"
#include <InclusiveJet.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libg4mbd.so)
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libInclusiveJet.so)

#endif
void Fun4All_test(int nEvents = 100)
{
  
  int runnumber = 53517;
  int segment = 0;
  
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  CDBInterface::instance()->Verbosity(1);

  std::string outfilename = "/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/output_" + to_string(runnumber) + "_" + to_string(segment) + ".root";

  MbdReco *mbdreco = new MbdReco();
  //mbdreco->Verbosity(5);
  se->registerSubsystem(mbdreco);

  GlobalVertexReco *gvertex = new GlobalVertexReco();
  //gvertex->Verbosity(5);
  se->registerSubsystem(gvertex);
  
  Process_Calo_Calib();

  HIJetReco();

  InclusiveJet *myJetVal = new InclusiveJet("AntiKt_Tower_r04", "AntiKt_Truth_r04", outfilename.c_str());
  myJetVal->doJetTriggerCut(false);
  myJetVal->setLeadPtCut(5.0);
  myJetVal->doJetLeadPtCut(false);
  myJetVal->setPtRange(2, 100);
  myJetVal->setEtaRange(-0.7, 0.7);
  myJetVal->doTruth(0);
  myJetVal->doSeeds(0);
  myJetVal->doTowers(0);
  myJetVal->doTruthParticles(0);
  myJetVal->doTracks(0);
  myJetVal->doTopoclusters(0);
  myJetVal->doEmcalClusters(0);
  se->registerSubsystem(myJetVal);
  
  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DSTcalo");
  in2->AddFile("/direct/sphenix+lustre01/sphnxpro/production2/run2pp/physics/caloyfitting/ana509_2024p022_v001A/run_00053500_00053600/DST_CALOFITTING_run2pp_ana509_2024p022_v001-00053517-00000.root");
  se->registerInputManager(in2);

  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("DSTtracks");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00000.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00001.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00002.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00003.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00004.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00005.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00006.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00007.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00008.root");
  in3->AddFile("/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053500_00053600/dst/DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001-00053517-00009.root");
  se->registerInputManager(in3);

  //Fun4AllInputManager *intrue2 = new Fun4AllRunNodeInputManager("DST_GEO");
  //std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
  //intrue2->AddFile(geoLocation);
  //se->registerInputManager(intrue2);

  se->run(nEvents);
  se->End();
  CDBInterface::instance()->Print();  // print used DB files
  gSystem->Exit(0);
  return;

}
