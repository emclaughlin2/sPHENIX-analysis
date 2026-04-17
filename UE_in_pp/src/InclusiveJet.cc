//module for producing a TTree with jet information for doing jet validation studies
// for questions/bugs please contact Virginia Bailey vbailey13@gsu.edu
#include <fun4all/Fun4AllBase.h>
#include <InclusiveJet.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <jetbase/JetMap.h>
#include <jetbase/JetContainer.h>
#include <jetbase/Jetv2.h>
#include <jetbase/Jetv1.h>
#include <centrality/CentralityInfo.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <ffarawobjects/Gl1Packet.h>
#include <jetbackground/TowerBackground.h>
#include <calobase/RawTowerDefs.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <ffaobjects/EventHeaderv1.h>

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdOut.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>

// G4 truth includes
#include <g4eval/SvtxEvalStack.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <ffaobjects/RunHeader.h>

#include <TTree.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include <TTree.h>
#include <TVector3.h>

//____________________________________________________________________________..
InclusiveJet::InclusiveJet(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename):
  SubsysReco("InclusiveJet_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-0.7, 0.7)
  , m_ptRange(5, 100)
  , m_doTruthJets(0)
  , m_doTruth(0)
  , m_doTowers(0)
  , m_doEmcalClusters(0)
  , m_doTopoclusters(0)
  , m_doSeeds(0)
  , m_doTracks(0)
  , m_T(nullptr)
  , m_event(-1)
  , m_nTruthJet(-1)
  , m_nJet(-1)
  , m_triggerVector()
  , m_liveTriggerVector()
  , m_nComponent()
  , m_eta()
  , m_phi()
  , m_e()
  , m_pt()
  , m_calibpt()
  , m_jetEmcalE()
  , m_jetIhcalE()
  , m_jetOhcalE()
  , m_truthNComponent()
  , m_truthEta()
  , m_truthPhi()
  , m_truthE()
  , m_truthPt()
  , m_eta_rawseed()
  , m_phi_rawseed()
  , m_pt_rawseed()
  , m_e_rawseed()
  , m_rawseed_cut()
  , m_eta_subseed()
  , m_phi_subseed()
  , m_pt_subseed()
  , m_e_subseed()
  , m_subseed_cut()
  , runnumber(-9999)
{
  std::cout << "InclusiveJet::InclusiveJet(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
InclusiveJet::~InclusiveJet()
{
  std::cout << "InclusiveJet::~InclusiveJet() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
//____________________________________________________________________________..
int InclusiveJet::Init(PHCompositeNode *topNode)
{
  std::cout << "InclusiveJet::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  std::cout << "InclusiveJet::Init - Output to " << m_outputFileName << std::endl;

  // configure Tree
  m_T = new TTree("T", "MyJetAnalysis Tree");
  m_T->Branch("m_event", &m_event, "event/I");
  m_T->Branch("nJet", &m_nJet, "nJet/I");
  m_T->Branch("zvtx", &m_zvtx);
  m_T->Branch("mbd_t0",&m_mbd_t0,"mbd_t0/F");
  m_T->Branch("mbd_ts",&m_mbd_ts,"mbd_ts/F");
  m_T->Branch("mbd_tn",&m_mbd_tn,"mbd_tn/F");
  if (m_doMBDeff && m_doTruth) {
    m_T->Branch("zsvtx", &m_zsvtx);
    m_T->Branch("zsiliconvtx",&m_zsiliconvtx);
    m_T->Branch("ztruthvtx", &m_ztruthvtx);
    m_T->Branch("svtxVector",&m_svtxVector);
    m_T->Branch("svtxBcoVector",&m_svtxBcoVector);
    m_T->Branch("siliconVector",&m_siliconVector);
    m_T->Branch("siliconBcoVector",&m_siliconBcoVector);
  }
  else if (m_doMBDeffsyst && m_doTruth) {
    m_T->Branch("ztruthvtx", &m_ztruthvtx);
  }
  m_T->Branch("nComponent", &m_nComponent);
  m_T->Branch("triggerVector", &m_triggerVector);
  m_T->Branch("liveTriggerVector", &m_liveTriggerVector);

  m_T->Branch("eta", &m_eta);
  m_T->Branch("phi", &m_phi);
  m_T->Branch("e", &m_e);
  m_T->Branch("pt", &m_pt);
  if (m_doCalibJet) {
    m_T->Branch("calibpt", &m_calibpt);
  }

  m_T->Branch("jetEmcalE", &m_jetEmcalE);
  m_T->Branch("jetIhcalE", &m_jetIhcalE);
  m_T->Branch("jetOhcalE", &m_jetOhcalE);
  m_T->Branch("jettime", &m_jettime);
  m_T->Branch("jetunweighttime", &m_jetunweighttime);

  if(m_doTruthJets){
    m_T->Branch("nTruthJet", &m_nTruthJet);
    m_T->Branch("truthNComponent", &m_truthNComponent);
    m_T->Branch("truthEta", &m_truthEta);
    m_T->Branch("truthPhi", &m_truthPhi);
    m_T->Branch("truthE", &m_truthE);
    m_T->Branch("truthPt", &m_truthPt);
  }

  if(m_doSeeds){
    m_T->Branch("rawseedEta", &m_eta_rawseed);
    m_T->Branch("rawseedPhi", &m_phi_rawseed);
    m_T->Branch("rawseedPt", &m_pt_rawseed);
    m_T->Branch("rawseedE", &m_e_rawseed);
    m_T->Branch("rawseedCut", &m_rawseed_cut);
    m_T->Branch("subseedEta", &m_eta_subseed);
    m_T->Branch("subseedPhi", &m_phi_subseed);
    m_T->Branch("subseedPt", &m_pt_subseed);
    m_T->Branch("subseedE", &m_e_subseed);
    m_T->Branch("subseedCut", &m_subseed_cut);
  }

  if(m_doTowers) {
    m_T->Branch("emcaln",&m_emcaln,"emcaln/I");
    m_T->Branch("emcale",m_emcale,"emcale[emcaln]/F");
    m_T->Branch("emcalchi2", m_emcalchi2,"emcalchi2[emcaln]/F");
    m_T->Branch("emcaleta",m_emcaleta,"emcaleta[emcaln]/F");
    m_T->Branch("emcalphi",m_emcalphi,"emcalphi[emcaln]/F");
    m_T->Branch("emcalieta",m_emcalieta,"emcalieta[emcaln]/I");
    m_T->Branch("emcaliphi",m_emcaliphi,"emcaliphi[emcaln]/I");
    m_T->Branch("emcalstatus",m_emcalstatus,"emcalstatus[emcaln]/I");
    //m_T->Branch("emcaltime",m_emcaltime,"emcaltime[emcaln]/F");

    m_T->Branch("ihcaln",&m_ihcaln,"ihcaln/I");
    m_T->Branch("ihcale",m_ihcale,"ihcale[ihcaln]/F");
    m_T->Branch("ihcalchi2", m_ihcalchi2,"ihcalchi2[ihcaln]/F");
    m_T->Branch("ihcaleta",m_ihcaleta,"ihcaleta[ihcaln]/F");
    m_T->Branch("ihcalphi",m_ihcalphi,"ihcalphi[ihcaln]/F");
    m_T->Branch("ihcalieta",m_ihcalieta,"ihcalieta[ihcaln]/I");
    m_T->Branch("ihcaliphi",m_ihcaliphi,"ihcaliphi[ihcaln]/I");
    m_T->Branch("ihcalstatus",m_ihcalstatus,"ihcalstatus[ihcaln]/I");
    //m_T->Branch("ihcaltime",m_ihcaltime,"ihcaltime[ihcaln]/F");

    m_T->Branch("ohcaln",&m_ohcaln,"ohcaln/I");
    m_T->Branch("ohcale",m_ohcale,"ohcale[ohcaln]/F");
    m_T->Branch("ohcalchi2", m_ohcalchi2,"ohcalchi2[ohcaln]/F");
    m_T->Branch("ohcaleta",m_ohcaleta,"ohcaleta[ohcaln]/F");
    m_T->Branch("ohcalphi",m_ohcalphi,"ohcalphi[ohcaln]/F");
    m_T->Branch("ohcalieta",m_ohcalieta,"ohcalieta[ohcaln]/I");
    m_T->Branch("ohcaliphi",m_ohcaliphi,"ohcaliphi[ohcaln]/I");
    m_T->Branch("ohcalstatus",m_ohcalstatus,"ohcalstatus[ohcaln]/I");
    //m_T->Branch("ohcaltime",m_ohcaltime,"ohcaltime[ohcaln]/F");
  }

  if(m_doEmcalClusters) {
    m_T->Branch("emcal_clsmult",&m_emcal_clsmult,"emcal_clsmult/I");
    m_T->Branch("emcal_cluster_e",m_emcal_cluster_e,"emcal_cluster_e[emcal_clsmult]/F");
    m_T->Branch("emcal_cluster_eta",m_emcal_cluster_eta,"emcal_cluster_eta[emcal_clsmult]/F");
    m_T->Branch("emcal_cluster_phi",m_emcal_cluster_phi,"emcal_cluster_phi[emcal_clsmult]/F");
  }

  if(m_doTopoclusters) {
    m_T->Branch("clsmult",&m_clsmult,"clsmult/I");
    m_T->Branch("cluster_e",m_cluster_e,"cluster_e[clsmult]/F");
    m_T->Branch("cluster_eta",m_cluster_eta,"cluster_eta[clsmult]/F");
    m_T->Branch("cluster_phi",m_cluster_phi,"cluster_phi[clsmult]/F");
    m_T->Branch("cluster_ntowers",m_cluster_ntowers,"cluster_ntowers[clsmult]/I");
    m_T->Branch("cluster_tower_e",m_cluster_tower_e,"cluster_tower_e[clsmult][500]/F");
    m_T->Branch("cluster_tower_calo",m_cluster_tower_calo,"cluster_tower_calo[clsmult][500]/I");
    m_T->Branch("cluster_tower_ieta",m_cluster_tower_ieta,"cluster_tower_ieta[clsmult][500]/I");
    m_T->Branch("cluster_tower_iphi",m_cluster_tower_iphi,"cluster_tower_iphi[clsmult][500]/I");
    
    m_T->Branch("clsmult2",&m_clsmult2,"clsmult2/I");
    m_T->Branch("cluster2_e",m_cluster2_e,"cluster2_e[clsmult2]/F");
    m_T->Branch("cluster2_eta",m_cluster2_eta,"cluster2_eta[clsmult2]/F");
    m_T->Branch("cluster2_phi",m_cluster2_phi,"cluster2_phi[clsmult2]/F");
    m_T->Branch("cluster2_ntowers",m_cluster2_ntowers,"cluster2_ntowers[clsmult2]/I");
    m_T->Branch("cluster2_tower_e",m_cluster2_tower_e,"cluster2_tower_e[clsmult2][500]/F");
    m_T->Branch("cluster2_tower_calo",m_cluster2_tower_calo,"cluster2_tower_calo[clsmult2][500]/I");
    m_T->Branch("cluster2_tower_ieta",m_cluster2_tower_ieta,"cluster2_tower_ieta[clsmult2][500]/I");
    m_T->Branch("cluster2_tower_iphi",m_cluster2_tower_iphi,"cluster2_tower_iphi[clsmult2][500]/I");

    m_T->Branch("clsmult4",&m_clsmult4,"clsmult4/I");
    m_T->Branch("cluster4_e",m_cluster4_e,"cluster4_e[clsmult4]/F");
    m_T->Branch("cluster4_eta",m_cluster4_eta,"cluster4_eta[clsmult4]/F");
    m_T->Branch("cluster4_phi",m_cluster4_phi,"cluster4_phi[clsmult4]/F");
    m_T->Branch("cluster4_ntowers",m_cluster4_ntowers,"cluster4_ntowers[clsmult4]/I");
    m_T->Branch("cluster4_tower_e",m_cluster4_tower_e,"cluster4_tower_e[clsmult4][500]/F");
    m_T->Branch("cluster4_tower_calo",m_cluster4_tower_calo,"cluster4_tower_calo[clsmult4][500]/I");
    m_T->Branch("cluster4_tower_ieta",m_cluster4_tower_ieta,"cluster4_tower_ieta[clsmult4][500]/I");
    m_T->Branch("cluster4_tower_iphi",m_cluster4_tower_iphi,"cluster4_tower_iphi[clsmult4][500]/I");
  }

  if(m_doTruth) {
    m_T->Branch("truthpar_n",&truthpar_n,"truthpar_n/I");
    m_T->Branch("truthpar_pz",truthpar_pz,"truthpar_pz[truthpar_n]/F");
    m_T->Branch("truthpar_pt",truthpar_pt,"truthpar_pt[truthpar_n]/F");
    m_T->Branch("truthpar_e",truthpar_e,"truthpar_e[truthpar_n]/F");
    m_T->Branch("truthpar_eta",truthpar_eta,"truthpar_eta[truthpar_n]/F");
    m_T->Branch("truthpar_phi",truthpar_phi,"truthpar_phi[truthpar_n]/F");
    m_T->Branch("truthpar_pid",truthpar_pid,"truthpar_pid[truthpar_n]/I");
    m_T->Branch("old_truthpar_n",&old_truthpar_n,"old_truthpar_n/I");
    m_T->Branch("old_truthpar_pz",old_truthpar_pz,"old_truthpar_pz[old_truthpar_n]/F");
    m_T->Branch("old_truthpar_pt",old_truthpar_pt,"old_truthpar_pt[old_truthpar_n]/F");
    m_T->Branch("old_truthpar_e",old_truthpar_e,"old_truthpar_e[old_truthpar_n]/F");
    m_T->Branch("old_truthpar_eta",old_truthpar_eta,"old_truthpar_eta[old_truthpar_n]/F");
    m_T->Branch("old_truthpar_phi",old_truthpar_phi,"old_truthpar_phi[old_truthpar_n]/F");
    m_T->Branch("old_truthpar_pid",old_truthpar_pid,"old_truthpar_pid[old_truthpar_n]/I");
  }

  if(m_doTracks) {
    m_T->Branch("trkmult",&m_trkmult,"trkmult/I");
    m_T->Branch("tr_p",m_tr_p,"tr_p[trkmult]/F");
    m_T->Branch("tr_pt",m_tr_pt,"tr_pt[trkmult]/F");
    m_T->Branch("tr_eta",m_tr_eta,"tr_eta[trkmult]/F");
    m_T->Branch("tr_phi",m_tr_phi,"tr_phi[trkmult]/F");
    m_T->Branch("tr_charge",m_tr_charge,"tr_charge[trkmult]/I");
    m_T->Branch("tr_chisq",m_tr_chisq,"tr_chisq[trkmult]/F");
    m_T->Branch("tr_ndf",m_tr_ndf,"tr_ndf[trkmult]/I");
    m_T->Branch("tr_nintt",m_tr_nintt,"tr_nintt[trkmult]/I");
    m_T->Branch("tr_nmaps",m_tr_nmaps,"tr_nmaps[trkmult]/I");
    m_T->Branch("tr_ntpc",m_tr_ntpc,"tr_ntpc[trkmult]/I");
    m_T->Branch("tr_quality",m_tr_quality,"tr_quality[trkmult]/F");
    m_T->Branch("tr_crossing",m_tr_crossing,"tr_crossing[trkmult]/I");
    m_T->Branch("tr_cemc_eta",m_tr_cemc_eta,"tr_cemc_eta[trkmult]/F");
    m_T->Branch("tr_cemc_phi",m_tr_cemc_phi,"tr_cemc_phi[trkmult]/F");
    m_T->Branch("tr_ihcal_eta",m_tr_ihcal_eta,"tr_ihcal_eta[trkmult]/F");
    m_T->Branch("tr_ihcal_phi",m_tr_ihcal_phi,"tr_ihcal_phi[trkmult]/F");
    m_T->Branch("tr_ohcal_eta",m_tr_ohcal_eta,"tr_ohcal_eta[trkmult]/F");
    m_T->Branch("tr_ohcal_phi",m_tr_ohcal_phi,"tr_ohcal_phi[trkmult]/F");
    m_T->Branch("tr_outer_cemc_eta",m_tr_outer_cemc_eta,"tr_outer_cemc_eta[trkmult]/F");
    m_T->Branch("tr_outer_cemc_phi",m_tr_outer_cemc_phi,"tr_outer_cemc_phi[trkmult]/F");
    m_T->Branch("tr_outer_ihcal_eta",m_tr_outer_ihcal_eta,"tr_outer_ihcal_eta[trkmult]/F");
    m_T->Branch("tr_outer_ihcal_phi",m_tr_outer_ihcal_phi,"tr_outer_ihcal_phi[trkmult]/F");
    m_T->Branch("tr_outer_ohcal_eta",m_tr_outer_ohcal_eta,"tr_outer_ohcal_eta[trkmult]/F");
    m_T->Branch("tr_outer_ohcal_phi",m_tr_outer_ohcal_phi,"tr_outer_ohcal_phi[trkmult]/F");
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::InitRun(PHCompositeNode *topNode)
{
  std::cout << "InclusiveJet::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  
  runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!runheader) {
    std::cout << "can't find runheader" << std::endl;
    return 1;
  }
  runnumber = runheader->get_RunNumber();
  std::ifstream file = std::ifstream("/sphenix/user/samfred/projects/mbdt0/histmaking/MbdOut.corr");
  std::string line;
  while (getline(file,line)) {
    int irunnum;
    float t0;
    std::istringstream iss(line);
    iss >> irunnum >> t0;
    if (irunnum == runnumber) {m_t0corr = t0; break; }
  }


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::process_event(PHCompositeNode *topNode)
{
  //  std::cout << "InclusiveJet::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  ++m_event;

  // interface to reco jets
  JetContainer* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  if (!jets)
    {
      std::cout
	<< "MyJetAnalysis::process_event - Error can not find DST Reco JetContainer node "
	<< m_recoJetName << std::endl;
      exit(-1);
    }

  JetContainer* calibjets = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_r04_calib");
  if (m_doCalibJet && !calibjets)
    {
      std::cout
  << "MyJetAnalysis::process_event - Error can not find DST Reco Calib JetContainer node: AntiKt_Tower_r04_calib"
  << std::endl;
      exit(-1);
    }

  //interface to truth jets
  //JetMap* jetsMC = findNode::getClass<JetMap>(topNode, m_truthJetName);
  JetContainer* jetsMC = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  if (!jetsMC && m_doTruthJets)
    {
      std::cout
	<< "MyJetAnalysis::process_event - Error can not find DST Truth JetMap node "
	<< m_truthJetName << std::endl;
      exit(-1);
    }
  
  // interface to jet seeds
  JetContainer* seedjetsraw = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
  if (!seedjetsraw && m_doSeeds)
    {
      std::cout
	<< "MyJetAnalysis::process_event - Error can not find DST raw seed jets "
<< std::endl;
      exit(-1);
    }

  JetContainer* seedjetssub = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsSub_r02");
  if (!seedjetssub && m_doSeeds)
    {
      std::cout
<< "MyJetAnalysis::process_event - Error can not find DST subtracted seed jets "
<< std::endl;
      exit(-1);
    }
  
  
  //zvertex
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
    {
      std::cout
        << "MyJetAnalysis::process_event - Error can not find global vertex  node "
        << std::endl;
      exit(-1);
    }
  if (vertexmap->empty())
    {
      std::cout
        << "MyJetAnalysis::process_event - global vertex node is empty "
        << std::endl;
        m_zvtx = -9999;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  else
    {
      for (auto giter = vertexmap->begin(); giter != vertexmap->end(); giter++) {
        GlobalVertex *globalVertex = giter->second;

        //std::cout << "Event: " << m_event << std::endl;
        //globalVertex->identify();

        auto mbdStartIter = globalVertex->find_vertexes(GlobalVertex::MBD);
        auto mbdEndIter = globalVertex->end_vertexes();

        for (auto iter = mbdStartIter; iter != mbdEndIter; ++iter)
        {
          const auto &[type, vertexVec] = *iter;
          if (type != GlobalVertex::MBD) continue;
          for (const auto *vertex : vertexVec)
          {
            if (!vertex) continue; 
            //std::cout << "mbd vertex beam crossing: " << vertex->get_beam_crossing() << std::endl;
            m_zvtx = vertex->get_z();
          }
        }
        if (m_doMBDeff) {
          auto svtxStartIter = globalVertex->find_vertexes(GlobalVertex::SVTX);
          auto svtxEndIter = globalVertex->end_vertexes();

          for (auto iter = svtxStartIter; iter != svtxEndIter; ++iter)
          {
            const auto &[type, vertexVec] = *iter;
            if (type != GlobalVertex::SVTX) continue;
            for (const auto *vertex : vertexVec)
            {
              if (!vertex) continue; 
              //std::cout << vertex->get_z() << " " << vertex->get_beam_crossing() << std::endl;
              m_svtxVector.push_back(vertex->get_z());
              m_svtxBcoVector.push_back(vertex->get_beam_crossing());
              //std::cout << "svtx vertex beam crossing: " << vertex->get_beam_crossing() << std::endl;
              if (vertex->get_beam_crossing() != 0) { continue; }
              m_zsvtx = vertex->get_z();
            }
          }
        }
        if ((m_doMBDeff || m_doMBDeffsyst) && m_doTruth) {
          auto truthStartIter = globalVertex->find_vertexes(GlobalVertex::TRUTH);
          auto truthEndIter = globalVertex->end_vertexes();

          for (auto iter = truthStartIter; iter != truthEndIter; ++iter)
          {
            const auto &[type, vertexVec] = *iter;
            if (type != GlobalVertex::TRUTH) continue;
            for (const auto *vertex : vertexVec)
            {
              if (!vertex) continue; 
              //std::cout << "truth vertex beam crossing: " << vertex->get_beam_crossing() << std::endl;
              //if (vertex->get_beam_crossing() != 0) { continue; }
              m_ztruthvtx = vertex->get_z();
            }
          }
        }
      }
    }

    if (m_doMBDeff) {
      SvtxVertexMap *siliconvertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SiliconVertexMap");
      if (!siliconvertexmap)
      {
        std::cout << PHWHERE<< "SiliconVertexMap node is missing, can't collect tracks"<< std::endl;
        return Fun4AllReturnCodes::EVENT_OK;
      }

      // note : in one F4A tracker event, the BCO coverage can be up to 500 BCOs, so there can be several vertices in one F4A tracker event, so we loop over them
      for (auto it = siliconvertexmap->begin(); it != siliconvertexmap->end(); ++it)
      {
          auto vertex = it->second;
          if (!vertex) continue;

          m_siliconVector.push_back(vertex->get_z());
          m_siliconBcoVector.push_back(vertex->get_beam_crossing());

          if (vertex->get_beam_crossing() == 0) {
              m_zsiliconvtx = vertex->get_z();
          }
      }
    }

    //if (m_doMBDeff) {
      //PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      //if (!truthinfo) std::cout << PHWHERE << "PHG4TruthInfoContainer node is missing, can't collect truth vertex"<< std::endl;
      //PHG4TruthInfoContainer::VtxRange vtxrange = truthinfo->GetVtxRange();
      //for (PHG4TruthInfoContainer::ConstVtxIterator iter = vtxrange.first; iter != vtxrange.second; ++iter) {
        //PHG4VtxPoint *vtx = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
        //if (vtx->get_id() != 0) { continue; }
        //std::cout << " truth vertex: x = " << vtx->get_x() << " y = " << vtx->get_y() << " z = " << vtx->get_z() << " id = " << vtx->get_id() << std::endl;
      //}
    //}

    if (m_doMBDeff) { std::cout << "using vertex " << m_zvtx << " " << m_zsvtx << " " << m_ztruthvtx << " " << m_zsiliconvtx << std::endl; }
    else { std::cout << "using vertex " << m_zvtx << std::endl; }

 /*
  MbdVertexMap *mvertexmap = findNode::getClass<MbdVertexMap>(topNode,"MbdVertexMap");
  if (!mvertexmap)
  {
    std::cout << "MbdVertexMap node is missing" << std::endl;
  }
  if (mvertexmap && !mvertexmap->empty())
  {
    MbdVertex *mvtx = mvertexmap->begin()->second;
    std::cout << "Event: " << m_event << std::endl;
    mvtx->identify();
    if (mvtx)
    {
	    std::cout << "mbd vertex: " << mvtx->get_z() << std::endl;
    }
    else
    {
	    std::cout << "no mbd vertex" << std::endl;
    }
  }
  */

  MbdOut * mbdout = static_cast<MbdOut*>(findNode::getClass<MbdOut>(topNode,"MbdOut"));
   if(mbdout){
     m_mbd_t0 = mbdout->get_t0() - m_t0corr;
     m_mbd_ts = mbdout->get_time(0) - m_t0corr; // south side
     m_mbd_tn = mbdout->get_time(1) - m_t0corr; // north side
   }
   /*
  MbdPmtContainer* mbdpmt = findNode::getClass<MbdPmtContainer>(topNode,"MbdPmtContainer");
  m_mbdavgt[0] = 0;
  m_mbdavgt[1] = 0;
  m_mbdhit[0] = 0;
  m_mbdhit[1] = 0;
  if(mbdpmt) {
    for(int i=0; i<128; ++i) {
      MbdPmtHit* pmt = mbdpmt->get_pmt(i);
      if(pmt) {
        if(pmt->get_q() > 0.4) {
          ++m_mbdhit[i/64];
          m_mbdavgt[i/64] += pmt->get_time();
        }
      }
      }
      if (m_mbdhit[0] != 0) m_mbdavgt[0]/=m_mbdhit[0];
      if (m_mbdhit[1] != 0) m_mbdavgt[1]/=m_mbdhit[1];
    }
    */

  if (fabs(m_zvtx) > 60 && !m_doMBDeff && !m_doMBDeffsyst) {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  if (m_doMBDeffsyst && m_doTruth && fabs(m_ztruthvtx) > 60) {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  //calorimeter towers
  TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if(!towersEM3 || !towersIH3 || !towersOH3){
    std::cout
      <<"MyJetAnalysis::process_event - Error can not find raw tower node "
      << std::endl;
    exit(-1);
  }

  TowerInfoContainer *EMtowers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainer *IHtowers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *OHtowers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  if(!tower_geom || !tower_geomOH){
    std::cout
      <<"MyJetAnalysis::process_event - Error can not find raw tower geometry "
      << std::endl;
    exit(-1);
  }

  RawClusterContainer *topoclusters = findNode::getClass<RawClusterContainer>(topNode,"TOPOCLUSTER_ALLCALO"); 
  if (m_doTopoclusters && !topoclusters) {
        std::cout
    <<"MyJetAnalysis::process_event - Error can not find topoclusters "
    << std::endl;
  exit(-1);
  }

  RawClusterContainer *topoclusters2 = findNode::getClass<RawClusterContainer>(topNode,"TOPOCLUSTER_ALLCALO_2SIGMA"); 
  if (m_doTopoclusters && !topoclusters2) {
        std::cout
    <<"MyJetAnalysis::process_event - Error can not find 2 sigma topoclusters "
    << std::endl;
  exit(-1);
  }

  RawClusterContainer *topoclusters4 = findNode::getClass<RawClusterContainer>(topNode,"TOPOCLUSTER_ALLCALO_4SIGMA"); 
  if (m_doTopoclusters && !topoclusters4) {
        std::cout
    <<"MyJetAnalysis::process_event - Error can not find 2 sigma topoclusters "
    << std::endl;
  exit(-1);
  }

  RawClusterContainer *emcalclusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTERINFO_CEMC");
  if (m_doEmcalClusters && !emcalclusters) {
      std::cout
    <<"MyJetAnalysis::process_event - Error can not find EMCal clusters "
    << std::endl;
  exit(-1);
  }

  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (m_doTruth && !truthinfo) {
    std::cout << PHWHERE << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"<< std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if(m_doTracks){
    RawTowerGeomContainer_Cylinderv1* cemcGeomContainer = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC");
    RawTowerGeomContainer_Cylinderv1* ihcalGeomContainer = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALIN");
    RawTowerGeomContainer_Cylinderv1* ohcalGeomContainer = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALOUT");

    if(!cemcGeomContainer){
      std::cout << "ERROR: TOWERGEOM_CEMC not found" << std::endl;
    }
    if(!ihcalGeomContainer){
      std::cout << "ERROR: TOWERGEOM_HCALIN not found" << std::endl;
    }
    if(!ohcalGeomContainer){
      std::cout << "ERROR: TOWERGEOM_HCALOUT not found" << std::endl;
    }

    m_cemcRadius = cemcGeomContainer->get_radius();
    m_ihcalRadius = ihcalGeomContainer->get_radius();
    m_ohcalRadius = ohcalGeomContainer->get_radius();
    m_cemcOuterRadius = cemcGeomContainer->get_radius() + cemcGeomContainer->get_thickness();
    m_ihcalOuterRadius = ihcalGeomContainer->get_radius() + ihcalGeomContainer->get_thickness();
    m_ohcalOuterRadius = ohcalGeomContainer->get_radius() + ohcalGeomContainer->get_thickness();
  }

  //get reco jets
  m_nJet = 0;
  float leadpt = 0;
  if (m_doCalibJet && jets->size() != calibjets->size())
  {
    std::cout << "Jet containers have different sizes!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  for (unsigned int i = 0; i < jets->size(); ++i)
  {
      Jet* jet = jets->get_jet(i);
      Jet* calibjet = calibjets->get_jet(i);

      if (!jet || !calibjet) continue;
      bool eta_cut = (jet->get_eta() >= m_etaRange.first) and (jet->get_eta() <= m_etaRange.second);
      bool pt_cut = (jet->get_pt() >= m_ptRange.first) and (jet->get_pt() <= m_ptRange.second);
      if ((not eta_cut) or (not pt_cut)) continue;
      if(jet->get_e() < 0) {
        return Fun4AllReturnCodes::EVENT_OK; // currently applied to deal with cold EMCal IB
      }
      m_nComponent.push_back(jet->size_comp());
      m_eta.push_back(jet->get_eta());
      m_phi.push_back(jet->get_phi());
      m_e.push_back(jet->get_e());
      m_pt.push_back(jet->get_pt());
      m_calibpt.push_back(calibjet->get_pt());

      //if (m_calibpt.back() > leadpt) { leadpt = m_calibpt.back(); }
      if (m_pt.back() > leadpt) { leadpt = m_pt.back(); }

      float emcalE = 0;
      float ihcalE = 0;
      float ohcalE = 0;
      float jet_time = 0;
      float jet_unweighttime = 0;
      float jet_time_count = 0;
      for (auto comp: jet->get_comp_vec())
      {
        TowerInfo *tower;
        unsigned int channel = comp.second;
        if (comp.first == 14 || comp.first == 29 || comp.first == 28) {
          tower = towersEM3->get_tower_at_channel(channel);
          if(!tower) { continue; }
          emcalE += tower->get_energy();
          float jet_tower_e = tower->get_energy();
          float jet_tower_time = tower->get_time();
          if (jet_tower_e > 0.1) {
            jet_time += jet_tower_time * jet_tower_e;
            jet_unweighttime += jet_tower_time;
            jet_time_count += jet_tower_e;
          }
        }    
        if (comp.first == 15 ||  comp.first == 30 || comp.first == 26)
        {
          tower = towersIH3->get_tower_at_channel(channel);
          if(!tower) { continue; }
          ihcalE += tower->get_energy();
          float jet_tower_e = tower->get_energy();
          float jet_tower_time = tower->get_time();
          if (jet_tower_e > 0.1) {
            jet_time += jet_tower_time * jet_tower_e;
            jet_unweighttime += jet_tower_time;
            jet_time_count += jet_tower_e;
          }
        }

        if (comp.first == 16 ||  comp.first == 31 || comp.first == 27)
        {
          tower = towersOH3->get_tower_at_channel(channel);
          if(!tower) { continue; }
          ohcalE += tower->get_energy();
          float jet_tower_e = tower->get_energy();
          float jet_tower_time = tower->get_time();
          if (jet_tower_e > 0.1) {
            jet_time += jet_tower_time * jet_tower_e;
            jet_unweighttime += jet_tower_time;
            jet_time_count += jet_tower_e;
          }
        }
      }

      m_jetEmcalE.push_back(emcalE);
      m_jetIhcalE.push_back(ihcalE);
      m_jetOhcalE.push_back(ohcalE);
      if (jet_time_count > 0) {
        jet_time = jet_time / jet_time_count;
        jet_unweighttime = jet_unweighttime / jet_time_count;
      }
      else {
        jet_time = -9999;
        jet_unweighttime = -9999;
      }
      m_jettime.push_back(jet_time);
      m_jetunweighttime.push_back(jet_unweighttime);

      if (ohcalE/m_e.back() > 0.9) {
        std::cout << "event " << m_event << " emcalE " << emcalE << " ohcalE " << ohcalE << " totalE " << emcalE + ihcalE + ohcalE << " Jet E " << m_e.back() << std::endl;
        //return Fun4AllReturnCodes::EVENT_OK;
      }
      m_nJet++;
    }

  //get truth jets
  float truthleadpt = 0;
  if(m_doTruthJets)
    {
      m_nTruthJet = 0;
      //for (JetMap::Iter iter = jetsMC->begin(); iter != jetsMC->end(); ++iter)
      for (auto truthjet : *jetsMC)
	{
	  //Jet* truthjet = iter->second;
	    
	  bool eta_cut = (truthjet->get_eta() >= m_etaRange.first) and (truthjet->get_eta() <= m_etaRange.second);
	  bool pt_cut = (truthjet->get_pt() >= m_ptRange.first) and (truthjet->get_pt() <= m_ptRange.second);
	  if ((not eta_cut) or (not pt_cut)) continue;
	  m_truthNComponent.push_back(truthjet->size_comp());
	  m_truthEta.push_back(truthjet->get_eta());
	  m_truthPhi.push_back(truthjet->get_phi());
	  m_truthE.push_back(truthjet->get_e());
	  m_truthPt.push_back(truthjet->get_pt());
    if (m_truthPt.back() > truthleadpt) { truthleadpt = m_truthPt.back(); }
	  m_nTruthJet++;
	}
    }

  if (((m_doTruthJets && m_doTruthLeadPtCut && truthleadpt < m_truthLeadPtCut) || !m_doTruthJets) && (m_doLeadPtCut && leadpt < m_leadPtCut)) {
    return Fun4AllReturnCodes::EVENT_OK; 
  }
  
  //get seed jets
  if(m_doSeeds)
    {
      for (auto jet : *seedjetsraw)
	{
	  int passesCut = jet->get_property(seedjetsraw->property_index(Jet::PROPERTY::prop_SeedItr));
	  m_eta_rawseed.push_back(jet->get_eta());
	  m_phi_rawseed.push_back(jet->get_phi());
	  m_e_rawseed.push_back(jet->get_e());
	  m_pt_rawseed.push_back(jet->get_pt());
	  m_rawseed_cut.push_back(passesCut);
	}
      
      for (auto jet : *seedjetssub)
	{
	  int passesCut = jet->get_property(seedjetssub->property_index(Jet::PROPERTY::prop_SeedItr));
	  m_eta_subseed.push_back(jet->get_eta());
	  m_phi_subseed.push_back(jet->get_phi());
	  m_e_subseed.push_back(jet->get_e());
	  m_pt_subseed.push_back(jet->get_pt());
	  m_subseed_cut.push_back(passesCut);
	}
    }

  //grab the gl1 data
  Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "14001");
  if (m_doTriggerCut && !gl1PacketInfo)
    {
      std::cout << PHWHERE << "caloTreeGen::process_event: GL1Packet node is missing. Output related to this node will be empty" << std::endl;
    }

  if (gl1PacketInfo) {
    bool jettrig = false;
    uint64_t triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++) {
	    bool trig_decision = ((triggervec & 0x1U) == 0x1U);
      if (trig_decision) {
        m_triggerVector.push_back(i);
        if ((i >= 16 && i <= 23) || (i >= 32 && i <= 35)) { jettrig = true; }
      }
	    triggervec = triggervec >> 1U;
	  }

    uint64_t livetriggervec = gl1PacketInfo->getLiveVector();
    for (int i = 0; i < 64; i++) {
      bool trig_decision = ((livetriggervec & 0x1U) == 0x1U);
      if (trig_decision) { m_liveTriggerVector.push_back(i); }
      livetriggervec = livetriggervec >> 1U;
    }

    if (m_doTriggerCut && !jettrig) {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  if (m_doTowers) {
    if(EMtowers) { 
      m_emcaln = 0;
      int nchannels = 24576; 
      for(int i=0; i<nchannels; ++i) {
        TowerInfo *tower = EMtowers->get_tower_at_channel(i); //get EMCal tower
        //if(!tower->get_isGood()) { continue; }
        int key = EMtowers->encode_key(i);
        int etabin = EMtowers->getTowerEtaBin(key);
        int phibin = EMtowers->getTowerPhiBin(key);
        m_emcale[m_emcaln] = tower->get_energy(); 
        m_emcalchi2[m_emcaln] = tower->get_chi2();
        if (!tower->get_isBadChi2() && tower->get_chi2() > 10000) {
          std::cout << "event " << m_event << " ieta " << etabin << " iphi " << phibin << " chi2 " << tower->get_chi2() << " e " << tower->get_energy() << std::endl;
        }
        const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, etabin, phibin);
        RawTowerGeom *geom = tower_geomEM->get_tower_geometry(geomkey); 
        TVector3 tower_pos;
        tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - m_zvtx);
        if ((m_doMBDeff || m_doMBDeffsyst) && fabs(m_zvtx) > 200) {
          tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - 0.0); 
        }
        m_emcaleta[m_emcaln] = tower_pos.Eta();
        m_emcalphi[m_emcaln] = tower_pos.Phi();
        m_emcalieta[m_emcaln] = etabin;
        m_emcaliphi[m_emcaln] = phibin;
        m_emcalstatus[m_emcaln] = tower->get_status();
        //m_emcaltime[m_emcaln] = tower->get_time_float();
        m_emcaln++;
      }
    }

    if(IHtowers) { 
      m_ihcaln = 0;
      int nchannels = 1536; 
      for(int i=0; i<nchannels; ++i) {
        TowerInfo *tower = IHtowers->get_tower_at_channel(i); //get IHCal tower
        //if(!tower->get_isGood()) { continue; }
        int key = IHtowers->encode_key(i);
        int etabin = IHtowers->getTowerEtaBin(key);
        int phibin = IHtowers->getTowerPhiBin(key);
        m_ihcale[m_ihcaln] = tower->get_energy(); 
        m_ihcalchi2[m_ihcaln] = tower->get_chi2();
        const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, etabin, phibin);
        RawTowerGeom *geom = tower_geom->get_tower_geometry(geomkey); 
        TVector3 tower_pos;
        tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - m_zvtx);
        if ((m_doMBDeff || m_doMBDeffsyst) && fabs(m_zvtx) > 200) {
          tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - 0.0); 
        }
        m_ihcaleta[m_ihcaln] = tower_pos.Eta();
        m_ihcalphi[m_ihcaln] = tower_pos.Phi();
        m_ihcalieta[m_ihcaln] = etabin;
        m_ihcaliphi[m_ihcaln] = phibin;
        m_ihcalstatus[m_ihcaln] = tower->get_status();
        //m_ihcaltime[m_ihcaln] = tower->get_time_float();
        m_ihcaln++;
      }
    }

    if(OHtowers) { 
      m_ohcaln = 0;
      int nchannels = 1536; 
      for(int i=0; i<nchannels; ++i) {
        TowerInfo *tower = OHtowers->get_tower_at_channel(i); //get OHCal tower
        //if(!tower->get_isGood()) { continue; }
        int key = OHtowers->encode_key(i);
        int etabin = OHtowers->getTowerEtaBin(key);
        int phibin = OHtowers->getTowerPhiBin(key);
        m_ohcale[m_ohcaln] = tower->get_energy(); 
        m_ohcalchi2[m_ohcaln] = tower->get_chi2();
        const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, etabin, phibin);
        RawTowerGeom *geom = tower_geomOH->get_tower_geometry(geomkey); 
        TVector3 tower_pos;
        tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - m_zvtx);
        if ((m_doMBDeff || m_doMBDeffsyst) && fabs(m_zvtx) > 200) {
          tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - 0.0); 
        }
        m_ohcaleta[m_ohcaln] = tower_pos.Eta();
        m_ohcalphi[m_ohcaln] = tower_pos.Phi();
        m_ohcalieta[m_ohcaln] = etabin;
        m_ohcaliphi[m_ohcaln] = phibin;
        m_ohcalstatus[m_ohcaln] = tower->get_status();
        //m_ohcaltime[m_ohcaln] = tower->get_time_float();
        m_ohcaln++;
      }
    }
  }
  if (m_doEmcalClusters) {
    if (emcalclusters) {
      RawClusterContainer::Map clusterMap = emcalclusters->getClustersMap();
      m_emcal_clsmult = 0;
      for(auto entry : clusterMap){
        RawCluster* cluster = entry.second;
        CLHEP::Hep3Vector origin(0, 0, m_zvtx);
        m_emcal_cluster_e[m_emcal_clsmult] = cluster->get_energy();
        m_emcal_cluster_eta[m_emcal_clsmult] = RawClusterUtility::GetPseudorapidity(*cluster, origin);
        m_emcal_cluster_phi[m_emcal_clsmult] = RawClusterUtility::GetAzimuthAngle(*cluster, origin);
        m_emcal_clsmult++;
        if (m_emcal_clsmult == 10000) { break; }
      } 
    }
  }
  if (m_doTopoclusters) {
    if (topoclusters) {
      RawClusterContainer::Map clusterMap = topoclusters->getClustersMap();
      m_clsmult = 0;
      float cluster_vertex = m_zvtx;
      if ((m_doMBDeff || m_doMBDeffsyst) && fabs(m_zvtx) > 200) { cluster_vertex = 0.0; }
      for(auto entry : clusterMap){
        RawCluster* cluster = entry.second;
        CLHEP::Hep3Vector origin(0, 0, cluster_vertex);
        m_cluster_e[m_clsmult] = cluster->get_energy();
        m_cluster_eta[m_clsmult] = RawClusterUtility::GetPseudorapidity(*cluster, origin);
        m_cluster_phi[m_clsmult] = RawClusterUtility::GetAzimuthAngle(*cluster, origin);
        m_cluster_ntowers[m_clsmult] = (int)cluster->getNTowers();
        int m_clstower = 0;
        for (const auto& [tower_id, tower_e] : cluster->get_towermap())
        {
          /*
            if (m_cluster_e[m_clsmult] < -1.0) {
              std::cout << "Tower ID: " << tower_id
                        << " | Tower E: " << tower_e
                        << " | Calorimeter ID: " << static_cast<int>(RawTowerDefs::decode_caloid(tower_id))
                        << " | Index1: " << RawTowerDefs::decode_index1(tower_id)
                        << " | Index2: " << RawTowerDefs::decode_index2(tower_id)
                        << std::endl;
            }
            */
            m_cluster_tower_calo[m_clsmult][m_clstower] = static_cast<int>(RawTowerDefs::decode_caloid(tower_id));
            m_cluster_tower_ieta[m_clsmult][m_clstower]  = RawTowerDefs::decode_index1(tower_id);
            m_cluster_tower_iphi[m_clsmult][m_clstower]  = RawTowerDefs::decode_index2(tower_id);
            m_cluster_tower_e[m_clsmult][m_clstower] = tower_e;
            /*
            if (m_cluster_e[m_clsmult] < -1.0) {
              std::cout << "Tower ID: " << tower_id
                        << " | Tower E: " << tower_e << " " << m_cluster_tower_e[m_clsmult][m_clstower]
                        << " | Calorimeter ID: " << m_cluster_tower_calo[m_clsmult][m_clstower]
                        << " | Index1: " << m_cluster_tower_ieta[m_clsmult][m_clstower]
                        << " | Index2: " << m_cluster_tower_iphi[m_clsmult][m_clstower]
                        << std::endl;
            }
            */
            m_clstower++;
        }
        //std::cout << std::endl;
        m_clsmult++;
        if (m_clsmult == 10000) { break; }
      } 
      //std::cout << "number of clusters: " << m_clsmult << std::endl;
    }

    if (topoclusters2) {
      RawClusterContainer::Map clusterMap = topoclusters2->getClustersMap();
      m_clsmult2 = 0;
      float cluster_vertex = m_zvtx;
      if ((m_doMBDeff || m_doMBDeffsyst) && fabs(m_zvtx) > 200) { cluster_vertex = 0.0; }
      for(auto entry : clusterMap){
        RawCluster* cluster = entry.second;
        CLHEP::Hep3Vector origin(0, 0, cluster_vertex);
        m_cluster2_e[m_clsmult2] = cluster->get_energy();
        m_cluster2_eta[m_clsmult2] = RawClusterUtility::GetPseudorapidity(*cluster, origin);
        m_cluster2_phi[m_clsmult2] = RawClusterUtility::GetAzimuthAngle(*cluster, origin);
        m_cluster2_ntowers[m_clsmult2] = (int)cluster->getNTowers();
        int m_clstower = 0;
        for (const auto& [tower_id, tower_e] : cluster->get_towermap())
        {
            m_cluster2_tower_calo[m_clsmult2][m_clstower] = static_cast<int>(RawTowerDefs::decode_caloid(tower_id));
            m_cluster2_tower_ieta[m_clsmult2][m_clstower]  = RawTowerDefs::decode_index1(tower_id);
            m_cluster2_tower_iphi[m_clsmult2][m_clstower]  = RawTowerDefs::decode_index2(tower_id);
            m_cluster2_tower_e[m_clsmult2][m_clstower] = tower_e;
            m_clstower++;
        }
        m_clsmult2++;
        if (m_clsmult2 == 10000) { break; }
      } 
    }

    if (topoclusters4) {
      RawClusterContainer::Map clusterMap = topoclusters4->getClustersMap();
      m_clsmult4 = 0;
      float cluster_vertex = m_zvtx;
      if ((m_doMBDeff || m_doMBDeffsyst) && fabs(m_zvtx) > 200) { cluster_vertex = 0.0; }
      for(auto entry : clusterMap){
        RawCluster* cluster = entry.second;
        CLHEP::Hep3Vector origin(0, 0, cluster_vertex);
        m_cluster4_e[m_clsmult4] = cluster->get_energy();
        m_cluster4_eta[m_clsmult4] = RawClusterUtility::GetPseudorapidity(*cluster, origin);
        m_cluster4_phi[m_clsmult4] = RawClusterUtility::GetAzimuthAngle(*cluster, origin);
        m_cluster4_ntowers[m_clsmult4] = (int)cluster->getNTowers();
        int m_clstower = 0;
        for (const auto& [tower_id, tower_e] : cluster->get_towermap())
        {
            m_cluster4_tower_calo[m_clsmult4][m_clstower] = static_cast<int>(RawTowerDefs::decode_caloid(tower_id));
            m_cluster4_tower_ieta[m_clsmult4][m_clstower]  = RawTowerDefs::decode_index1(tower_id);
            m_cluster4_tower_iphi[m_clsmult4][m_clstower]  = RawTowerDefs::decode_index2(tower_id);
            m_cluster4_tower_e[m_clsmult4][m_clstower] = tower_e;
            m_clstower++;
        }
        m_clsmult4++;
        if (m_clsmult4 == 10000) { break; }
      } 
    }
  }

  if (m_doTruth) {
    PHG4TruthInfoContainer::Range range = truthinfo->GetSPHENIXPrimaryParticleRange();
    truthpar_n = 0;
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter) {
      // Get truth particle
      const PHG4Particle *truth = iter->second;
      if (!truth) { std::cout << "missing particle" << std::endl; continue; }
      //if (!truthinfo->is_primary(truth)) continue;
      /// Get this particles momentum, etc.
      float eta = atanh(truth->get_pz() / sqrt(truth->get_px()*truth->get_px()+truth->get_py()*truth->get_py()+truth->get_pz()*truth->get_pz()));
      if (fabs(eta) > 1.1) { continue; }
      truthpar_pt[truthpar_n] = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
      truthpar_pz[truthpar_n] = truth->get_pz();
      truthpar_e[truthpar_n] = truth->get_e();
      truthpar_phi[truthpar_n] = atan2(truth->get_py(), truth->get_px());
      truthpar_eta[truthpar_n] = atanh(truth->get_pz() / sqrt(truth->get_px()*truth->get_px()+truth->get_py()*truth->get_py()+truth->get_pz()*truth->get_pz()));
      if (truthpar_eta[truthpar_n] != truthpar_eta[truthpar_n]) truthpar_eta[truthpar_n] = -999; // check for nans
      truthpar_pid[truthpar_n] = truth->get_pid();
      truthpar_n++;
      if(truthpar_n > 99999)
      {
        return Fun4AllReturnCodes::EVENT_OK;
      }
    }

    PHG4TruthInfoContainer::Range range2 = truthinfo->GetPrimaryParticleRange();
    old_truthpar_n = 0;
    for (PHG4TruthInfoContainer::ConstIterator iter2 = range2.first; iter2 != range2.second; ++iter2) {
      // Get truth particle
      const PHG4Particle *truth = iter2->second;
      if (!truth) { std::cout << "missing particle" << std::endl; continue; }
      //if (!truthinfo->is_primary(truth)) continue;
      /// Get this particles momentum, etc.
      float eta = atanh(truth->get_pz() / sqrt(truth->get_px()*truth->get_px()+truth->get_py()*truth->get_py()+truth->get_pz()*truth->get_pz()));
      if (fabs(eta) > 1.1) { continue; }
      old_truthpar_pt[old_truthpar_n] = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
      old_truthpar_pz[old_truthpar_n] = truth->get_pz();
      old_truthpar_e[old_truthpar_n] = truth->get_e();
      old_truthpar_phi[old_truthpar_n] = atan2(truth->get_py(), truth->get_px());
      old_truthpar_eta[old_truthpar_n] = atanh(truth->get_pz() / sqrt(truth->get_px()*truth->get_px()+truth->get_py()*truth->get_py()+truth->get_pz()*truth->get_pz()));
      if (old_truthpar_eta[old_truthpar_n] != old_truthpar_eta[old_truthpar_n]) old_truthpar_eta[old_truthpar_n] = -999; // check for nans
      old_truthpar_pid[old_truthpar_n] = truth->get_pid();
      old_truthpar_n++;
      if(old_truthpar_n > 99999)
      {
        return Fun4AllReturnCodes::EVENT_OK;
      }
    }
  }

  if(m_doTracks) {
    SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    // Check whether SvtxTrackMap is present in the node tree or not
    if (!trackmap) {
      std::cout << PHWHERE << "SvtxTrackMap node is missing, can't collect tracks" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    /*
    SvtxVertexMap *vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    if (!vertexmap) {
      std::cout << PHWHERE << "SvtxVertexMap node is missing, can't collect tracks" << std::endl;
      return;
    }
    */

    /*
    // EvalStack for truth track matching
    if(!m_svtxEvalStack){ m_svtxEvalStack = new SvtxEvalStack(topNode); }
    m_svtxEvalStack->next_event(topNode);
    SvtxTrackEval *trackeval = m_svtxEvalStack->get_track_eval();
    */

    m_trkmult = 0;
    // Looping over tracks in an event
    for(auto entry : *trackmap) {
      SvtxTrack *track = entry.second;

      m_tr_p[m_trkmult] = track->get_p();
      m_tr_pt[m_trkmult] = track->get_pt();
      m_tr_eta[m_trkmult] = track->get_eta();
      m_tr_phi[m_trkmult] = track->get_phi();
    
      m_tr_charge[m_trkmult] = track->get_charge();
      m_tr_chisq[m_trkmult] = track->get_chisq();
      m_tr_ndf[m_trkmult] = track->get_ndf();
      m_tr_quality[m_trkmult] = track->get_quality();
      m_tr_crossing[m_trkmult] = track->get_crossing();

      TrackSeed *silseed = track->get_silicon_seed();
      TrackSeed *tpcseed = track->get_tpc_seed();

      m_tr_nmaps[m_trkmult] = 0;
      m_tr_nintt[m_trkmult] = 0;
      m_tr_ntpc[m_trkmult] = 0;

      if(tpcseed) {
        for (TrackSeed::ConstClusterKeyIter iter = tpcseed->begin_cluster_keys(); iter != tpcseed->end_cluster_keys(); ++iter) {
          TrkrDefs::cluskey cluster_key = *iter;
          unsigned int layer = TrkrDefs::getLayer(cluster_key);
          if (_nlayers_maps > 0 && layer < _nlayers_maps) {
            m_tr_nmaps[m_trkmult]++;
          }
          if (_nlayers_intt > 0 && layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt) {
            m_tr_nintt[m_trkmult]++;
          }
          if (_nlayers_tpc > 0 && layer >= (_nlayers_maps + _nlayers_intt) && layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc)) {
            m_tr_ntpc[m_trkmult]++;
          }
        }
      } 

      if(silseed) {
        for (TrackSeed::ConstClusterKeyIter iter = silseed->begin_cluster_keys(); iter != silseed->end_cluster_keys(); ++iter) {
          TrkrDefs::cluskey cluster_key = *iter;
          unsigned int layer = TrkrDefs::getLayer(cluster_key);
          if (_nlayers_maps > 0 && layer < _nlayers_maps) {
            m_tr_nmaps[m_trkmult]++;
          }
          if (_nlayers_intt > 0 && layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt) {
            m_tr_nintt[m_trkmult]++;
          }
          if (_nlayers_tpc > 0 && layer >= (_nlayers_maps + _nlayers_intt) && layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc)) {
            m_tr_ntpc[m_trkmult]++;
            }
          }
        }

      // Project tracks to CEMC
      if(track->find_state(m_cemcRadius) != track->end_states()){
        m_tr_cemc_eta[m_trkmult] = calculateProjectionEta( track->find_state(m_cemcRadius)->second);
        m_tr_cemc_phi[m_trkmult] = calculateProjectionPhi( track->find_state(m_cemcRadius)->second);
      }
      else{
        m_tr_cemc_eta[m_trkmult] = -999;
        m_tr_cemc_phi[m_trkmult] = -999;
      }
      
      // Project tracks to inner HCAL
      if(track->find_state(m_ihcalRadius) != track->end_states()){
        m_tr_ihcal_eta[m_trkmult] = calculateProjectionEta( track->find_state(m_ihcalRadius)->second);
        m_tr_ihcal_phi[m_trkmult] = calculateProjectionPhi( track->find_state(m_ihcalRadius)->second);
      }
      else{
        m_tr_ihcal_eta[m_trkmult] = -999;
        m_tr_ihcal_phi[m_trkmult] = -999;
      }

      // Project tracks to outer HCAL
      if(track->find_state(m_ohcalRadius) != track->end_states()){
        m_tr_ohcal_eta[m_trkmult] = calculateProjectionEta( track->find_state(m_ohcalRadius)->second);
        m_tr_ohcal_phi[m_trkmult] = calculateProjectionPhi( track->find_state(m_ohcalRadius)->second);
      }
      else{
        m_tr_ohcal_eta[m_trkmult] = -999;
        m_tr_ohcal_phi[m_trkmult] = -999;
      }

      if(track->find_state(m_cemcOuterRadius) != track->end_states()){
        m_tr_outer_cemc_eta[m_trkmult] = calculateProjectionEta( track->find_state(m_cemcOuterRadius)->second);
        m_tr_outer_cemc_phi[m_trkmult] = calculateProjectionPhi( track->find_state(m_cemcOuterRadius)->second);
      }
      else{
        m_tr_outer_cemc_eta[m_trkmult] = -999;
        m_tr_outer_cemc_phi[m_trkmult] = -999;
      }
      
      // Project tracks to inner HCAL
      if(track->find_state(m_ihcalOuterRadius) != track->end_states()){
        m_tr_outer_ihcal_eta[m_trkmult] = calculateProjectionEta( track->find_state(m_ihcalOuterRadius)->second);
        m_tr_outer_ihcal_phi[m_trkmult] = calculateProjectionPhi( track->find_state(m_ihcalOuterRadius)->second);
      }
      else{
        m_tr_outer_ihcal_eta[m_trkmult] = -999;
        m_tr_outer_ihcal_phi[m_trkmult] = -999;
      }

      // Project tracks to outer HCAL
      if(track->find_state(m_ohcalOuterRadius) != track->end_states()){
        m_tr_outer_ohcal_eta[m_trkmult] = calculateProjectionEta( track->find_state(m_ohcalOuterRadius)->second);
        m_tr_outer_ohcal_phi[m_trkmult] = calculateProjectionPhi( track->find_state(m_ohcalOuterRadius)->second);
      }
      else{
        m_tr_outer_ohcal_eta[m_trkmult] = -999;
        m_tr_outer_ohcal_phi[m_trkmult] = -999;
      }

      /*
      // Matching SvtxTracks to G4 truth
      PHG4Particle *truthtrack = trackeval->max_truth_particle_by_nclusters(track);

      if(truthtrack){
        m_tr_truth_pid[m_trkmult] =truthtrack->get_pid();
        m_tr_truth_is_primary[m_trkmult] =truthinfo->is_primary(truthtrack);
    
        m_tr_truth_e[m_trkmult] =truthtrack->get_e();

        float px = truthtrack->get_px();
        float py = truthtrack->get_py();
        float pz = truthtrack->get_pz();

        m_tr_truth_pt[m_trkmult] = sqrt(px*px+py*py);
        m_tr_truth_phi[m_trkmult] = atan2(py, px);
        m_tr_truth_eta[m_trkmult] = atanh(pz/sqrt(px*px+py*py+pz*pz));
        m_tr_truth_track_id[m_trkmult] = truthtrack->get_track_id();
      } else{
        m_tr_truth_pid[m_trkmult] = -999;
        m_tr_truth_is_primary[m_trkmult] = -999;
        m_tr_truth_e[m_trkmult] = -999;
        m_tr_truth_pt[m_trkmult] = -999;
        m_tr_truth_eta[m_trkmult] = -999;
        m_tr_truth_phi[m_trkmult] = -999;
        m_tr_truth_track_id[m_trkmult] = -999;
      }
      */
      m_trkmult++;
    }

    /*
    m_vtxmult = 0;
    for(auto entry : *vertexmap)
    {
      SvtxVertex *vertex = entry.second;
      
      m_vertex_id[m_vtxmult] = vertex->get_id();
      m_vx[m_vtxmult] = vertex->get_x();
      m_vy[m_vtxmult] = vertex->get_y();
      m_vz[m_vtxmult] = vertex->get_z();
      m_vtxmult++;
    }
    */
  }
  
  //fill the tree
  m_T->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::ResetEvent(PHCompositeNode *topNode)
{
  //std::cout << "InclusiveJet::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  m_nComponent.clear();
  m_eta.clear();
  m_phi.clear();
  m_e.clear();
  m_pt.clear();
  m_calibpt.clear();

  m_jetEmcalE.clear();
  m_jetIhcalE.clear();
  m_jetOhcalE.clear();

  m_jettime.clear();
  m_jetunweighttime.clear();

  m_truthNComponent.clear();
  m_truthEta.clear();
  m_truthPhi.clear();
  m_truthE.clear();
  m_truthPt.clear();
  m_truthdR.clear();

  m_eta_subseed.clear();
  m_phi_subseed.clear();
  m_e_subseed.clear();
  m_pt_subseed.clear();
  m_subseed_cut.clear();

  m_eta_rawseed.clear();
  m_phi_rawseed.clear();
  m_e_rawseed.clear();
  m_pt_rawseed.clear();
  m_rawseed_cut.clear();
  
  m_triggerVector.clear();
  m_liveTriggerVector.clear();

  m_svtxVector.clear();
  m_svtxBcoVector.clear();
  m_siliconVector.clear();
  m_siliconBcoVector.clear();

  m_zvtx = -9999;
  m_zsvtx = -9999;
  m_ztruthvtx = -9999;
  m_zsiliconvtx = -9999;
  m_emcaln = 0;
  m_ihcaln = 0;
  m_ohcaln = 0;
  m_clsmult = 0;
  m_clsmult2 = 0;
  m_emcal_clsmult = 0;
  m_trkmult = 0;
  m_mbd_t0 = 0;
  m_mbd_ts = 0;
  m_mbd_tn = 0;

  for (int i = 0; i < 24576; i++) {
    m_emcale[i] = 0;
    m_emcalieta[i] = 0;
    m_emcaliphi[i] = 0;
    m_emcaleta[i] = 0;
    m_emcalphi[i] = 0;
    m_emcalchi2[i] = 0;
    m_emcalstatus[i] = 0;
    //m_emcaltime[i] = 0;
  }

  for (int i = 0; i < 1536; i++) {
    m_ihcale[i] = 0;
    m_ihcalieta[i] = 0;
    m_ihcaliphi[i] = 0;
    m_ihcaleta[i] = 0;
    m_ihcalphi[i] = 0;
    m_ihcalchi2[i] = 0;
    m_ihcalstatus[i] = 0;
    //m_ihcaltime[i] = 0;
    m_ohcale[i] = 0;
    m_ohcalieta[i] = 0;
    m_ohcaliphi[i] = 0;
    m_ohcaleta[i] = 0;
    m_ohcalphi[i] = 0;
    m_ohcalchi2[i] = 0;
    m_ohcalstatus[i] = 0;
    //m_ohcaltime[i] = 0;
  }

  for (int i = 0; i < 100000; i++) {
    truthpar_pt[i] = 0;
    truthpar_pz[i] = 0;
    truthpar_e[i] = 0;
    truthpar_phi[i] = 0;
    truthpar_eta[i] = 0;
    truthpar_pid[i] = 0;
  }

  for (int i = 0; i < 100000; i++) {
    old_truthpar_pt[i] = 0;
    old_truthpar_pz[i] = 0;
    old_truthpar_e[i] = 0;
    old_truthpar_phi[i] = 0;
    old_truthpar_eta[i] = 0;
    old_truthpar_pid[i] = 0;
  }

  for (int i = 0; i < 2000; i++) {
    m_cluster_e[i] = 0;
    m_cluster_eta[i] = 0;
    m_cluster_phi[i] = 0;
    m_cluster_ntowers[i] = 0;
    m_emcal_cluster_e[i] = 0;
    m_emcal_cluster_eta[i] = 0;
    m_emcal_cluster_phi[i] = 0;
    for (int j = 0; j < 500; j++) {
      m_cluster_tower_calo[i][j] = 0;
      m_cluster_tower_ieta[i][j] = 0;
      m_cluster_tower_iphi[i][j] = 0;
      m_cluster_tower_e[i][j] = 0;
    }
    m_cluster2_e[i] = 0;
    m_cluster2_eta[i] = 0;
    m_cluster2_phi[i] = 0;
    m_cluster2_ntowers[i] = 0;
    for (int j = 0; j < 500; j++) {
      m_cluster2_tower_calo[i][j] = 0;
      m_cluster2_tower_ieta[i][j] = 0;
      m_cluster2_tower_iphi[i][j] = 0;
      m_cluster2_tower_e[i][j] = 0;
    }
  }

  for (int i = 0; i < 2000; i++) {
    m_tr_p[i] = 0;
    m_tr_pt[i] = 0;
    m_tr_eta[i] = 0;
    m_tr_phi[i] = 0;
    m_tr_charge[i] = 0;
    m_tr_chisq[i] = 0;
    m_tr_ndf[i] = 0;
    m_tr_nintt[i] = 0;
    m_tr_nmaps[i] = 0;
    m_tr_ntpc[i] = 0;
    m_tr_quality[i] = 0;
    m_tr_crossing[i] = 0;
    m_tr_cemc_eta[i] = 0; // Projection of track to calorimeters
    m_tr_cemc_phi[i] = 0;
    m_tr_ihcal_eta[i] = 0;
    m_tr_ihcal_phi[i] = 0;
    m_tr_ohcal_eta[i] = 0;
    m_tr_ohcal_phi[i] = 0;
    m_tr_outer_cemc_eta[i] = 0;
    m_tr_outer_cemc_phi[i] = 0;
    m_tr_outer_ihcal_eta[i] = 0;
    m_tr_outer_ihcal_phi[i] = 0;
    m_tr_outer_ohcal_eta[i] = 0;
    m_tr_outer_ohcal_phi[i] = 0;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

float InclusiveJet::calculateProjectionEta(SvtxTrackState* projectedState){
  float x = projectedState->get_x();// - initialState->get_x();
  float y = projectedState->get_y();// - initialState->get_y();
  float z = projectedState->get_z();// - initialState->get_z();

  float theta = acos(z / sqrt(x*x + y*y + z*z) );

  return -log( tan(theta/2.0) );
}

float InclusiveJet::calculateProjectionPhi(SvtxTrackState* projectedState){
  float x = projectedState->get_x();// - initialState->get_x();
  float y = projectedState->get_y();// - initialState->get_y();
 
  return atan2(y,x);
}

//____________________________________________________________________________..
int InclusiveJet::EndRun(const int runnumber)
{
  std::cout << "InclusiveJet::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::End(PHCompositeNode *topNode)
{
  std::cout << "InclusiveJet::End - Output to " << m_outputFileName << std::endl;
  PHTFileServer::get().cd(m_outputFileName);

  m_T->Write();
  std::cout << "InclusiveJet::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::Reset(PHCompositeNode *topNode)
{
  std::cout << "InclusiveJet::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void InclusiveJet::Print(const std::string &what) const
{
  std::cout << "InclusiveJet::Print(const std::string &what) const Printing info for " << what << std::endl;
}
