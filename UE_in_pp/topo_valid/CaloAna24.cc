

#include "CaloAna24.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

// Fun4All
#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// ROOT stuff
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TTree.h>

// For clusters and geometry
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

// Tower stuff
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// GL1 Information
#include <ffarawobjects/Gl1Packet.h>

// for cluster vertex correction
#include <CLHEP/Geometry/Point3D.h>

// for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdGeom.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <TLorentzVector.h>
//____________________________________________________________________________..
CaloAna24::CaloAna24(const std::string &name) : SubsysReco(name)
{
  std::cout << "CaloAna24::CaloAna24(const std::string &name) Calling ctor"
            << std::endl;
}

//____________________________________________________________________________..
CaloAna24::~CaloAna24()
{
  std::cout << "CaloAna24::~CaloAna24() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CaloAna24::Init(PHCompositeNode *topNode)
{
  fout = new TFile("caloana.root", "RECREATE");

  T = new TTree("T", "T");

  T->Branch("nrealclusterpair", &m_nrealclusterpair, "nrealclusterpair/I");
  T->Branch("nmixedclusterpair", &m_nmixedclusterpair, "nmixedclusterpair/I");
  T->Branch("nsinglecluster", &m_nsinglecluster, "nsinglecluster/I");

  T->Print();

  T->Branch("singleclusterE", m_singleclusterE, "singleclusterE[nsinglecluster]/F");
  T->Branch("singleclustereta", m_singleclustereta, "singleclustereta[nsinglecluster]/F");
  T->Branch("singleclusterpT", m_singleclusterpT, "singleclusterpT[nsinglecluster]/F");
  T->Branch("singleclusterprob", m_singleclusterprob, "singleclusterprob[nsinglecluster]/F");

  T->Branch("realphoton1E", m_realphoton1E, "realphoton1E[nrealclusterpair]/F");
  //T->Branch("realphoton1prob", m_realphoton1prob, "realphoton1prob[nrealclusterpair]/F");
  T->Branch("realphoton2E", m_realphoton2E, "realphoton2E[nrealclusterpair]/F");
  //T->Branch("realphoton2prob", m_realphoton2prob, "realphoton2prob[nrealclusterpair]/F");
  T->Branch("realpairmass", m_realpairmass, "realpairmass[nrealclusterpair]/F");
  T->Branch("realpairpt", m_realpairpt, "realpairpt[nrealclusterpair]/F");
  T->Branch("realpairdR", m_realpairdR, "realpairdR[nrealclusterpair]/F");
  T->Branch("realpaireta", m_realpaireta, "realpaireta[nrealclusterpair]/F");
  T->Branch("realpairphi", m_realpairphi, "realpairphi[nrealclusterpair]/F");

  
  T->Branch("mixedphoton1E", m_mixedphoton1E, "mixedphoton1E[nmixedclusterpair]/F");
  //T->Branch("mixedphoton1prob", m_mixedphoton1prob, "mixedphoton1prob[nmixedclusterpair]/F");
  T->Branch("mixedphoton2E", m_mixedphoton2E, "mixedphoton2E[nmixedclusterpair]/F");
  //T->Branch("mixedphoton2prob", m_mixedphoton2prob, "mixedphoton2prob[nmixedclusterpair]/F");
  T->Branch("mixedpairmass", m_mixedpairmass, "mixedpairmass[nmixedclusterpair]/F");
  T->Branch("mixedpairpt", m_mixedpairpt, "mixedpairpt[nmixedclusterpair]/F");
  T->Branch("mixedpairdR", m_mixedpairdR, "mixedpairdR[nmixedclusterpair]/F");
  T->Branch("mixedpaireta", m_mixedpaireta, "mixedpaireta[nmixedclusterpair]/F");
  T->Branch("mixedpairphi", m_mixedpairphi, "mixedpairphi[nmixedclusterpair]/F");

  T->Print();

  h_mass_pt = new TH2F("h_mass_pt", "h_mass_pt", 100, 0, 20, 500, 0, 2);

  h_mass_pt_mixed = new TH2F("h_mass_pt_mixed", "h_mass_pt_mixed", 100, 0, 20, 500, 0, 2);

  h_dETdEta = new TProfile("h_dETdEta", "h_dETdEta", 24, -1, 1);

  h_unrebindETdEta = new TProfile("h_unrebindETdEta", "h_unrebindETdEta", 96, -1, 1);

  h_unscaledtriggercount = new TH1D("h_unscaledtriggercount",
                                    "h_unscaledtriggercount", 32, 0, 32);
  h_scaledtriggercount = new TH1D("h_scaledtriggercount",
                                  "h_scaledtriggercount", 32, 0, 32);
  h_MB_ncluster = new TH1D("h_MB_ncluster", "h_MB_ncluster", 100, 0, 100);

  h_MB_clusterprob = new TH1D("h_MB_clusterprob", "h_MB_clusterprob", 100, 0, 1);

  h_MB_clusterpt = new TH1D("h_MB_clusterpt", "h_MB_clusterpt", 100, 0, 20);

  h_MB_clusterpt_prob = new TH2F("h_MB_clusterpt_prob", "h_MB_clusterpt_prob", 100, 0, 20, 100, 0, 1);

  h_MB_maxclusterpt = new TH1D("h_MB_maxclusterpt", "h_MB_maxclusterpt", 100, -1, 20);

  h_MB_vertexz = new TH1D("h_MB_vertexz", "h_MB_vertexz", 100, -50, 50);

  h_MB_eventcount = new TH1D("h_MB_eventcount", "h_MB_eventcount", 4, 0, 4);

  h_clusterET = new TH1D("h_clusterET", "h_clusterET", 100, 0, 30);

  for (int i = 0; i < ntriggerbins; i++)
  {
    for (int j = 0; j < nassociatebins; j++)
    {
      int lowtriggerpt = triggerptbins[i].first;
      int lowassociatept = associateptbins[j].first;

      h_dphi_deta[i][j] = new TH2D(Form("h_dphi_deta_%d_%d", lowtriggerpt, lowassociatept),
                                   Form("h_dphi_deta_%d_%d", lowtriggerpt, lowassociatept),
                                   64, -M_PI / 2 + M_PI / 64., 3. / 2 * M_PI + M_PI / 64., 48, -2.2 + 2.2 / 48, 2.2 - 2.2 / 48);

      h_dphi_deta_mix[i][j] = new TH2D(Form("h_dphi_deta_mix_%d_%d", lowtriggerpt, lowassociatept),
                                       Form("h_dphi_deta_mix_%d_%d", lowtriggerpt, lowassociatept),
                                       64, -M_PI / 2 + M_PI / 64., 3. / 2 * M_PI + M_PI / 64., 48, -2.2 + 2.2 / 48, 2.2 - 2.2 / 48);
    }
  }

  for (int i = 0; i < ntriggerbins; i++)
  {
    rblisttrig[i] = new RollingBuffer(3);
  }

  for (int i = 0; i < nassociatebins; i++)
  {
    rblistassoc[i] = new RollingBuffer(50);
  }

  rblistpi0cluster = new RollingBuffer(50);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloAna24::InitRun(PHCompositeNode *topNode)
{
  std::cout
      << "CaloAna24::InitRun(PHCompositeNode *topNode) Initializing for Run XXX"
      << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloAna24::process_event(PHCompositeNode *topNode)
{

  TowerInfoContainer *CEMC_towers_calib =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWER_CALIB_HISTO_CEMC");

  int CEMCsize = CEMC_towers_calib->size();

  RawTowerGeomContainer *CEMC_geom =
      findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  const int netabins = 24;
  const int netabinsnorebin = 96;
  const float deta = 2. / netabins;
  const float detanorebin = 2. / netabinsnorebin;
  float toweretatotal[netabins] = {0};
  float toweretatotalnorebin[netabinsnorebin] = {0};
  for (int i = 0; i < CEMCsize; i++)
  {
    TowerInfo *tower = CEMC_towers_calib->get_tower_at_channel(i);

    unsigned int towerkey = CEMC_towers_calib->encode_key(i);
    int ieta = CEMC_towers_calib->getTowerEtaBin(towerkey);
    int iphi = CEMC_towers_calib->getTowerPhiBin(towerkey);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(
        RawTowerDefs::convert_name_to_caloid("CEMC"), ieta, iphi);
    float tower_eta = CEMC_geom->get_tower_geometry(key)->get_eta();

    float tower_E = tower->get_energy();
    if (tower_E > 100)
      continue;
    float tower_ET = tower_E / cosh(tower_eta);

    toweretatotal[(ieta / 4)] += tower_ET / deta;
    toweretatotalnorebin[ieta] += tower_ET / detanorebin;
  }

  for (int i = 0; i < netabins; i++)
  {
    float bincenter = h_dETdEta->GetBinCenter(i + 1);
    // check for nan and inf
    if (toweretatotal[i] != toweretatotal[i])
      toweretatotal[i] = 0;
    if (std::isinf(toweretatotal[i]))
      toweretatotal[i] = 0;
    h_dETdEta->Fill(bincenter, toweretatotal[i]);
  }

  for (int i = 0; i < netabinsnorebin; i++)
  {
    float bincenter = h_unrebindETdEta->GetBinCenter(i + 1);
    // check for nan and inf
    if (toweretatotalnorebin[i] != toweretatotalnorebin[i])
      toweretatotalnorebin[i] = 0;
    if (std::isinf(toweretatotalnorebin[i]))
      toweretatotalnorebin[i] = 0;
    h_unrebindETdEta->Fill(bincenter, toweretatotalnorebin[i]);
  }

  float m_vertex = 0;
  std::vector<TLorentzVector> goodcluster;
  std::vector<std::vector<TLorentzVector>> associatedclusterlist;
  associatedclusterlist.resize(nassociatebins);
  std::vector<std::vector<TLorentzVector>> triggerclusterlist;
  triggerclusterlist.resize(ntriggerbins);

  RawClusterContainer *clusterContainer =
      findNode::getClass<RawClusterContainer>(topNode,
                                              "CLUSTERINFO_CEMC");
  if (!clusterContainer)
  {
    std::cout << PHWHERE
              << "caloTreeGen::process_event: CLUSTER_CEMC node "
                 "is missing. "
                 "Output related to this node will be empty"
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator clusterIter;
  float maxclusterpt = -1;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second;
       clusterIter++)
  {
    RawCluster *recoCluster = clusterIter->second;
    float prob = recoCluster->get_prob();

    CLHEP::Hep3Vector vertex(0, 0, m_vertex);
    float E = recoCluster->get_energy();
    float eta = RawClusterUtility::GetPseudorapidity(*recoCluster, vertex);
    float phi = RawClusterUtility::GetAzimuthAngle(*recoCluster, vertex);
    float ET = E/cosh(eta);

    h_clusterET->Fill(ET);
    TLorentzVector clustercore;
    clustercore.SetPtEtaPhiE(fabs(ET), eta, phi, E);

    for (int i = 0; i < ntriggerbins; i++)
    {
      if (ET > triggerptbins[i].first && ET < triggerptbins[i].second)
      {
        triggerclusterlist[i].push_back(clustercore);
      }
    }
    for (int i = 0; i < nassociatebins; i++)
    {
      if (ET > associateptbins[i].first && ET < associateptbins[i].second)
      {
        associatedclusterlist[i].push_back(clustercore);
      }
    }

    if (ET > maxclusterpt)
    {
      maxclusterpt = ET;
    }

    h_MB_clusterpt->Fill(ET);
    if (ET > 1)
    {
      h_MB_clusterprob->Fill(prob);
      h_MB_clusterpt_prob->Fill(ET, prob);
      if (prob > 1)
        std::cout << "Prob: " << prob << std::endl;
      // eta cut for clusters
      if (abs(eta) > 1)
        continue;
      
      if (prob > 0.0)
      {
        TLorentzVector cluster4;
        cluster4.SetPtEtaPhiE(fabs(ET), eta, phi, E);
        goodcluster.push_back(cluster4);
      }
    }
  }
  std::cout << "Number of good clusters: " << goodcluster.size() << std::endl;
  h_MB_ncluster->Fill(goodcluster.size());
  h_MB_maxclusterpt->Fill(maxclusterpt);
  process_cluster(goodcluster);

  rblistpi0cluster->AddFront(goodcluster);

  auto process_mixed_event_pi0 = [&](RollingBuffer *rb)
  {
    // std::cout<<"Processing mixed event pi0"<<std::endl;
    if (rb->GetSize() < 2)
    {
      std::cout << "There must be at least 2 elements in the buffer."
                << std::endl;
      return;
    }
    std::vector<TLorentzVector> goodcluster1 = rb->GetElement(0);
    for (int i = 1; i < rb->GetSize(); i++)
    {
      std::vector<TLorentzVector> goodcluster2 = rb->GetElement(i);
      for (int i = 0; i < (int)goodcluster1.size(); i++)
      {
        for (int j = 0; j < (int)goodcluster2.size(); j++)
        {
          TLorentzVector cluster1 = goodcluster1[i];
          TLorentzVector cluster2 = goodcluster2[j];

          if (!CutPhoton(cluster1, cluster2))
            continue;
          TLorentzVector cluster = cluster1 + cluster2;
          float pt = cluster.Pt();
          float mass = cluster.M();
          float eta = cluster.Eta();
          float dR = cluster1.DeltaR(cluster2);
          float phi = cluster.Phi();

          m_mixedphoton1E[m_nmixedclusterpair] = cluster1.Energy();
          m_mixedphoton2E[m_nmixedclusterpair] = cluster2.Energy();
          m_mixedpairmass[m_nmixedclusterpair] = mass;
          m_mixedpairpt[m_nmixedclusterpair] = pt;
          m_mixedpairdR[m_nmixedclusterpair] = dR;
          m_mixedpaireta[m_nmixedclusterpair] = eta;
          m_mixedpairphi[m_nmixedclusterpair] = phi;

          m_nmixedclusterpair++;
          // std::cout<<"Mixed cluster pair: "<<nmixedclusterpair<<std::endl;
          if (m_nmixedclusterpair >= mixedclusterpairmaxlength)
          {
            std::cout << "calotrkana::process_event(PHCompositeNode *topNode) "
                         "nmixedclusterpair exceeds max length"
                      << std::endl;
            exit(1);
          }

          if (abs(eta) > 0.8)
            continue;

          h_mass_pt_mixed->Fill(pt, mass);
        }
      }
    }
  };
  process_mixed_event_pi0(rblistpi0cluster);

  for (int i = 0; i < ntriggerbins; i++)
  {
    rblisttrig[i]->AddFront(triggerclusterlist[i]);
  }

  for (int i = 0; i < nassociatebins; i++)
  {
    rblistassoc[i]->AddFront(associatedclusterlist[i]);
  }

  // two particle correlation
  for (int i = 0; i < ntriggerbins; i++)
  {
    for (int j = 0; j < nassociatebins; j++)
    {
      auto process_mixed_event_corr = [&](RollingBuffer *rbtrig, RollingBuffer *rbassoc)
      {
        // std::cout<<"Processing mixed event correlation"<<std::endl;
        if (rbassoc->GetSize() < 2)
        {
          std::cout << "There must be at least 2 elements associate buffer."
                    << std::endl;
          return;
        }
        std::vector<TLorentzVector> trigger1 = rbtrig->GetElement(0);
        for (int ia = 1; ia < rbassoc->GetSize(); ia++)
        {
          std::vector<TLorentzVector> associate1 = rbassoc->GetElement(ia);
          for (int it = 0; it < (int)trigger1.size(); it++)
          {
            for (int ja = 0; ja < (int)associate1.size(); ja++)
            {
              TLorentzVector trigger = trigger1[it];
              TLorentzVector associate = associate1[ja];
              float dr = trigger.DeltaR(associate);
              if (dr < 1e-6)
                continue;
              float dphi = trigger.DeltaPhi(associate);
              while (dphi < (-M_PI / 2 + M_PI / 64.))
                dphi += 2 * M_PI;
              float deta = trigger.Eta() - associate.Eta();
              h_dphi_deta_mix[i][j]->Fill(dphi, deta);
            }
          }
        }
      };
      process_mixed_event_corr(rblisttrig[i], rblistassoc[j]);

      for (int k = 0; k < (int)triggerclusterlist[i].size(); k++)
      {
        for (int l = 0; l < (int)associatedclusterlist[j].size(); l++)
        {
          TLorentzVector trigger = triggerclusterlist[i][k];
          TLorentzVector associate = associatedclusterlist[j][l];
          float dr = trigger.DeltaR(associate);
          if (dr < 1e-6)
            continue;
          float dphi = trigger.DeltaPhi(associate);
          while (dphi < (-M_PI / 2 + M_PI / 64.))
            dphi += 2 * M_PI;
          float deta = trigger.Eta() - associate.Eta();
          h_dphi_deta[i][j]->Fill(dphi, deta);
        }
      }
    }
  }
  T->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna24::process_cluster(std::vector<TLorentzVector> goodcluster)
{
  for (int i = 0; i < (int)goodcluster.size(); i++)
  {
    for (int j = 0; j < (int)goodcluster.size(); j++)
    {
      if (i == j)
        continue;
      TLorentzVector cluster1 = goodcluster[i];
      TLorentzVector cluster2 = goodcluster[j];

      if (!CutPhoton(cluster1, cluster2))
        continue;
      TLorentzVector cluster = cluster1 + cluster2;
      float pt = cluster.Pt();
      float mass = cluster.M();
      float eta = cluster.Eta();
      float dR = cluster1.DeltaR(cluster2);
      float phi = cluster.Phi();

      m_realphoton1E[m_nrealclusterpair] = cluster1.Energy();
      m_realphoton2E[m_nrealclusterpair] = cluster2.Energy();
      m_realpairmass[m_nrealclusterpair] = mass;
      m_realpairpt[m_nrealclusterpair] = pt;
      m_realpairdR[m_nrealclusterpair] = dR;
      m_realpaireta[m_nrealclusterpair] = eta;
      m_realpairphi[m_nrealclusterpair] = phi;

      m_nrealclusterpair++;

      if (m_nrealclusterpair >= trueclusterpairmaxlength)
      {
        std::cout << "calotrkana::process_event(PHCompositeNode *topNode) "
                     "nrealclusterpair exceeds max length"
                  << std::endl;
        exit(1);
      }

      if (abs(eta) > 0.8)
        continue;

      h_mass_pt->Fill(pt, mass);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloAna24::End(PHCompositeNode *topNode)
{

  fout->cd();
  T->Write();

  fout->Write();
  fout->Close();

  // flush
  for (int i = 0; i < ntriggerbins; i++)
  {
    delete rblisttrig[i];
  }

  for (int i = 0; i < nassociatebins; i++)
  {
    delete rblistassoc[i];
  }

  delete rblistpi0cluster;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____not sure why this has no use :(____________________________________________..
int CaloAna24::ResetEvent(PHCompositeNode *topNode)
{

  m_nrealclusterpair = 0;
  m_nmixedclusterpair = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void CaloAna24::Print(const std::string &what) const
{
  std::cout << "CaloAna24::Print(const std::string &what) const Printing "
               "info for "
            << what << std::endl;
}

RollingBuffer::RollingBuffer(int maxLength) : maxLength(maxLength) {}

RollingBuffer::~RollingBuffer()
{
  Flush();
}

void RollingBuffer::AddFront(const std::vector<TLorentzVector> &vec)
{
  buffer.insert(buffer.begin(), vec);
  if ((int)buffer.size() > maxLength)
  {
    PopBack();
  }
}

int RollingBuffer::GetSize() const
{
  return buffer.size();
}

void RollingBuffer::Flush()
{
  buffer.clear();
}

std::vector<TLorentzVector> RollingBuffer::GetElement(int index) const
{
  if (index >= 0 && index < (int)buffer.size())
  {
    return buffer[index];
  }
  return std::vector<TLorentzVector>(); // Return an empty vector if index is out of bounds
}

void RollingBuffer::PopBack()
{
  buffer.pop_back();
}
