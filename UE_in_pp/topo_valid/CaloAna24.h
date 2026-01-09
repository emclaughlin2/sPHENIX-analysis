// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOANA24_H
#define CALOANA24_H

#include <fun4all/SubsysReco.h>

#include <TFile.h>
#include <TTree.h>
#include <TH3.h>
#include <TProfile2D.h>
#include <string>
#include <TLorentzVector.h>

class PHCompositeNode;

class RollingBuffer {
public:
  RollingBuffer(int maxLength = 1000);

  // destructors
  ~RollingBuffer();

  void AddFront(const std::vector<TLorentzVector> &vec);
  int GetSize() const;
  void Flush();
  std::vector<TLorentzVector> GetElement(int index) const;

private:
  void PopBack();

  int maxLength;
  std::vector<std::vector<TLorentzVector>> buffer;
};

class CaloAna24 : public SubsysReco {
public:
  CaloAna24(const std::string &name = "CaloAna24");

  ~CaloAna24() override;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int ResetEvent(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

private:
  int ievent = 0;
  TFile *fout;
  TTree *T = nullptr;
  static const int nsingleclustermaxlength = 10000;
  int m_nsinglecluster = 0;
  float m_singleclusterE[nsingleclustermaxlength] = {0};
  float m_singleclustereta[nsingleclustermaxlength] = {0};
  float m_singleclusterpT[nsingleclustermaxlength] = {0};
  float m_singleclusterprob[nsingleclustermaxlength] = {0};


  static const int trueclusterpairmaxlength = 10000;
  
  float m_realphoton1E[trueclusterpairmaxlength] = {0};
  //float m_realphoton1prob[trueclusterpairmaxlength] = {0};
  float m_realphoton2E[trueclusterpairmaxlength] = {0};
  //float m_realphoton2prob[trueclusterpairmaxlength] = {0};
  float m_realpairmass[trueclusterpairmaxlength] = {0};
  float m_realpairpt[trueclusterpairmaxlength] = {0};
  float m_realpairdR[trueclusterpairmaxlength] = {0};
  float m_realpaireta[trueclusterpairmaxlength] = {0};
  float m_realpairphi[trueclusterpairmaxlength] = {0};
  int m_nrealclusterpair = 0;
  static const int mixedclusterpairmaxlength = 100000;
  
  float m_mixedphoton1E[mixedclusterpairmaxlength] = {0};
  //float m_mixedphoton1prob[mixedclusterpairmaxlength] = {0};
  float m_mixedphoton2E[mixedclusterpairmaxlength] = {0};
  //float m_mixedphoton2prob[mixedclusterpairmaxlength] = {0};
  float m_mixedpairmass[mixedclusterpairmaxlength] = {0};
  float m_mixedpairpt[mixedclusterpairmaxlength] = {0};
  float m_mixedpairdR[mixedclusterpairmaxlength] = {0};
  float m_mixedpaireta[mixedclusterpairmaxlength] = {0};
  float m_mixedpairphi[mixedclusterpairmaxlength] = {0};
  int m_nmixedclusterpair = 0;
  
  // 2d Tprofile
  TH2F *h_mass_pt{nullptr};
  TH2F *h_mass_pt_mixed{nullptr};
  TH1D *h_unscaledtriggercount{nullptr};
  TH1D *h_scaledtriggercount{nullptr};
  TH1D *h_MB_clusterprob{nullptr};
  TH1D *h_MB_ncluster{nullptr};
  TH1D *h_MB_clusterpt{nullptr};
  TH2F* h_MB_clusterpt_prob{nullptr};
  TH1D *h_MB_maxclusterpt{nullptr};
  TH1D *h_MB_vertexz{nullptr};
  TH1D *h_MB_eventcount{nullptr};
  TH1D *h_clusterET{nullptr};
  TProfile * h_dETdEta{nullptr};
  TProfile *h_unrebindETdEta{nullptr};
  static const int ntriggerbins = 4;
  static const int nassociatebins = 4;
  TH2D *h_dphi_deta[ntriggerbins][nassociatebins] = {nullptr};
  TH2D *h_dphi_deta_mix[ntriggerbins][nassociatebins] = {nullptr};
  

  RollingBuffer *rblistpi0cluster = {nullptr};
  RollingBuffer *rblisttrig[ntriggerbins] = {nullptr};
  RollingBuffer *rblistassoc[nassociatebins] = {nullptr};


  std::vector<std::pair<float, float> > associateptbins = {
      {1, 2}, {2, 3}, {3,4}, {4,6} };

  std::vector<std::pair<float, float> > triggerptbins = {
      {2, 3}, {3,4}, {4,6}, {6,8}};
  bool initilized = false;

  long long initscaler[32][3] = {0};
  long long currentscaler[32][3] = {0};
  bool scaledtrigger[32] = {false};
  bool livetrigger[32] = {false};
  int nscaledtrigger[32] = {0};
  int nlivetrigger[32] = {0};
  int process_cluster(std::vector<TLorentzVector> goodcluster);
  float DeltaR(TLorentzVector photon1, TLorentzVector photon2) {
    float deta = photon1.PseudoRapidity() - photon2.PseudoRapidity();
    float dphi = abs(photon1.Phi() - photon2.Phi());

    if (dphi > M_PI)
      dphi = 2 * M_PI - dphi;
    float dr = sqrt(deta * deta + dphi * dphi);
    return dr;
  }

  bool CutPhoton(TLorentzVector photon1, TLorentzVector photon2) {
    // lead photon energy cut
    float leadingphotoncut = 1;
    float subleadingphotoncut = 1;
    if(photon1.Energy() < photon2.Energy()) return false;
    float asymcut = 1.0;
    // drcut
    //float dr = DeltaR(photon1, photon2);
    // if (dr < 0.1) return false;
    //if (dr > 0.4) return false;

    if (photon1.Energy() < leadingphotoncut)
      return false;
    if (photon2.Energy() < subleadingphotoncut)
      return false;
    // asym cut
    float asym = abs(photon1.Energy() - photon2.Energy()) /
                 (photon1.Energy() + photon2.Energy());
    if (asym > asymcut)
      return false;
    // cut on the pt of the sum of the two photons
    float ptsum = (photon1 + photon2).Pt();
    if (ptsum < 2)
      return false;
    if (ptsum > 100)
      return false;
    return true;
  }
};

#endif // CALOANA24_H
