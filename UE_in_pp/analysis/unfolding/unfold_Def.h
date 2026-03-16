const int n_hcal_etabin = 24;
const int n_hcal_phibin = 64;

const float PI = TMath::Pi();
const float jet_radius = 0.4;

const float dijet_min_phi = 3*PI/4.;

const float deltaeta = 0.0916667;
const float deltaphi = 0.0981748;
const float secteta = 2.2;
const float sectphi = (2.0*M_PI)/3.0;

// first set of optimized binning
//double calibptbins[] = {18, 22.5, 28, 34, 40, 47, 55, 63.5, 72};
//double truthptbins[] = {14, 18, 22.5, 28, 34, 40, 47, 55, 63.5, 72, 82};
//double calibetbins[] = {0.0, 0.1, 0.31001396179199225, 0.4800130683506723, 0.6800045540707071, 0.8900158955437308, 1.1399858191862862, 1.4400130785669445, 1.8100128184303708, 2.3399763215468474, 3.2499956903596625, 6.182471266853046, 11, 20};
//double truthetbins[] = {0.0, 0.1, 0.44997577667236327, 0.7099975290529982, 0.9600378567654406, 1.2099803190052243, 1.4599559400810458, 1.7199664593967476, 1.9899819546215256, 2.270032089561787, 2.5799696440961695, 2.9099702816398665, 3.2999990541901827, 3.729994228266298, 4.249971156932052, 4.9100057772473065, 5.829997597189914, 7.4100324517596015, 11.029102638932947, 16.137312894715734, 25, 35};

// binning used for 4 jet bin new et bin studies with Run 22 dataset
//double calibptbins[] = {21, 26, 32.5, 40.5, 63.5}; // 50.5
//double truthptbins[] = {17, 21, 26, 32.5, 40.5, 63.5, 82}; // 72
//int calibnpt = sizeof(calibptbins) / sizeof(calibptbins[0]) - 1;
//int truthnpt = sizeof(truthptbins) / sizeof(truthptbins[0]) - 1;
//double calibetbins[] = {0.0, 0.1, 0.44997577667236327, 0.7099975290529982, 1.2099803190052243, 1.7199664593967476, 2.270032089561787, 2.9099702816398665, 3.729994228266298, 4.9100057772473065, 6, 7.3, 15};
//double truthetbins[] = {0.0, 0.1, 0.44997577667236327, 0.7099975290529982, 1.2099803190052243, 1.7199664593967476, 2.270032089561787, 2.9099702816398665, 3.729994228266298, 4.9100057772473065, 6, 7.3, 15, 35};
//int calibnet = sizeof(calibetbins) / sizeof(calibetbins[0]) - 1;
//int truthnet = sizeof(truthetbins) / sizeof(truthetbins[0]) - 1;

// 7.8.25 - new jet and et bins for Run 21 dataset
// new jet pt binning {17.0, 20.0, 22.0, 25.0, 29.5, 35.0, 42.0, 52.0, 64.0, 82.0};
// new calo et binning {0.0, 0.1, 0.64, 1.08, 1.52, 1.97, 2.47, 3.05, 3.74, 4.68, 6.2, 15.0, 35.0};
//double calibptbins[] = {20.0, 22.5, 25.5, 29.5, 35.0, 42.0, 52.0, 64.0}; // 52.0
//double truthptbins[] = {17.0, 20.0, 22.5, 25.5, 29.5, 35.0, 42.0, 52.0, 64.0, 82.0}; // 52.0

// 8.8.25 - update calibetbins to include 2 negative ET bins
//double calibetbins[] = {0.0, 0.1, 0.64, 1.08, 1.52, 1.97, 2.47, 3.05, 3.74, 4.68, 6.2, 15.0};
//double calibetbins[] = {-1.08, -0.64, -0.1, 0.0, 0.1, 0.64, 1.08, 1.52, 1.97, 2.47, 3.05, 3.74, 4.68, 6.2, 15.0};
//double truthetbins[] = {0.0, 0.1, 0.64, 1.08, 1.52, 1.97, 2.47, 3.05, 3.74, 4.68, 6.2, 15.0, 35.0};

//double calibptbins[] = {22, 27, 34, 40.5, 52.0, 63.5}; // 50.5
//double truthptbins[] = {17, 22, 27, 34, 40.5, 52.0, 63.5, 82}; // 72
//double calibetbins[] = {-0.64, -0.1, 0.0, 0.1, 0.64, 1.52, 2.47, 3.74, 6.2, 15.0}; // calib bins set 2
//double truthetbins[] = {0.0, 0.1, 0.64, 1.52, 2.47, 3.74, 6.2, 15.0, 35.0};

// using these bins Sept. 2025 - Dec. 2025
//double calibptbins[] = {21, 26, 32.5, 40.5, 63.5}; // 52.0
//double truthptbins[] = {17, 21, 26, 32.5, 40.5, 63.5, 82}; // 72
//double calibetbins[] = {-1.08, -0.1, 0.0, 0.1, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0}; // calib bins set 1
//double truthetbins[] = {0.0, 0.1, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0, 35.0};


// try bins Dec. 12th 2025
//double calibptbins[] = {21, 26, 32.5, 40.5, 63.5}; // 52.0
//double truthptbins[] = {17, 21, 26, 32.5, 40.5, 63.5, 82}; // 72
//double calibetbins[] = {-1.08, -0.1, 0.0, 0.1, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0}; // calib bins set 1
//double truthetbins[] = {0.0, 0.5, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0, 35.0};

// try bins Feb. 25th 2025
double calibptbins[] = {21, 26, 32.5, 40.5, 63.5}; // 52.0
double truthptbins[] = {17, 21, 26, 32.5, 40.5, 63.5, 82}; // 72
double calibetbins[] = {-1.08, -0.1, 0.1, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0}; // calib bins set 1
double truthetbins[] = {0.0, 0.5, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0, 35.0};


//double calibetbins[] = {-1.0, -0.1, 0.0, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0}; // yeonju request 
//double truthetbins[] = {0.0, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 35.0}; // yeonju request 

int calibnet = sizeof(calibetbins) / sizeof(calibetbins[0]) - 1;
int truthnet = sizeof(truthetbins) / sizeof(truthetbins[0]) - 1;

int calibnpt = sizeof(calibptbins) / sizeof(calibptbins[0]) - 1;
int truthnpt = sizeof(truthptbins) / sizeof(truthptbins[0]) - 1;

// other binnings
//double calibptbins[] = {15, 19, 23.5, 28.5, 34, 40, 47, 54.5, 63, 72}; // edited to not include MB dataset
//double truthptbins[] = {7, 11, 15, 19, 23.5, 28.5, 34, 40, 47, 54.5, 63, 72, 82}; // edited to not include MB dataset
//1.237,1.213,1.193,1.176,1.175,1.16, 1.156,1.143
//1.25, 1.244,1.231,1.176,1.175,1.17, 1.154,1.134
//double calibetbins[] = {0, 0.1, 0.3116779512076485, 0.38832204879235144, 0.6085863393124066, 0.7882980254420496, 1.008011325310775, 1.1887648668876865, 1.4074941056668524, 1.681781795903305, 2.0062606953561173, 2.373979086620967, 2.979359992203088, 3.8947394990204995, 7.221553415466785, 19.11578012478359};
//double truthetbins[] = {0, 0.1, 0.4526558845425177, 0.6502244308615936, 1.050339564333189, 1.247996705376193, 1.6496178605282503, 1.847070810367734, 2.248378298259465, 2.4477982462134174, 2.943660447418454, 3.247786360560215, 3.835005941330415, 4.339726386570874, 5.213822858113433, 6.27999143815229, 8.565044669116304, 15.074950655966962, 25, 35};
// other ATLAS binning {0.0,0.191,0.382,0.573,0.7639,0.9549,1.146,1.337,1.528,1.719,1.91,2.101,2.292,2.483,2.674,2.865,3.82};
//double calibetbins[] = {0.0, 0.1, 0.31001396179199225, 0.4800130683506723, 0.6800045540707071, 0.8900158955437308, 1.1399858191862862, 1.4400130785669445, 1.8100128184303708, 2.3399763215468474, 3.2499956903596625, 6.182471266853046, 20};
//double truthetbins[] = {0.0, 0.1, 0.44997577667236327, 0.7099975290529982, 0.9600378567654406, 1.2099803190052243, 1.4599559400810458, 1.7199664593967476, 1.9899819546215256, 2.270032089561787, 2.5799696440961695, 2.9099702816398665, 3.2999990541901827, 3.729994228266298, 4.249971156932052, 4.9100057772473065, 5.829997597189914, 7.3, 15, 35};

double Jet5GeV_cross_section = 1.3878e-4; // 1.3878e+0.8 pb
double Jet12GeV_cross_section = 1.4903e-6; // 1.4903e+06 pb
double Jet20GeV_cross_section = 6.2623e-8; // 6.2623e+04 pb
double Jet30GeV_cross_section = 2.5298e-9; // 2.5298e+03 pb
double Jet40GeV_cross_section = 1.3553e-10; // 1.3553e+02 pb
double Jet50GeV_cross_section = 7.3113e-12; // 7.3113 pb
double Jet60GeV_cross_section = 3.3261e-13; // 3.3261e-01 pb

double Herwig_Jet10GeV_cross_section = 4.252e-6; // 4.252e+06 pb
double Herwig_Jet30GeV_cross_section = 2.2831e-9; // 2.2831e+03

double Herwig_Jet10GeV_scale = Herwig_Jet10GeV_cross_section/(double)(Herwig_Jet10GeV_cross_section);
double Herwig_Jet30GeV_scale = Herwig_Jet30GeV_cross_section/(double)(Herwig_Jet10GeV_cross_section);

// Dec. 11 2025 CORRECT scaling uses just the cross section
double Jet5GeV_scale = Jet5GeV_cross_section/Jet12GeV_cross_section;
double Jet12GeV_scale = Jet12GeV_cross_section/Jet12GeV_cross_section;
double Jet20GeV_scale = Jet20GeV_cross_section/Jet12GeV_cross_section;
double Jet30GeV_scale = Jet30GeV_cross_section/Jet12GeV_cross_section;
double Jet40GeV_scale = Jet40GeV_cross_section/Jet12GeV_cross_section;
double Jet50GeV_scale = Jet50GeV_cross_section/Jet12GeV_cross_section;
double Jet60GeV_scale = Jet60GeV_cross_section/Jet12GeV_cross_section;

float get_deta(float eta1, float eta2) {
  return eta1 - eta2;
}

float get_dphi(float phi1, float phi2) {
  float dphi1 = phi1 - phi2;
  float dphi2 = phi1 - phi2 + 2*TMath::Pi();
  float dphi3 = phi1 - phi2 - 2*TMath::Pi();
  if (fabs(dphi1) > fabs(dphi2)) {
    dphi1 = dphi2;
  }
  if (fabs(dphi1) > fabs(dphi3)) {
    dphi1 = dphi3;
  }
  return dphi1;
}

float get_dR(float eta1, float phi1, float eta2, float phi2) {
  float deta = get_deta(eta1, eta2);
  float dphi = get_dphi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

float get_emcal_mineta_zcorrected(float zvertex) {
  float minz_EM = -130.23;
  float radius_EM = 93.5;
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float maxz_EM = 130.23;
  float radius_EM = 93.5;
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float minz_IH = -170.299;
  float radius_IH = 127.503;
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float maxz_IH = 170.299;
  float radius_IH = 127.503;
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float minz_OH = -301.683;
  float radius_OH = 225.87;
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float maxz_OH = 301.683;
  float radius_OH = 225.87;
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}
