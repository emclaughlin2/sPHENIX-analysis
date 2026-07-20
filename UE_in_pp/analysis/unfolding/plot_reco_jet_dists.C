#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TPad.h>
#include <vector>
#include <string>

// 1.4903e-6 - Jet 12 cross section -> scale 1.4903
// 3.997e-6 - Jet 10 cross section -> scale 3.997

int rgb[20][3] = {{0, 0, 0}, {230, 25, 75}, {60, 180, 75}, {255, 225, 25}, {0, 130, 200}, {245, 130, 48}, {145, 30, 180}, {70, 240, 240}, {240, 50, 230}, {210, 245, 60}, {250, 190, 212}, {0, 128, 128}, {220, 190, 255}, {170, 110, 40}, {128, 128, 128}, {128, 0, 0}, {128, 128, 0}, {255, 215, 180}, {0, 0, 128}, {34, 139, 34}};
int colors[20];
for (int ij = 0; ij < 20; ij++) colors[ij] = TColor::GetColor(rgb[ij][0], rgb[ij][1], rgb[ij][2]);

void plot_reco_jet_dists() {
  std::vector<std::string> inputFiles = {
      //"analysis_sim_run21_output/output_efrac_bkg_cut_sim_iter_3.root",
      //"analysis_sim_run28_output/test_run21_output_analysis_clean_none_bkg_cut_sim_iter_3.root",
      //"analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3_run21_event_weights_jet10_20_30_50_zvtx_lt_30_eta_lt_0.7.root",
      "analysis_sim_run28_output/output_none_bkg_cut_sim_iter_1.root",
      //"analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3.root",
      "analysis_sim_ppg09_test/output_none_bkg_cut_sim_iter_1.root",
      //"analysis_sim_ppg09_test/output_dijet_bkg_cut_sim_iter_1_no_jet5.root",
      "analysis_sim_run21_test_output/output_none_bkg_cut_sim_iter_1.root",
      //"analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3_run28_event_weights_alljets_zvtx_lt_60_calo_accept_eta.root",
      //"analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut_zvtx_lt_30_eta_lt_0.7.root",
      //"analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root",
      //"analysis_data_run28_output/output_pu_correct_data_vz_lt_30_dijet_bkg_cut.root",
      //"analysis_sim_run28_output/output_zvertex_noreweight_none_bkg_cut_sim_iter_1.root",
      //"analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut_zvtx_lt_60_calo_accept_eta_no_jet_bkg_eff_correction.root"

  };

  std::vector<float> scales = {1.4903 * 0.81608598, 1.4903 * 0.81608598, 3.997, 0, 0};

  std::vector<std::string> labels = {
      "MC Run28 dijet w Jet 5",
      //"MC Run28 dijet w/o Jet 5",
      "MC Run28 Hanpu ttree",
      //"MC Run28 Hanpu ttree w/o Jet 5",
      "MC Run21 dijet",
      //"MC |zvtx| < 30 cm & |eta| < 0.7",
      //"MC |zvtx| < 30 cm & |eta| < 0.7 check",
      //"MC |zvtx| < 30 cm & |eta| < 0.7 check run28",
      //"MC |zvtx| < 60 cm & vtx reweight",
      //"MC |zvtx| < 60 cm & |eta| in calo accept",
      //"Data |zvtx| < 30 cm & |eta| < 0.7",
      "Data dijet",
      "Data |zvtx| < 30 cm dijet"};

  std::vector<TH1D*> truth_projX, truth_projY;
  std::vector<TH1D*> measure_projX, measure_projY;
  std::vector<TH1D*> zvertex;
  std::vector<TH1D*> truth_jet, truth_et;
  std::vector<TH1D*> calib_jet, calib_et;
  std::vector<TH1D*> uncalib_jet;

  // Process simulation files (0-2)
  for (int i = 0; i < 3; ++i)
  {
    TFile* file = TFile::Open(inputFiles[i].c_str());
    if (!file || file->IsZombie())
    {
      std::cerr << "Failed to open file: " << inputFiles[i] << std::endl;
      if (file) file->Close();
      continue;
    }

    TH1D* hZVtx = (TH1D*) file->Get("h_zvertex");
    if (hZVtx)
    {
      TH1D* hZVtx_clone = (TH1D*) hZVtx->Clone(("zvertex_" + std::to_string(i)).c_str());
      hZVtx_clone->SetDirectory(0);
      zvertex.push_back(hZVtx_clone);
    }
    else
    {
      std::cerr << "Histogram h_zvertex NOT FOUND in " << inputFiles[i] << std::endl;
    }

    // Verify the response object and its histograms exist and have the expected types
    // TObject* obj = file->Get("h_respmatrix_calib_dijet_reweight_trim_10");
    TObject* obj = file->Get("h_respmatrix_calib_dijet");
    RooUnfoldResponse* response_check = dynamic_cast<RooUnfoldResponse*>(obj);
    if (!response_check)
    {
      std::cerr << "Response 'h_respmatrix_calib_dijet_reweight_trim_10' not found or wrong type in " << inputFiles[i] << std::endl;
      file->Close();
      continue;
    }

    TH2D* hTruth_check = dynamic_cast<TH2D*>(response_check->Htruth());
    TH2D* hMeasured_check = dynamic_cast<TH2D*>(response_check->Hmeasured());
    if (!hTruth_check || !hMeasured_check)
    {
      std::cerr << "Htruth or Hmeasured missing/wrong type in response in " << inputFiles[i] << std::endl;
      file->Close();
      continue;
    }

    TH1D* px = (TH1D*) hTruth_check->ProjectionX();
    TH1D* py = (TH1D*) hTruth_check->ProjectionY();
    TH1D* mx = (TH1D*) hMeasured_check->ProjectionX();
    TH1D* my = (TH1D*) hMeasured_check->ProjectionY();

    TH1D* px_clone = (TH1D*) px->Clone(("truth_projX_" + std::to_string(i)).c_str());
    TH1D* py_clone = (TH1D*) py->Clone(("truth_projY_" + std::to_string(i)).c_str());
    TH1D* mx_clone = (TH1D*) mx->Clone(("measure_projX_" + std::to_string(i)).c_str());
    TH1D* my_clone = (TH1D*) my->Clone(("measure_projY_" + std::to_string(i)).c_str());
    px_clone->SetDirectory(0);
    py_clone->SetDirectory(0);
    mx_clone->SetDirectory(0);
    my_clone->SetDirectory(0);

    truth_projX.push_back(px_clone);
    truth_projY.push_back(py_clone);
    measure_projX.push_back(mx_clone);
    measure_projY.push_back(my_clone);

    TH1D* htruthjet = (TH1D*) file->Get("h_truth_jet_pt_tight");
    //TH1D* htruthjet = (TH1D*) file->Get("h_lead_truth_spectra_record");
    TH1D* htruthjet_clone = (TH1D*) htruthjet->Clone(("truth_jet_" + std::to_string(i)).c_str());
    htruthjet_clone->SetDirectory(0);
    truth_jet.push_back(htruthjet_clone);
    TH1D* hcalibjet = (TH1D*) file->Get("h_calib_jet_pt_tight");
    //TH1D* hcalibjet = (TH1D*) file->Get("h_lead_spectra_record");
    TH1D* hcalibjet_clone = (TH1D*) hcalibjet->Clone(("calib_jet_" + std::to_string(i)).c_str());
    hcalibjet_clone->SetDirectory(0);
    calib_jet.push_back(hcalibjet_clone);
    TH1D* htruthet = (TH1D*) file->Get("h_truth_calo_et_tight");
    TH1D* htruthet_clone = (TH1D*) htruthet->Clone(("truth_et_" + std::to_string(i)).c_str());
    htruthet_clone->SetDirectory(0);
    truth_et.push_back(htruthet_clone);
    TH1D* hcalibet = (TH1D*) file->Get("h_calib_calo_et_tight");
    TH1D* hcalibet_clone = (TH1D*) hcalibet->Clone(("calib_et_" + std::to_string(i)).c_str());
    hcalibet_clone->SetDirectory(0);
    calib_et.push_back(hcalibet_clone);
    TH1D* huncalibjet = (TH1D*) file->Get("h_uncalib_lead_spectra_record");
    TH1D* huncalibjet_clone = (TH1D*) huncalibjet->Clone(("uncalib_jet_" + std::to_string(i)).c_str());
    huncalibjet_clone->SetDirectory(0);
    uncalib_jet.push_back(huncalibjet_clone);
    file->Close();
    }
    
    // Process data files (3-4)
    /*
    for (int i = 3; i < 5; ++i) {
        TFile* file = TFile::Open(inputFiles[i].c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Failed to open file: " << inputFiles[i] << std::endl;
            continue;
        }

        TH1D* hZVtx = (TH1D*)file->Get("h_zvertex");
        if (hZVtx) {
            TH1D* hZVtx_clone = (TH1D*)hZVtx->Clone(("zvertex_" + std::to_string(i)).c_str());
            hZVtx_clone->SetDirectory(0);
            zvertex.push_back(hZVtx_clone);
        } else {
            std::cerr << "Histogram h_zvertex NOT FOUND in " << inputFiles[i] << std::endl;
        }

        TH2D* hData = (TH2D*)file->Get("h_calibjet_pt_dijet_eff");
        if (!hData) {
            std::cerr << "Histogram h_calibjet_pt_dijet_eff NOT FOUND in " 
                    << inputFiles[i] << std::endl;
            file->Close();
            continue;
        }

        TH1D* mx = (TH1D*)hData->ProjectionX();
        TH1D* my = (TH1D*)hData->ProjectionY();
        TH1D* mx_clone = (TH1D*)mx->Clone(("measure_projX_" + std::to_string(i)).c_str());
        TH1D* my_clone = (TH1D*)my->Clone(("measure_projY_" + std::to_string(i)).c_str());  
        mx_clone->SetDirectory(0);
        my_clone->SetDirectory(0);

        measure_projX.push_back(mx_clone);
        measure_projY.push_back(my_clone);
        
        TH1D* hcalibjet = (TH1D*)file->Get("h_calib_jet_pt_tight");
        TH1D* hcalibjet_clone = (TH1D*)hcalibjet->Clone(("calib_jet_" + std::to_string(i)).c_str());
        hcalibjet_clone->SetDirectory(0);
        calib_jet.push_back(hcalibjet_clone);
        TH1D* hcalibet = (TH1D*)file->Get("h_calib_calo_et_tight");
        TH1D* hcalibet_clone = (TH1D*)hcalibet->Clone(("calib_et_" + std::to_string(i)).c_str());
        hcalibet_clone->SetDirectory(0);
        calib_et.push_back(hcalibet_clone);
        
        file->Close();
    }
    */    

    std::cout << "Finished loading histograms." << std::endl;
    std::cout << "Number of measure projections X: " << measure_projX.size() << std::endl;
    std::cout << "Number of measure projections Y: " << measure_projY.size() << std::endl;
    std::cout << "Number of truth projections X: " << truth_projX.size() << std::endl;
    std::cout << "Number of truth projections Y: " << truth_projY.size() << std::endl;
    std::cout << "Labels size: " << labels.size() << std::endl;

    for (size_t i = 0; i < measure_projX.size(); i++) {
        std::cout << "measure_projX[" << i << "] = " << measure_projX[i]  << " measure_projY[" << i << "] = " << measure_projY[i]  << " label = " << labels[i] << std::endl;
    }

    for (size_t i = 0; i < zvertex.size(); i++) {
        zvertex[i]->SetTitle(labels[i].c_str());
        zvertex[i]->Scale(1.0 / zvertex[i]->Integral());
    }

    TCanvas* zcanvas = new TCanvas("zcanvas", "Z Vertex Distribution", 1200, 800);
    for (size_t i = 0; i < zvertex.size(); ++i) {
        zvertex[i]->SetLineColor(colors[i]);
        zvertex[i]->SetLineWidth(2);
        zvertex[i]->SetMarkerStyle(20);
        zvertex[i]->SetMarkerColor(colors[i]);
        zvertex[i]->GetYaxis()->SetRangeUser(0, 0.03);
        if (i == 0) {
            zvertex[i]->Draw();
        } else {
            zvertex[i]->Draw("SAME");
        }
    }
    zcanvas->BuildLegend();
    zcanvas->SaveAs("compare_run28_run21_test_plots/zvertex_dists.png");

    for (size_t i = 0; i < truth_jet.size(); i++) {
        truth_jet[i]->SetTitle(labels[i].c_str());
        if (scales[i] != 0) {
            truth_jet[i]->Scale(scales[i]);
        } else {
            truth_jet[i]->Scale(1.0 / truth_jet[i]->Integral());
        }
        truth_et[i]->SetTitle(labels[i].c_str());
        truth_et[i]->Scale(1.0 / truth_et[i]->Integral());
    }

    TCanvas* tcanvas = new TCanvas("tcanvas", "Truth Jet Distribution", 1200, 900);
    TPad* tpad_top = new TPad("tpad_top", "truth jet top pad", 0.0, 0.30, 1.0, 1.0);
    TPad* tpad_bottom = new TPad("tpad_bottom", "truth jet ratio pad", 0.0, 0.0, 1.0, 0.30);
    tpad_top->SetBottomMargin(0.02);
    tpad_bottom->SetTopMargin(0.02);
    tpad_bottom->SetBottomMargin(0.30);
    tpad_top->Draw();
    tpad_bottom->Draw();

    tpad_top->cd(); 
    for (size_t i = 0; i < truth_jet.size(); ++i) {
        truth_jet[i]->SetLineColor(colors[i]);
        truth_jet[i]->SetLineWidth(2);  
        truth_jet[i]->SetMarkerStyle(20);
        truth_jet[i]->SetMarkerColor(colors[i]);
        truth_jet[i]->SetXTitle("p_{T}^{Truth} [GeV]");
        //truth_jet[i]->GetYaxis()->SetRangeUser(0, 0.03);
        if (i == 0) {
            truth_jet[i]->Draw();
        } else {
            truth_jet[i]->Draw("SAME");
        }
    }
    tpad_top->SetLogy(1);
    tpad_top->BuildLegend();

    tpad_bottom->cd();
    TH1D* truth_ref = (TH1D*)truth_jet[0]->Clone("truth_ref");
    truth_ref->SetDirectory(0);
    std::vector<TH1D*> truth_ratio_hists;
    for (size_t i = 0; i < truth_jet.size(); ++i) {
        TH1D* ratio = (TH1D*)truth_jet[i]->Clone(("truth_jet_ratio_" + std::to_string(i)).c_str());
        ratio->SetDirectory(0);
        ratio->Divide(truth_ref);
        ratio->SetLineColor(colors[i]);
        ratio->SetMarkerColor(colors[i]);
        ratio->SetMarkerStyle(20);
        ratio->GetYaxis()->SetTitle("Ratio");
        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetTitleOffset(0.45);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetXaxis()->SetTitle("p_{T}^{Truth} [GeV]");
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->GetYaxis()->SetRangeUser(0.4, 1.6);
        if (i == 0) {
            ratio->Draw("EP");
        } else {
            ratio->Draw("EP SAME");
        }
        truth_ratio_hists.push_back(ratio);
    }
    TLine* truth_unity = new TLine(truth_ref->GetXaxis()->GetXmin(), 1.0,
                                   truth_ref->GetXaxis()->GetXmax(), 1.0);
    truth_unity->SetLineStyle(2);
    truth_unity->Draw("SAME");
    tcanvas->SaveAs("compare_run28_run21_test_plots/truth_jet_dists_xsection_scale.png");
    //tcanvas->SaveAs("compare_run28_run21_test_plots/truth_jet_dists.png");
/*
    TCanvas* tecanvas = new TCanvas("tecanvas", "Truth ET Distribution", 1200, 800);
    for (size_t i = 0; i < truth_et.size(); ++i) {
        truth_et[i]->SetLineColor(colors[i]);
        truth_et[i]->SetLineWidth(2);
        truth_et[i]->SetMarkerStyle(20);
        truth_et[i]->SetMarkerColor(colors[i]);
        truth_et[i]->SetXTitle("E_{T}^{Truth} [GeV]");
        //truth_et[i]->GetYaxis()->SetRangeUser(0, 0.03);
        if (i == 0) {
            truth_et[i]->Draw();
        } else {
            truth_et[i]->Draw("SAME");
        }
    }
    tecanvas->BuildLegend();
    tecanvas->SaveAs("compare_run28_run21_test_plots/truth_et_dists.png");
*/
    for (size_t i = 0; i < calib_jet.size(); i++) {
        calib_jet[i]->SetTitle(labels[i].c_str());
        if (scales[i] != 0) {
            calib_jet[i]->Scale(scales[i]);
        } else {
            calib_jet[i]->Scale(1.0 / calib_jet[i]->Integral());
        }
        calib_et[i]->SetTitle(labels[i].c_str());
        calib_et[i]->Scale(1.0 / calib_et[i]->Integral());
    }

    TCanvas* ccanvas = new TCanvas("ccanvas", "Calib Jet Distribution", 1200, 900);
    TPad* cpad_top = new TPad("cpad_top", "calib jet top pad", 0.0, 0.30, 1.0, 1.0);
    TPad* cpad_bottom = new TPad("cpad_bottom", "calib jet ratio pad", 0.0, 0.0, 1.0, 0.30);
    cpad_top->SetBottomMargin(0.02);
    cpad_bottom->SetTopMargin(0.02);
    cpad_bottom->SetBottomMargin(0.30);
    cpad_top->Draw();
    cpad_bottom->Draw();

    cpad_top->cd();
    for (size_t i = 0; i < calib_jet.size(); ++i) {
        calib_jet[i]->SetLineColor(colors[i]);
        calib_jet[i]->SetLineWidth(2);
        calib_jet[i]->SetMarkerStyle(20);
        calib_jet[i]->SetMarkerColor(colors[i]);
        calib_jet[i]->SetXTitle("p_{T}^{Calib} [GeV]");
        //calib_jet[i]->GetYaxis()->SetRangeUser(, 0.6);
        //calib_jet[i]->GetXaxis()->SetRangeUser(20, 65);
        if (i == 0) {
            calib_jet[i]->Draw();
        } else {
            calib_jet[i]->Draw("SAME");
        }
    }
    cpad_top->SetLogy(1);
    cpad_top->BuildLegend();

    cpad_bottom->cd();
    TH1D* calib_ref = (TH1D*)calib_jet[0]->Clone("calib_ref");
    calib_ref->SetDirectory(0);
    std::vector<TH1D*> calib_ratio_hists;
    for (size_t i = 0; i < calib_jet.size(); ++i) {
        TH1D* ratio = (TH1D*)calib_jet[i]->Clone(("calib_jet_ratio_" + std::to_string(i)).c_str());
        ratio->SetDirectory(0);
        ratio->Divide(calib_ref);
        ratio->SetLineColor(colors[i]);
        ratio->SetMarkerColor(colors[i]);
        ratio->SetMarkerStyle(20);
        ratio->GetYaxis()->SetTitle("Ratio");
        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetTitleOffset(0.45);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetXaxis()->SetTitle("p_{T}^{Calib} [GeV]");
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->GetYaxis()->SetRangeUser(0.4, 1.6);
        if (i == 0) {
            ratio->Draw("EP");
        } else {
            ratio->Draw("EP SAME");
        }
        calib_ratio_hists.push_back(ratio);
    }
    TLine* calib_unity = new TLine(calib_ref->GetXaxis()->GetXmin(), 1.0,
                                   calib_ref->GetXaxis()->GetXmax(), 1.0);
    calib_unity->SetLineStyle(2);
    calib_unity->Draw("SAME");
    ccanvas->SaveAs("compare_run28_run21_test_plots/calib_jet_dists_xsection_scale.png");
    //ccanvas->SaveAs("compare_run28_run21_test_plots/calib_jet_dists.png");

    for (size_t i = 0; i < uncalib_jet.size(); i++) {
        uncalib_jet[i]->SetTitle(labels[i].c_str());
        if (scales[i] != 0) {
            uncalib_jet[i]->Scale(scales[i]);
        } else {
            uncalib_jet[i]->Scale(1.0 / uncalib_jet[i]->Integral());
        }
    }

    TCanvas* unccanvas = new TCanvas("unccanvas", "Calib Jet Distribution", 1200, 900);
    TPad* uncpad_top = new TPad("uncpad_top", "uncalib jet top pad", 0.0, 0.30, 1.0, 1.0);
    TPad* uncpad_bottom = new TPad("uncpad_bottom", "uncalib jet ratio pad", 0.0, 0.0, 1.0, 0.30);
    uncpad_top->SetBottomMargin(0.02);
    uncpad_bottom->SetTopMargin(0.02);
    uncpad_bottom->SetBottomMargin(0.30);
    uncpad_top->Draw();
    uncpad_bottom->Draw();

    uncpad_top->cd();
    for (size_t i = 0; i < uncalib_jet.size(); ++i) {
        uncalib_jet[i]->SetLineColor(colors[i]);
        uncalib_jet[i]->SetLineWidth(2);
        uncalib_jet[i]->SetMarkerStyle(20);
        uncalib_jet[i]->SetMarkerColor(colors[i]);
        uncalib_jet[i]->SetXTitle("p_{T}^{Calib} [GeV]");
        //uncalib_jet[i]->GetYaxis()->SetRangeUser(, 0.6);
        //uncalib_jet[i]->GetXaxis()->SetRangeUser(20, 65);
        if (i == 0) {
            uncalib_jet[i]->Draw();
        } else {
            uncalib_jet[i]->Draw("SAME");
        }
    }
    uncpad_top->SetLogy(1);
    uncpad_top->BuildLegend();

    uncpad_bottom->cd();
    TH1D* uncalib_ref = (TH1D*)uncalib_jet[0]->Clone("uncalib_ref");
    uncalib_ref->SetDirectory(0);
    std::vector<TH1D*> uncalib_ratio_hists;
    for (size_t i = 0; i < uncalib_jet.size(); ++i) {
        TH1D* ratio = (TH1D*)uncalib_jet[i]->Clone(("uncalib_jet_ratio_" + std::to_string(i)).c_str());
        ratio->SetDirectory(0);
        ratio->Divide(uncalib_ref);
        ratio->SetLineColor(colors[i]);
        ratio->SetMarkerColor(colors[i]);
        ratio->SetMarkerStyle(20);
        ratio->GetYaxis()->SetTitle("Ratio");
        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetTitleOffset(0.45);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetXaxis()->SetTitle("p_{T}^{Calib} [GeV]");
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->GetYaxis()->SetRangeUser(0.4, 1.6);
        if (i == 0) {
            ratio->Draw("EP");
        } else {
            ratio->Draw("EP SAME");
        }
        uncalib_ratio_hists.push_back(ratio);
    }
    TLine* uncalib_unity = new TLine(uncalib_ref->GetXaxis()->GetXmin(), 1.0,
                                   uncalib_ref->GetXaxis()->GetXmax(), 1.0);
    uncalib_unity->SetLineStyle(2);
    uncalib_unity->Draw("SAME");
    unccanvas->SaveAs("compare_run28_run21_test_plots/uncalib_jet_dists_xsection_scale.png");
    //ccanvas->SaveAs("compare_run28_run21_test_plots/uncalib_jet_dists.png");
    /*
    TCanvas* cecanvas = new TCanvas("cecanvas", "Calib ET Distribution", 1200, 800);
    for (size_t i = 0; i < calib_et.size(); ++i) {
        calib_et[i]->SetLineColor(colors[i]);
        calib_et[i]->SetLineWidth(2);
        calib_et[i]->SetMarkerStyle(20);
        calib_et[i]->SetMarkerColor(colors[i]);
        calib_et[i]->SetXTitle("E_{T}^{Reco} [GeV]");
        //calib_et[i]->GetYaxis()->SetRangeUser(0, 0.03);
        if (i == 0) {
            calib_et[i]->Draw();
        } else {
            calib_et[i]->Draw("SAME");
        }
    }
    cecanvas->SetLogy(1);
    cecanvas->BuildLegend();
    cecanvas->SaveAs("compare_run28_run21_test_plots/calib_et_dists.png");

    // Plot histograms with labels
    for (size_t i = 0; i < measure_projX.size(); ++i) {
        std::cout << i << " " << labels[i] << " " << measure_projX[i] << std::endl;
        measure_projX[i]->SetTitle(labels[i].c_str());
        measure_projY[i]->SetTitle(labels[i].c_str());
        //measure_projX[i]->Scale(scales[i]);
        //measure_projY[i]->Scale(scales[i]);
        measure_projX[i]->Scale(1.0 / measure_projX[i]->Integral());
        measure_projY[i]->Scale(1.0 / measure_projY[i]->Integral());
    }

    for (size_t i = 0; i < truth_projX.size(); ++i) {
        truth_projX[i]->SetTitle(labels[i].c_str());
        truth_projY[i]->SetTitle(labels[i].c_str());
        //truth_projX[i]->Scale(scales[i]);
        //truth_projY[i]->Scale(scales[i]);
        truth_projX[i]->Scale(1.0 / truth_projX[i]->Integral());
        truth_projY[i]->Scale(1.0 / truth_projY[i]->Integral());
    }

    TCanvas* canvasX = new TCanvas("canvasX", "Projection X", 1200, 800);
    for (size_t i = 0; i < measure_projX.size(); ++i) {
        //measure_projX[i]->SetLineColor(i + 1);
        //measure_projX[i]->SetLineWidth(2);
        measure_projX[i]->SetMarkerStyle(20);
        measure_projX[i]->SetMarkerColor(colors[i]);
        measure_projX[i]->GetYaxis()->SetRangeUser(0,1);
        if (i == 0) {
            measure_projX[i]->Draw("p");
        } else {
            measure_projX[i]->Draw("p,SAME");
        }
    }
    canvasX->BuildLegend();
    //canvasX->SaveAs("compare_run28_run21_test_plots/reco_dists_projX_xsection_scale.png");
    canvasX->SaveAs("compare_run28_run21_test_plots/reco_dists_projX.png");

    TCanvas* canvasY = new TCanvas("canvasY", "Projection Y", 1200, 800);
    for (size_t i = 0; i < measure_projY.size(); ++i) {
        //measure_projY[i]->SetLineColor(i + 1);
        //measure_projY[i]->SetLineWidth(2);
        measure_projY[i]->SetMarkerStyle(20);
        measure_projY[i]->SetMarkerColor(colors[i]);
        if (i == 0) {
            measure_projY[i]->Draw("p");
        } else {
            measure_projY[i]->Draw("p,SAME");
        }
    }
    canvasY->BuildLegend();
    //canvasY->SaveAs("compare_run28_run21_test_plots/reco_dists_projY_xsection_scale.png");
    canvasY->SaveAs("compare_run28_run21_test_plots/reco_dists_projY.png");

    TCanvas* truthcanvasX = new TCanvas("truthcanvasX", "Projection X", 1200, 800);
    for (size_t i = 0; i < truth_projX.size(); ++i) {
        //truth_projX[i]->SetLineColor(i + 1);
        //truth_projX[i]->SetLineWidth(2);
        truth_projX[i]->SetMarkerStyle(20);
        truth_projX[i]->SetMarkerColor(colors[i]);
        if (i == 0) {
            truth_projX[i]->Draw("p");
        } else {
            truth_projX[i]->Draw("p,SAME");
        }
    }
    truthcanvasX->BuildLegend();
    //truthcanvasX->SaveAs("compare_run28_run21_test_plots/truth_dists_projX_xsection_scale.png");
    truthcanvasX->SaveAs("compare_run28_run21_test_plots/truth_dists_projX.png");

    TCanvas* truthcanvasY = new TCanvas("truthcanvasY", "Projection Y", 1200, 800);
    for (size_t i = 0; i < truth_projY.size(); ++i) {
        //truth_projY[i]->SetLineColor(i + 1);
        //truth_projY[i]->SetLineWidth(2);
        truth_projY[i]->SetMarkerStyle(20);
        truth_projY[i]->SetMarkerColor(colors[i]);
        if (i == 0) {
            truth_projY[i]->Draw("p");
        } else {
            truth_projY[i]->Draw("p,SAME");
        }
    }
    truthcanvasY->BuildLegend();
    //truthcanvasY->SaveAs("compare_run28_run21_test_plots/truth_dists_projY_xsection_scale.png");
    truthcanvasY->SaveAs("compare_run28_run21_test_plots/truth_dists_projY.png");
    */
}