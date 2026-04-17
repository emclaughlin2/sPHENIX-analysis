#!/usr/bin/env python3
import sys
import subprocess

# Steps definition with optional working directories
steps = {
    1: [
        ("condor_submit analyze_sim_1_condor.sub", None)
    ],
    2: [
        #("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet5.root output_dijet_bkg_cut_sim_iter_1_jet5_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet12.root output_dijet_bkg_cut_sim_iter_1_jet12_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet20.root output_dijet_bkg_cut_sim_iter_1_jet20_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet30.root output_dijet_bkg_cut_sim_iter_1_jet30_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet40.root output_dijet_bkg_cut_sim_iter_1_jet40_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet50.root output_dijet_bkg_cut_sim_iter_1_jet50_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet60.root output_dijet_bkg_cut_sim_iter_1_jet60_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1.root output_dijet_bkg_cut_sim_iter_1_jet12.root output_dijet_bkg_cut_sim_iter_1_jet20.root output_dijet_bkg_cut_sim_iter_1_jet30.root output_dijet_bkg_cut_sim_iter_1_jet40.root output_dijet_bkg_cut_sim_iter_1_jet50.root output_dijet_bkg_cut_sim_iter_1_jet60.root", "sphenix_primary_analysis_sim_run28_output"),
        #("hadd -f -k output_dijet_bkg_cut_sim_iter_1.root output_dijet_bkg_cut_sim_iter_1_jet5.root output_dijet_bkg_cut_sim_iter_1_jet12.root output_dijet_bkg_cut_sim_iter_1_jet20.root output_dijet_bkg_cut_sim_iter_1_jet30.root output_dijet_bkg_cut_sim_iter_1_jet40.root output_dijet_bkg_cut_sim_iter_1_jet50.root output_dijet_bkg_cut_sim_iter_1_jet60.root", "sphenix_primary_analysis_sim_run28_output"),
        #("hadd -f -k output_none_bkg_cut_sim_iter_1_jet5.root output_none_bkg_cut_sim_iter_1_jet5_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet12.root output_none_bkg_cut_sim_iter_1_jet12_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet20.root output_none_bkg_cut_sim_iter_1_jet20_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet30.root output_none_bkg_cut_sim_iter_1_jet30_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet40.root output_none_bkg_cut_sim_iter_1_jet40_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet50.root output_none_bkg_cut_sim_iter_1_jet50_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet60.root output_none_bkg_cut_sim_iter_1_jet60_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1.root output_none_bkg_cut_sim_iter_1_jet12.root output_none_bkg_cut_sim_iter_1_jet20.root output_none_bkg_cut_sim_iter_1_jet30.root output_none_bkg_cut_sim_iter_1_jet40.root output_none_bkg_cut_sim_iter_1_jet50.root output_none_bkg_cut_sim_iter_1_jet60.root", "sphenix_primary_analysis_sim_run28_output"),
        #("hadd -f -k output_none_bkg_cut_sim_iter_1.root output_none_bkg_cut_sim_iter_1_jet5.root output_none_bkg_cut_sim_iter_1_jet12.root output_none_bkg_cut_sim_iter_1_jet20.root output_none_bkg_cut_sim_iter_1_jet30.root output_none_bkg_cut_sim_iter_1_jet40.root output_none_bkg_cut_sim_iter_1_jet50.root output_none_bkg_cut_sim_iter_1_jet60.root", "sphenix_primary_analysis_sim_run28_output"),
    ],
    3: [
        ("condor_submit analyze_sim_2_condor.sub", None)
    ],
    4: [
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_2.root output_dijet_bkg_cut_sim_iter_2_jet*_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_2.root output_none_bkg_cut_sim_iter_2_jet*_*.root", "sphenix_primary_analysis_sim_run28_output"),
    ],
    5: [
        ("condor_submit analyze_data_condor.sub", None)
    ],
    6: [
        ("hadd -f -k output_pu_correct_data_dijet_bkg_cut.root output_pu_correct_dijet_bkg_cut_*.root", "analysis_data_run28_output"),
        ("hadd -f -k output_pu_correct_data_efrac_bkg_cut.root output_pu_correct_efrac_bkg_cut_*.root", "analysis_data_run28_output"),
    ],
    7: [
        ("root -l -q -b 'get_reweightmatrix.C(\"sphenix_primary_analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_2.root\",\"analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root\",\"sphenix_primary_run28_output_files/output_reweighted_respmatrix_run28_iter_2_dijet_bkg_cut.root\")'", None),
        ("root -l -q -b 'get_reweightmatrix.C(\"sphenix_primary_analysis_sim_run28_output/output_none_bkg_cut_sim_iter_2.root\",\"analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut.root\",\"sphenix_primary_run28_output_files/output_reweighted_respmatrix_run28_iter_2_none_bkg_cut.root\")'", None),  
    ],
    8: [
        ("condor_submit analyze_sim_3_condor.sub", None)
    ],
    9: [
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_3.root output_dijet_bkg_cut_sim_iter_3_jet*_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_3.root output_none_bkg_cut_sim_iter_3_jet*_*.root", "sphenix_primary_analysis_sim_run28_output"),
    ],
    10: [
        ("root -l -q -b 'do_unfold.C(\"sphenix_primary_analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3.root\",\"analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root\",\"sphenix_primary_run28_output_files/output_unfolded_data_8calibetbin_dijet_bkg_cut_\",11)'", None),
        ("root -l -q -b 'do_unfold.C(\"sphenix_primary_analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3.root\",\"analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut.root\",\"sphenix_primary_run28_output_files/output_unfolded_data_8calibetbin_efrac_bkg_cut_\",11)'", None),
        ("hadd -f -k output_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_2sigma_noise_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_4sigma_noise_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_clus_smear_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_half1_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_half2_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_jerdown_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_jerup_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_jesdown_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_jesup_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_ohcal_mc_data_var_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_run28_iter_3_1000toys.root", "sphenix_primary_run28_output_files"),   
        ("hadd -f -k output_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_2sigma_noise_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_4sigma_noise_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_clus_smear_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_half1_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_half2_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_jerdown_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_jerup_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_jesdown_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_jesup_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_ohcal_mc_data_var_run28_iter_3_1000toys.root output_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_run28_iter_3_1000toys.root", "sphenix_primary_run28_output_files"),   
    ],
    11: [
        ("python3 apply_mbd_vertex_efficiency.py ../mbd_vertex_eff/sim_based_mbd_vertex_eff_w_syst.root sphenix_primary_run28_output_files/output_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root", None),
        ("python3 apply_mbd_vertex_efficiency.py ../mbd_vertex_eff/sim_based_mbd_vertex_eff_w_syst.root sphenix_primary_run28_output_files/output_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root", None),
        ("python3 mbd_efficiency_effects.py sphenix_primary_run28_output_files/output_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root sphenix_primary_mbd_efficiency_effects_dijet_bkg_cut.root", None),
        ("python3 mbd_efficiency_effects.py sphenix_primary_run28_output_files/output_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root sphenix_primary_mbd_efficiency_effects_efrac_bkg_cut.root", None),
    ],
    12 : [
        ("root -l -q -b 'plot_result.C(0)'", None),
        ("root -l -q -b 'plot_result.C(1)'", None),
    ],
    13: [
        ('root -l -q -b \'do_closure.C("sphenix_primary_analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3.root","sphenix_primary_run28_output_files/output_closure_dijet_bkg_cut_sim_run28_iter_3_1000toys.root")\'', None),
        ('root -l -q -b \'do_closure.C("sphenix_primary_analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3.root","sphenix_primary_run28_output_files/output_closure_none_bkg_cut_sim_run28_iter_3_1000toys.root")\'', None),
        ('root -l -q -b \'plot_closure.C("sphenix_primary_run28_output_files/output_closure_dijet_bkg_cut_sim_run28_iter_3_1000toys.root")\'', None),   
        ('root -l -q -b \'plot_closure.C("sphenix_primary_run28_output_files/output_closure_none_bkg_cut_sim_run28_iter_3_1000toys.root")\'', None), 
    ],
    14: [
        ("root -l -q -b 'plot_unfold.C()'", None),
        ("root -l -q -b 'plot_reweight.C()'", None),
    ]

}


def run_step(step_number: int):
    if step_number not in steps:
        print(f"❌ Step {step_number} is not defined.")
        return

    print(f"➡️ Running Step {step_number}...")
    for cmd, cwd in steps[step_number]:
        print(f"   $ {cmd}   (in {cwd or 'current directory'})")
        subprocess.run(cmd, shell=True, check=True, cwd=cwd)
    print(f"✅ Finished Step {step_number}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_steps.py <step_number>")
        sys.exit(1)

    try:
        step = int(sys.argv[1])
        run_step(step)
    except ValueError:
        print("❌ Step number must be an integer.")
