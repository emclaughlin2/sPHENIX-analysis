#!/usr/bin/env python3
import sys
import subprocess

# Steps definition with optional working directories
steps = {
    1: [
        ("condor_submit analyze_herwig_1_condor.sub", None)
    ],
    2: [
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_herwig_jet12.root output_dijet_bkg_cut_sim_iter_1_herwig_jet12_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_herwig_jet20.root output_dijet_bkg_cut_sim_iter_1_herwig_jet20_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_herwig_jet30.root output_dijet_bkg_cut_sim_iter_1_herwig_jet30_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_herwig_jet40.root output_dijet_bkg_cut_sim_iter_1_herwig_jet40_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_herwig_jet50.root output_dijet_bkg_cut_sim_iter_1_herwig_jet50_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_herwig.root output_dijet_bkg_cut_sim_iter_1_herwig_jet12.root output_dijet_bkg_cut_sim_iter_1_herwig_jet20.root output_dijet_bkg_cut_sim_iter_1_herwig_jet30.root output_dijet_bkg_cut_sim_iter_1_herwig_jet40.root output_dijet_bkg_cut_sim_iter_1_herwig_jet50.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_herwig_jet12.root output_none_bkg_cut_sim_iter_1_herwig_jet12_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_herwig_jet20.root output_none_bkg_cut_sim_iter_1_herwig_jet20_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_herwig_jet30.root output_none_bkg_cut_sim_iter_1_herwig_jet30_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_herwig_jet40.root output_none_bkg_cut_sim_iter_1_herwig_jet40_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_herwig_jet50.root output_none_bkg_cut_sim_iter_1_herwig_jet50_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_herwig.root output_none_bkg_cut_sim_iter_1_herwig_jet12.root output_none_bkg_cut_sim_iter_1_herwig_jet20.root output_none_bkg_cut_sim_iter_1_herwig_jet30.root output_none_bkg_cut_sim_iter_1_herwig_jet40.root output_none_bkg_cut_sim_iter_1_herwig_jet50.root", "sphenix_primary_analysis_sim_run28_output"),
    ],
    3: [
        ("condor_submit analyze_herwig_2_condor.sub", None)
    ],
    4: [
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_2_herwig.root output_dijet_bkg_cut_sim_iter_2_herwig_jet*_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_2_herwig.root output_none_bkg_cut_sim_iter_2_herwig_jet*_*.root", "sphenix_primary_analysis_sim_run28_output"),
    ],
    5: [
        ("condor_submit analyze_data_condor.sub", None)
    ],
    6: [
        ("hadd -f -k output_pu_correct_data_dijet_bkg_cut.root output_pu_correct_dijet_bkg_cut_*.root", "analysis_data_run28_output"),
        ("hadd -f -k output_pu_correct_data_efrac_bkg_cut.root output_pu_correct_efrac_bkg_cut_*.root", "analysis_data_run28_output"),
    ],
    7: [
        ("root -l -q -b 'get_reweightmatrix.C(\"sphenix_primary_analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_2_herwig.root\",\"analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root\",\"sphenix_primary_run28_output_files/output_herwig_reweighted_respmatrix_run28_iter_2_dijet_bkg_cut.root\")'", None),
        ("root -l -q -b 'get_reweightmatrix.C(\"sphenix_primary_analysis_sim_run28_output/output_none_bkg_cut_sim_iter_2_herwig.root\",\"analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut.root\",\"sphenix_primary_run28_output_files/output_herwig_reweighted_respmatrix_run28_iter_2_none_bkg_cut.root\")'", None),  
    ],
    8: [
        ("condor_submit analyze_herwig_3_condor.sub", None)
    ],
    9: [
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_3_herwig.root output_dijet_bkg_cut_sim_iter_3_herwig_jet*_*.root", "sphenix_primary_analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_3_herwig.root output_none_bkg_cut_sim_iter_3_herwig_jet*_*.root", "sphenix_primary_analysis_sim_run28_output"),
    ],
    10: [
        ("root -l -q -b 'do_unfold.C(\"sphenix_primary_analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3_herwig.root\",\"analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root\",\"sphenix_primary_run28_output_files/output_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_\",1)'", None),
        ("root -l -q -b 'do_unfold.C(\"sphenix_primary_analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3_herwig.root\",\"analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut.root\",\"sphenix_primary_run28_output_files/output_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_\",1)'", None),
    ],
    11: [
        ("python3 apply_mbd_vertex_efficiency.py ../mbd_vertex_eff/sim_based_mbd_vertex_eff_w_syst.root sphenix_primary_run28_output_files/output_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_run28_iter_3_1000toys.root sphenix_primary_run28_output_files/output_mbd_correct_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root", None),
        ("python3 apply_mbd_vertex_efficiency.py ../mbd_vertex_eff/sim_based_mbd_vertex_eff_w_syst.root sphenix_primary_run28_output_files/output_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_run28_iter_3_1000toys.root sphenix_primary_run28_output_files/output_mbd_correct_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root", None),
        ("python3 mbd_efficiency_effects.py sphenix_primary_run28_output_files/output_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_run28_iter_3_1000toys.root sphenix_primary_run28_output_files/output_mbd_correct_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root sphenix_primary_mbd_efficiency_effects_herwig_dijet_bkg_cut.root", None),
        ("python3 mbd_efficiency_effects.py sphenix_primary_run28_output_files/output_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_calib_dijet_run28_iter_3_1000toys.root sphenix_primary_run28_output_files/output_mbd_correct_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root sphenix_primary_mbd_efficiency_effects_herwig_efrac_bkg_cut.root", None),
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
