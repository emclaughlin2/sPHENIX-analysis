#!/usr/bin/env python3
import sys
import subprocess

# Steps definition with optional working directories
steps = {
    1: [
        ("condor_submit analyze_sim_1_condor.sub", None)
    ],
    2: [
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet10.root output_dijet_bkg_cut_sim_iter_1_jet10_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet15.root output_dijet_bkg_cut_sim_iter_1_jet15_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet20.root output_dijet_bkg_cut_sim_iter_1_jet20_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet30.root output_dijet_bkg_cut_sim_iter_1_jet30_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet50.root output_dijet_bkg_cut_sim_iter_1_jet50_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1_jet70.root output_dijet_bkg_cut_sim_iter_1_jet70_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_1.root output_dijet_bkg_cut_sim_iter_1_jet10.root output_dijet_bkg_cut_sim_iter_1_jet15.root output_dijet_bkg_cut_sim_iter_1_jet20.root output_dijet_bkg_cut_sim_iter_1_jet30.root output_dijet_bkg_cut_sim_iter_1_jet50.root output_dijet_bkg_cut_sim_iter_1_jet70.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet10.root output_none_bkg_cut_sim_iter_1_jet10_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet15.root output_none_bkg_cut_sim_iter_1_jet15_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet20.root output_none_bkg_cut_sim_iter_1_jet20_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet30.root output_none_bkg_cut_sim_iter_1_jet30_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet50.root output_none_bkg_cut_sim_iter_1_jet50_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1_jet70.root output_none_bkg_cut_sim_iter_1_jet70_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_1.root output_none_bkg_cut_sim_iter_1_jet10.root output_none_bkg_cut_sim_iter_1_jet15.root output_none_bkg_cut_sim_iter_1_jet20.root output_none_bkg_cut_sim_iter_1_jet30.root output_none_bkg_cut_sim_iter_1_jet50.root output_none_bkg_cut_sim_iter_1_jet70.root", "analysis_sim_run28_output"),
    ],
    3: [
        ("condor_submit analyze_sim_2_condor.sub", None)
    ],
    4: [
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_2.root output_dijet_bkg_cut_sim_iter_2_jet*_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_2.root output_none_bkg_cut_sim_iter_2_jet*_*.root", "analysis_sim_run28_output"),
    ],
    5: [
        ("condor_submit analyze_data_condor.sub", None)
    ],
    6: [
        ("hadd -f -k output_pu_correct_data_dijet_bkg_cut.root output_pu_correct_dijet_bkg_cut_*.root", "analysis_data_run28_output"),
        ("hadd -f -k output_pu_correct_data_efrac_bkg_cut.root output_pu_correct_efrac_bkg_cut_*.root", "analysis_data_run28_output"),
    ],
    7: [
        ("root -l -q -b 'get_reweightmatrix.C(\"analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_2.root\",\"analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root\",\"run28_output_files/output_reweighted_respmatrix_run28_iter_2_dijet_bkg_cut.root\")'", None),
        ("root -l -q -b 'get_reweightmatrix.C(\"analysis_sim_run28_output/output_none_bkg_cut_sim_iter_2.root\",\"analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut.root\",\"run28_output_files/output_reweighted_respmatrix_run28_iter_2_none_bkg_cut.root\")'", None),  
    ],
    8: [
        ("condor_submit analyze_sim_3_condor.sub", None)
    ],
    9: [
        ("hadd -f -k output_dijet_bkg_cut_sim_iter_3.root output_dijet_bkg_cut_sim_iter_3_jet*_*.root", "analysis_sim_run28_output"),
        ("hadd -f -k output_none_bkg_cut_sim_iter_3.root output_none_bkg_cut_sim_iter_3_jet*_*.root", "analysis_sim_run28_output"),
        ],
    10: [
        ("root -l -q -b 'do_unfold.C(\"analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3.root\",\"analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root\",\"run28_output_files/output_unfolded_data_dijet_bkg_cut_\",1)'", None),
        ("root -l -q -b 'do_unfold.C(\"analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3.root\",\"analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut.root\",\"run28_output_files/output_unfolded_data_efrac_bkg_cut_\",1)'", None),
        #("hadd -f -k output_unfolded_data_dijet_bkg_cut_run28_iter_3_1000toys.root output_unfolded_data_dijet_bkg_cut_calib_dijet_*run28_iter_3_1000toys.root", "run28_output_files"),
        #("hadd -f -k output_unfolded_data_efrac_bkg_cut_run28_iter_3_1000toys.root output_unfolded_data_efrac_bkg_cut_calib_dijet_*run28_iter_3_1000toys.root", "run28_output_files"),

    ],
    11: [
        ("root -l -q -b 'plot_unfold.C()'", None),
        ("root -l -q -b 'plot_result.C()'", None),
        ("root -l -q -b 'plot_reweight.C()'", None),
    ],
    12: [
        ('root -l -q -b \'do_closure.C("analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3.root","run28_output_files/output_closure_dijet_bkg_cut_sim_run28_iter_3_1000toys.root")\'', None),
        ('root -l -q -b \'do_closure.C("analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3.root","run28_output_files/output_closure_none_bkg_cut_sim_run28_iter_3_1000toys.root")\'', None),
        ('root -l -q -b \'plot_closure.C("run28_output_files/output_closure_dijet_bkg_cut_sim_run28_iter_3_1000toys.root")\'', None),   
        ('root -l -q -b \'plot_closure.C("run28_output_files/output_closure_none_bkg_cut_sim_run28_iter_3_1000toys.root")\'', None), 
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
