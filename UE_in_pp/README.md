UE\_in\_pp analysis code is separated into 6 categories:
1. src - analysis ttree making code
2. macro - fun4all scripts for generating analysis ttrees
3. analysis - components for performing analysis on ttrees
4. plotting - initial plotting code
5. final\_plotting - final plotting macros used for figures in paper
6. topo\_valid - ttree making code for topocluster validation studies

The src directory should be made with build "ana.536" and all scripts run interactively should also be run with build "ana.536".

In the macro directory:
- Main F4A macro used for analysis is "./Fun4All\_InclusiveJet.C", this macro is run via the "./InclusiveJet.sh" script which is called by the codnor submission scripts
- The condor submission scripts include:
	- datacondor.sub
	- mc\_run28\_jet\#\_condor.sub
	- mc\_run28\_herwig_jet\#\_condor.sub
- List files that are read by the F4A macro can be created for MC via the scripts available in the "./mc\_dst\_files" folder

Already produced ttrees from the preliminary analysis can be found at "/sphenix/tg/tg01/jets/egm2153/UEinppOutput" and have the following general file format:
- Data: output_ana509\_v2\_{runnumber}\_{segment}.root
- Pythia: sim\_run28\_jet{trigger number}\_3sigma\_output\_{segment}.root
- Herwig: sim\_run28\_herwig\_jet{trigger number}\_3sigma\_output\_{segment}.root
This is the expected data format for the input ttree files to all of the analysis scripts 

In the analysis directory the main components (directories) are:
1. jet\_background\_eff:
	- perform jet background cut efficiency calculation for timing cuts (dijet/efrac cuts accounted for in unfolding), main output is efficiency factors included in "../unfolding/analysis\_data.C" lines 207-216
2. jet\_trig\_eff:
	- perform jet trigger efficiency calcuation for 3 different run ranges, 0mrad early, 0mrad mid and 1.5mrad running, main output is fit files "./0mrad\_early\_tefficiency\_fits.root", "./0mrad\_mid\_tefficiency\_fits.root" and "./1.5mrad\_only\_tefficiency\_fits.root" used in "../unfolding/analysis\_data.C" lines 180-198
3. mbd\_vertex\_eff:
	- perform simulation-based mbd vertex efficiency calcuation and data-sim mbd vertex efficiency systematic, main output is "./sim\_based\_mbd\_vertex\_eff\_w\_syst.root" which is used to correct unfolded results in step 11 of "../unfolding/run\_analysis.py"
4. pileup:
	- perform pileup analysis, main output is the pileup correction factors included in "../unfolding/analysis\_data.C" line 26
5. unfolding:
	- use the "./run\_analysis.py" script to run through the analysis, runs "./analysis\_sim.C" and "./analysis\_data.C" scripts, outputs are the pythia truth results and the unfolded data results using the pythia response matrices
	- use the "./run\_herwig\_analysis.py" script to run the full analysis using the herwig samples, outputs are the herwig truth results and the unfolded data results using the herwig response matrices
Unless there are big changes to the analysis, only the scripts in the unfolding directory should need to be rerun to get different results. The results from all other directories should be fairly static.

Running "run\_analysis.py" from within "analysis/unfolding":
- Each step should be run in order as "python3 run\_analysis.py \#"
- Steps:
	1. Run over pythia simulation data with "analysis\_sim.C" to create truth hists and response matrices (iter\_1)
	2. hadd simulation root files (iter\_1)
	3. Run over pythia simulation data again with "analysis\_sim.C" to create truth hists and response matrices with sample trimming applied (iter\_2)
	4. hadd simulation root files (iter\_2)
	5. Run over data with "analysis\_data.C" to create measured hists 
	6. hadd data root files 
	7. Get prior reweighting histograms using data root files and simulation (iter\_2) root files
	8. Run over pythia simulation data again with "analysis\_sim.C" to create truth hists and response matries with sample trimming and prior reweighting applied (iter\_3)
	9. hadd simulation root files (iter\_3)
	10. Do unfolding procedure for nominal and syst variation histograms
	11. Apply mbd vertex efficiency to unfolded results 
	12. Plot results
	13. OPTIONAL Perform closure tests
	14. OPTIONAL Perform iteration optmization
- Herwig samples are run via "python3 run\_herwig\_analysis.py \#"
- Note it makes most sense to run both "run\_analysis.py" and "run\_herwig\_analysis.py" steps 1-11 and then run "run\_analysis.py" step 12 since the plot result step uses the results from both the pythia and herwig analysis. Additionally, the output root files from step 11 for both the run_analysis and run_herwig_analysis are the files used for the result figure macro in the final_plotting directory just copied over to there.
