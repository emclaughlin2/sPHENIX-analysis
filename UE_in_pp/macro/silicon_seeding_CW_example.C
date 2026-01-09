#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Trkr_Clustering.C>
#include <Trkr_Reco.C> // note : here fore the trackreco/*.h
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>
#include <intt/InttOdbcQuery.h>

#include <globalvertex/GlobalVertexReco.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <tpccalib/PHTpcResiduals.h>

#include <trackingqa/SiliconSeedsQA.h>
#include <trackingqa/TpcSeedsQA.h>
#include <trackingqa/TpcSiliconQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/InttClusterQA.h>
#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>
#include <calotrigger/TriggerRunInfoReco.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
//R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
//R__LOAD_LIBRARY(libtpc.so)
//R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libtrackingqa.so)
//R__LOAD_LIBRARY(libmvtxrawhitqa.so)
//R__LOAD_LIBRARY(libinttrawhitqa.so)

//R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libglobalvertex.so)
//R__LOAD_LIBRARY(libSeedChecker.so)
//R__LOAD_LIBRARY(libffamodules.so)
//R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libffarawmodules.so)


#include <siseeding_check/SiSeeding_check.h>
R__LOAD_LIBRARY(libSiSeeding_check.so)
R__LOAD_LIBRARY(libtrack_reco.so)

int silicon_seeding_test(
    int process_id = 0,
    int run_num = 53018,
    int nevents = 20,
    string output_directory = "/sphenix/tg/tg01/commissioning/INTT/work/cwshih/sPH_silicon_seeding/test/run53018",
  
    // note : for macro
    string input_directory = "/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana504_2024p023_v001/DST_TRKR_CLUSTER/run_00053000_00053100/dst",
    string input_filename = "DST_TRKR_CLUSTER_run2pp_ana504_2024p023_v001-00053018-00000.root",

    string output_file_name_suffix_in = "MVTXStrobe",

    bool UseGL1_in = true,
    bool ShowEvtDisplay_in = false,
    bool SearchInIntt_in = false,
    bool IsCheckTiming_in = false,
    bool ShowEvtDisplay2_in = false,
    bool checkLargeDeltaphi_evtDisplay_in = false,
    bool checkClusterDuplicate_in = true,

    std::pair<double,double> avg_vertexXY = {-0.009827, 0.1075},
    int ActsSeeding_Verbosity = 0,
    bool convertSeeds = true
){
    TStopwatch* watch = new TStopwatch();
    watch->Start();
  
    Fun4AllServer* se = Fun4AllServer::instance();
    
    Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_track");
    in->Verbosity(2);
    in->AddFile( input_directory + "/" + input_filename); 
    // in->AddListFile( input_directory + "/" + input_filename );
    se->registerInputManager(in);

    auto rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", run_num);
    Enable::CDB = true;
    rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
    //  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
    rc->set_uint64Flag("TIMESTAMP", run_num);

    // std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
    std::string geofile = "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/Tracking_Geometry/93/42/9342254d8393560d92364f4c048fcf03_TrackingGeometry-ana500.root";


    SyncReco *sync = new SyncReco();
    se->registerSubsystem(sync);
    
    HeadReco *head = new HeadReco();
    se->registerSubsystem(head);
    
    FlagHandler *flag = new FlagHandler();
    se->registerSubsystem(flag);

    G4TRACKING::SC_CALIBMODE = false;
    TRACKING::pp_mode = true;
    
    Enable::MVTX_APPLYMISALIGNMENT = true;
    ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

    Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
    ingeo->AddFile(geofile);
    se->registerInputManager(ingeo);
    
    //  auto eventhit = new Fun4AllDstInputManager("EventInputManager");
    //  eventhit->AddListFile(inputeventFile);
    //se->registerInputManager(eventhit);
    

    // auto hitsinclus = new Fun4AllDstInputManager("ClusterInputManager");
    // hitsinclus->AddListFile(inputclusterRawHitFile);
    // se->registerInputManager(hitsinclus);
    cout<<"i am here"<<endl;

    G4MAGNET::magfield_rescale = 1;

    TrackingInit();


    auto silicon_Seeding = new PHActsSiliconSeeding;
    //silicon_Seeding->Verbosity(1);
    silicon_Seeding->Verbosity(ActsSeeding_Verbosity);
    silicon_Seeding->setStrobeRange(-5,5);
    silicon_Seeding->setinttRPhiSearchWindow(0.2);
    silicon_Seeding->setinttZSearchWindow(1.0);  
    silicon_Seeding->seedAnalysis(false);

    if (SearchInIntt_in) {silicon_Seeding->searchInIntt();}
    
    //  silicon_Seeding->seedAnalysis(true);
    se->registerSubsystem(silicon_Seeding);

    cout<<"i am here"<<endl;

    auto merger = new PHSiliconSeedMerger;
    merger->Verbosity(0);
    se->registerSubsystem(merger);

    G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
    std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;
    if (G4TRACKING::convert_seeds_to_svtxtracks)
    {
        cout<<"i am here"<<endl;
        auto converter = new TrackSeedTrackMapConverter("SiliconSeedConverter");
        converter->setTrackSeedName("SiliconTrackSeedContainer");
        converter->setTrackMapName("SiliconSvtxTrackMap");
        converter->setFieldMap(G4MAGNET::magfield_tracking);
        converter->Verbosity(0);
        se->registerSubsystem(converter);
    }

    cout<<"i am here after merger"<<endl;
    
    auto finder = new PHSimpleVertexFinder;
    finder->Verbosity(0);
    //  finder->setDcaCut(0.075);
    finder->setDcaCut(0.05);
    finder->setTrackPtCut(0.1);
    finder->setBeamLineCut(1);
    finder->setTrackQualityCut(300); 
    finder->setNmvtxRequired(3);
    //  finder->setOutlierPairCut(0.075);
    finder->setOutlierPairCut(0.1);

    finder->setTrackMapName("SiliconSvtxTrackMap");
    finder->setVertexMapName("SiliconSvtxVertexMap");
    se->registerSubsystem(finder);

    cout<<"i am here after finder"<<endl;



    SiSeeding_check * IBR = new SiSeeding_check(
        "SiSeeding_check",
        process_id,
        run_num,
        output_directory,

        output_file_name_suffix_in,

        UseGL1_in,
        ShowEvtDisplay_in,
        SearchInIntt_in,
        IsCheckTiming_in,
        ShowEvtDisplay2_in,
        checkLargeDeltaphi_evtDisplay_in,
        checkClusterDuplicate_in
    );
    string final_output_file_name = IBR->GetOutputFileName();
    IBR -> SetAvgVtxXY(avg_vertexXY.first, avg_vertexXY.second);

    cout<<"final_output_file_name: "<<final_output_file_name<<endl;

  
    system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));  


    se->registerSubsystem(IBR);

    se->run(nevents);
    se->End();

    watch->Stop();
    cout << "----------------------------------" << endl;
    cout << "Real time: " << watch->RealTime() << endl;
    cout << "CPU time:  " << watch->CpuTime() << endl;  

    system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));

    std::cout<<"Done!"<<std::endl;

    delete se;

    return 888;
}