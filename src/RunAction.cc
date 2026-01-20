#include "RunAction.hh"
#include "Run.hh"
#include "G4Run.hh" 
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include <assert.h>

RunAction::RunAction()
{
  //theSDName.push_back(G4String("IonPro"));
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  man->SetNtupleMerging(true);
  man->SetVerboseLevel( 1 );
    

  man->CreateNtuple("RightData1", "Right sensors (sensor_Vol1)");  
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleIColumn("event");
  man->CreateNtupleIColumn("copyNo");
  man->FinishNtuple(0);

  man->CreateNtuple("LeftData2", "Left sensors (sensor_Vol2)");  
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleIColumn("event");
  man->CreateNtupleIColumn("copyNo");
man->FinishNtuple(1);

  man->CreateNtuple("TopData3", "TopData3");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleIColumn("event");
  man->CreateNtupleIColumn("copyNo");
  man->FinishNtuple(2);

  man->CreateNtuple("BottomData4", "BottomData4");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleIColumn("event");
  man->CreateNtupleIColumn("copyNo");
  man->FinishNtuple(3);

    
  man->CreateH1("CopyNo_Right", "Number of photons detected in right sensor", 25, 0, 25);    
  man->CreateH1("CopyNo_Top", "Number of photons detected in top sensor", 25, 0, 25);       
  man->CreateH1("CopyNo_Left", "Number of photons detected in left sensor ", 25, 0, 25);     
  man->CreateH1("CopyNo_Bottom", "Number of photons detected in bottom sensor", 25, 0, 25); 
  
  man->CreateNtuple("xy_reconstruction", "Position Reconstruction");
man->CreateNtupleDColumn("x");  // Column 0
man->CreateNtupleDColumn("y");  // Column 1
man->FinishNtuple(4);  // This becomes ntuple 4

man->CreateNtuple("cog_reconstruction", "Center of Gravity");
man->CreateNtupleDColumn("x_cog");  // Column 0
man->CreateNtupleDColumn("y_cog");  // Column 1  
man->FinishNtuple(5);  // This becomes ntuple 5
  
 
  man->CreateH2("xy2 ","xy2", 100, -30, 30, 100, -30, 30);

/*
  man->CreateNtuple("yposition","yposition");
  man->CreateNtupleDColumn("y");
  man->FinishNtuple(1);


  man->CreateNtuple("SiPM","siPM");
  man->CreateNtupleDColumn("Number");
  man->FinishNtuple(2);

  man->CreateNtuple("time","time");
  man->CreateNtupleDColumn("time");
  man->FinishNtuple(3);

  man->CreateNtuple("Edep","Edep");
  man->CreateNtupleDColumn("Edep");
  man->FinishNtuple(4);


  man->CreateH2("xy1 ","xy1", 100, -60, 60, 100, -60, 60);*/




}

RunAction::~RunAction()
{
  //theSDName.clear();
}
/*
G4Run* RunAction::GenerateRun()
{ 
 // return new Run(theSDName);
}*/

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
 G4cout << "*** Run " << aRun->GetRunID() << " start." << G4endl;

 G4AnalysisManager *man = G4AnalysisManager::Instance();

 G4int runID = aRun->GetRunID();

 std::stringstream strRunID;
 strRunID << runID;

 
 G4String filename = "output" + strRunID.str() + ".root";
 man->OpenFile(filename);
 G4cout << "Output file opened: " << filename <<G4endl; // Debugging line

}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
 // static G4String regName = "IonPro";
  //Run* Run_verbose = (Run*)aRun;

  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->Write();
  man->CloseFile();
  G4cout << "Output file closed." << G4endl; // Debugging line
}
