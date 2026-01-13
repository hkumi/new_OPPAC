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
    

  man->CreateNtuple("LeftData1", "LeftData1");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleIColumn("copyNo");
  man->FinishNtuple(0);

  man->CreateNtuple("RightData2", "RightData2");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleIColumn("copyNo");
  man->FinishNtuple(1);

  man->CreateNtuple("BottomData3", "BottomData3");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleIColumn("copyNo");
  man->FinishNtuple(2);

  man->CreateNtuple("TopData4", "TopData4");
  man->CreateNtupleDColumn("x");
  man->CreateNtupleDColumn("y");
  man->CreateNtupleDColumn("z");
  man->CreateNtupleIColumn("copyNo");
  man->FinishNtuple(3);

    
 
  man->CreateNtuple("x_reconstruction","x_reconstruction");
  man->CreateNtupleDColumn("x");
  man->FinishNtuple(4);

  man->CreateNtuple("y_reconstruction","y_reconstruction");
  man->CreateNtupleDColumn("y");
  man->FinishNtuple(5);

  
 
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
