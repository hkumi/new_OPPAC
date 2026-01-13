#ifndef StepAction_h
#define StepAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>
#include "DC.hh"  // Include the DetectorConstruction header
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "event.hh"

class StepAction : public G4UserSteppingAction
{

public:
  //StepAction(G4String data_file, double lunghezza_collimatore);
  StepAction(MyEventAction* eventAction);    // Constructor with MyEventAction parameter
  
  ~StepAction();

  virtual void UserSteppingAction(const G4Step*);  // Implemented by Geant4 for stepping
  double hysto(double valore);                     // Custom method

private:
  G4String datai_file;
  int s_x1[33];
  int s_x2[33];
  int s_y1[33];
  int s_y2[33];
  double collimatore;
  MyEventAction* fEventAction;  // Pointer to the event action
};

#endif
