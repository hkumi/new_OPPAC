#ifndef STEPACTION_HH
#define STEPACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

class StepAction : public G4UserSteppingAction
{
public:
    StepAction();  
    virtual ~StepAction();
    
    virtual void UserSteppingAction(const G4Step* step);
    
private:
    static const G4int MAX_PHOTON_STEPS = 10000;
};

#endif
