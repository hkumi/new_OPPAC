#include "StepAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
StepAction::StepAction() {}

StepAction::~StepAction() {}

void StepAction::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    
    // Only process optical photons
    if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return;
    
    // Kill photons that have taken too many steps (trapped in loops)
    if (track->GetCurrentStepNumber() > MAX_PHOTON_STEPS) {
        track->SetTrackStatus(fStopAndKill);
        
        // Optional debug output (only print occasionally)
        static G4int killedCount = 0;
        killedCount++;
        if (killedCount % 100 == 0) {
            G4cout << "Killed " << killedCount << " stuck optical photons (too many steps)" << G4endl;
        }
        return;
    }
    
    // Also kill photons with extremely small energy (optional)
    if (track->GetKineticEnergy() < 0.1 * eV) {
        track->SetTrackStatus(fStopAndKill);
        return;
    }
}
