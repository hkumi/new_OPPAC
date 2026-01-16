#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "SensorHit.hh"
#include "G4AnalysisManager.hh"

#include <vector>

class MySensitiveDetector : public G4VSensitiveDetector {
public:
    MySensitiveDetector(G4String name);
    ~MySensitiveDetector() override;

    void Initialize(G4HCofThisEvent* HCE) override;
    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) override;
    void EndOfEvent(G4HCofThisEvent* HCE) override;

private:
    SensorHitsCollection* SensorCollection = nullptr;
    G4int HitID = -1;

    void SaveToCSV(G4ThreeVector posPhotons, G4int copyNo);
    void RecordSensorData(const std::string& volumeName, int ntupleIndex, 
                         double posX, double posY, int event, int copyNo, 
                         G4AnalysisManager* man);
    
    double calculateWeightedMeanX(double P_x1, double N_x1, double sigma_x1, 
                                  double P_x2, double N_x2, double sigma_x2);
    double calculateWeightedMeanY(double P_y1, double N_y1, double sigma_y1, 
                                  double P_y2, double N_y2, double sigma_y2);
    
    // UPDATED: Added eventID parameter
    void FitHistogram(const std::vector<double>& xPosition1, 
                      const std::vector<double>& xPosition2, 
                      const std::vector<double>& yPosition1, 
                      const std::vector<double>& yPosition2,
                      G4int eventID);
};

#endif
