#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "SensorHit.hh"
#include <map>
#include <vector>
#include <TH1D.h>
#include "TF1.h"
#include "TFile.h"

class MySensitiveDetector : public G4VSensitiveDetector {
public:
    MySensitiveDetector(G4String name);
    ~MySensitiveDetector();
    
    virtual void Initialize(G4HCofThisEvent* HCE) override;
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* R0hist) override;
    virtual void EndOfEvent(G4HCofThisEvent* HCE) override;
    
    void SaveToCSV(G4ThreeVector posPhotons, G4int copyNo);
    G4double GetSensorPositionFromCopyNo(G4int copyNo, G4bool getX);
    void RecordSensorData(G4String volumeName, int ntupleIndex, 
                         double posX, double posY, int event, int copyNo);
    
    // Updated function signatures with position vectors instead of hit counts
    std::pair<double, double> CalculateCenterOfGravity(
        const std::map<G4int, std::vector<G4double>>& rightXPos,
        const std::map<G4int, std::vector<G4double>>& rightYPos,
        const std::map<G4int, std::vector<G4double>>& leftXPos,
        const std::map<G4int, std::vector<G4double>>& leftYPos,
        const std::map<G4int, std::vector<G4double>>& topXPos,
        const std::map<G4int, std::vector<G4double>>& topYPos,
        const std::map<G4int, std::vector<G4double>>& bottomXPos,
        const std::map<G4int, std::vector<G4double>>& bottomYPos);
    
    std::pair<double, double> PerformSensorGaussianFit(
        const std::map<G4int, std::vector<G4double>>& rightXPos,
        const std::map<G4int, std::vector<G4double>>& rightYPos,
        const std::map<G4int, std::vector<G4double>>& leftXPos,
        const std::map<G4int, std::vector<G4double>>& leftYPos,
        const std::map<G4int, std::vector<G4double>>& topXPos,
        const std::map<G4int, std::vector<G4double>>& topYPos,
        const std::map<G4int, std::vector<G4double>>& bottomXPos,
        const std::map<G4int, std::vector<G4double>>& bottomYPos,
        G4int eventID);
    
    // Weighted mean calculation methods from second code
    double calculateWeightedMeanX(double P_x1, double N_x1, double sigma_x1, 
                                  double P_x2, double N_x2, double sigma_x2);
    
    double calculateWeightedMeanY(double P_y1, double N_y1, double sigma_y1, 
                                  double P_y2, double N_y2, double sigma_y2);
    
private:
    SensorHitsCollection* SensorCollection;
    
   
    // Static members for analysis file
    static TFile* analysisFile;
    static bool analysisFileCreated;
    static std::mutex analysisFileMutex;  // Add mutex for thread safety
    
    // Constants (should match your geometry)
    static constexpr G4double SENSOR_X_POS = 28.5;  // mm
    static constexpr G4double SENSOR_Y_POS = 28.5;  // mm
    static constexpr G4double CELL_SIZE = 1.84;     // mm
};

#endif // DETECTOR_HH
