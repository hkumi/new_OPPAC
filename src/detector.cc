#include "detector.hh"
#include <numeric>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include "TF1.h"
#include "G4SDManager.hh"
#include "G4String.hh"
#include "G4HCofThisEvent.hh"
#include "SensorHit.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
    G4String HCname = "SensorCollection";
    collectionName.insert(HCname);
}

MySensitiveDetector::~MySensitiveDetector()
{}

void MySensitiveDetector::Initialize(G4HCofThisEvent* HCE) {
    // Create the hits collection
    SensorCollection = new SensorHitsCollection(SensitiveDetectorName, collectionName[0]);
    
    // Register the collection with the event
    static G4int HCID = -1;
    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    
    if (HCE) {
        HCE->AddHitsCollection(HCID, SensorCollection);
    } else {
        G4cerr << "Error: HCE is null in Initialize!" << G4endl;
    }
}

void MySensitiveDetector::SaveToCSV(G4ThreeVector posPhotons, G4int copyNo) {
    std::ofstream outputFile("photon_positions_4.csv", std::ios::app);
    
    if (outputFile.is_open()) {
        outputFile << posPhotons.x()/mm << "," 
                   << posPhotons.y()/mm << "," 
                   << posPhotons.z()/mm << "," 
                   << copyNo << "\n";
    } else {
        G4cerr << "Failed to open file!" << G4endl;
    }
    
    outputFile.close();
}

void MySensitiveDetector::RecordSensorData(const std::string& volumeName, int ntupleIndex, 
                                          double posX, double posY, int event, int copyNo, 
                                          G4AnalysisManager* man) {
    
    if (ntupleIndex >= 0 && ntupleIndex < 4) {
        // Fill ntuple columns
        man->FillNtupleDColumn(ntupleIndex, 0, posX);
        man->FillNtupleDColumn(ntupleIndex, 1, posY);
        man->FillNtupleIColumn(ntupleIndex, 2, event);
        man->FillNtupleIColumn(ntupleIndex, 3, copyNo);
        man->AddNtupleRow(ntupleIndex);
        
        //  MAPPING based on ACTUAL geometry:
        // sensor_Vol1 = RIGHT side (positive X) - copyNo 0-24
        // sensor_Vol3 = TOP side (positive Y) - copyNo 25-49
        // sensor_Vol2 = LEFT side (negative X) - copyNo 50-74
        // sensor_Vol4 = BOTTOM side (negative Y) - copyNo 75-99
        
        if (volumeName == "sensor_Vol1") {  //  RIGHT side
            man->FillH1(0, copyNo);  // Right sensors (0-24)
        }
        else if (volumeName == "sensor_Vol3") {  // TOP side
            man->FillH1(1, copyNo - 25);  //  Top sensors (0-24)
        }
        else if (volumeName == "sensor_Vol2") {  //  LEFT side
            man->FillH1(2, copyNo - 50);  //  Left sensors (0-24)
        }
        else if (volumeName == "sensor_Vol4") {  // BOTTOM side
            man->FillH1(3, copyNo - 75);  // Bottom sensors (0-24)
        }
    }
}

double MySensitiveDetector::calculateWeightedMeanX(double P_x1, double N_x1, double sigma_x1, 
                                                   double P_x2, double N_x2, double sigma_x2) {
    if (sigma_x1 == 0 || sigma_x2 == 0) {
        G4cerr << "Error: Sigma values cannot be zero!" << G4endl;
        return 0.0;
    }
    
    double numerator = (P_x1 * N_x1 / sigma_x1) + (P_x2 * N_x2 / sigma_x2);
    double denominator = (N_x1 / sigma_x1) + (N_x2 / sigma_x2);
    
    return numerator / denominator;
}

double MySensitiveDetector::calculateWeightedMeanY(double P_y1, double N_y1, double sigma_y1, 
                                                   double P_y2, double N_y2, double sigma_y2) {
    if (sigma_y1 == 0 || sigma_y2 == 0) {
        G4cerr << "Error: Sigma values cannot be zero!" << G4endl;
        return 0.0;
    }
    
    double numerator = (P_y1 * N_y1 / sigma_y1) + (P_y2 * N_y2 / sigma_y2);
    double denominator = (N_y1 / sigma_y1) + (N_y2 / sigma_y2);
    
    return numerator / denominator;
}

void MySensitiveDetector::FitHistogram(const std::vector<double>& xTop,      
                                       const std::vector<double>& xBottom,   
                                       const std::vector<double>& yRight,    
                                       const std::vector<double>& yLeft,     
                                       G4int eventID) {
    // FIXED: Use the new parameter names
    if (xTop.empty() || xBottom.empty() || 
        yRight.empty() || yLeft.empty()) {
        G4cout << "Event " << eventID << ": Not enough data for fitting" << G4endl;
        return;
    }
    
    int bins = 50;
    
    // Create ROOT file for histograms if needed
    static TFile* rootFile = nullptr;
    static bool firstCall = true;
    if (firstCall) {
        rootFile = new TFile("sensor_hists.root", "RECREATE");
        firstCall = false;
    }
    
    // Process X positions from TOP and BOTTOM sensors
    double minXTop = *std::min_element(xTop.begin(), xTop.end());
    double maxXTop = *std::max_element(xTop.begin(), xTop.end());
    double minXBottom = *std::min_element(xBottom.begin(), xBottom.end());
    double maxXBottom = *std::max_element(xBottom.begin(), xBottom.end());
    
    // Ensure reasonable histogram ranges
    if (maxXTop - minXTop < 0.1) { maxXTop = minXTop + 1.0; minXTop = minXTop - 1.0; }
    if (maxXBottom - minXBottom < 0.1) { maxXBottom = minXBottom + 1.0; minXBottom = minXBottom - 1.0; }
    
    TH1D* histXTop = new TH1D(Form("histXTop_evt%d", eventID), "X Top Distribution (sensor_Vol3)", bins, minXTop, maxXTop);
    TH1D* histXBottom = new TH1D(Form("histXBottom_evt%d", eventID), "X Bottom Distribution (sensor_Vol4)", bins, minXBottom, maxXBottom);
    
    for (double x : xTop) histXTop->Fill(x);
    for (double x : xBottom) histXBottom->Fill(x);
    
    TF1* gaussFitXTop = new TF1("gaussFitXTop", "gaus", minXTop, maxXTop);
    gaussFitXTop->SetLineColor(kRed);
    histXTop->Fit(gaussFitXTop, "QR0");  // "0" for no graphics
    
    TF1* gaussFitXBottom = new TF1("gaussFitXBottom", "gaus", minXBottom, maxXBottom);
    gaussFitXBottom->SetLineColor(kRed);
    histXBottom->Fit(gaussFitXBottom, "QR0");
    
    double meanXTop = gaussFitXTop->GetParameter(1);
    double sigmaXTop = gaussFitXTop->GetParameter(2);
    double amplitudeXTop = gaussFitXTop->GetParameter(0);
    
    double meanXBottom = gaussFitXBottom->GetParameter(1);
    double sigmaXBottom = gaussFitXBottom->GetParameter(2);
    double amplitudeXBottom = gaussFitXBottom->GetParameter(0);
    
    // Process Y positions from RIGHT and LEFT sensors
    double minYRight = *std::min_element(yRight.begin(), yRight.end());
    double maxYRight = *std::max_element(yRight.begin(), yRight.end());
    double minYLeft = *std::min_element(yLeft.begin(), yLeft.end());
    double maxYLeft = *std::max_element(yLeft.begin(), yLeft.end());
    
    // Ensure reasonable histogram ranges
    if (maxYRight - minYRight < 0.1) { maxYRight = minYRight + 1.0; minYRight = minYRight - 1.0; }
    if (maxYLeft - minYLeft < 0.1) { maxYLeft = minYLeft + 1.0; minYLeft = minYLeft - 1.0; }
    
    TH1D* histYRight = new TH1D(Form("histYRight_evt%d", eventID), "Y Right Distribution (sensor_Vol1)", bins, minYRight, maxYRight);
    TH1D* histYLeft = new TH1D(Form("histYLeft_evt%d", eventID), "Y Left Distribution (sensor_Vol2)", bins, minYLeft, maxYLeft);
    
    for (double y : yRight) histYRight->Fill(y);
    for (double y : yLeft) histYLeft->Fill(y);
    
    TF1* gaussFitYRight = new TF1("gaussFitYRight", "gaus", minYRight, maxYRight);
    gaussFitYRight->SetLineColor(kBlue);
    histYRight->Fit(gaussFitYRight, "QR0");
    
    TF1* gaussFitYLeft = new TF1("gaussFitYLeft", "gaus", minYLeft, maxYLeft);
    gaussFitYLeft->SetLineColor(kBlue);
    histYLeft->Fit(gaussFitYLeft, "QR0");
    
    double meanYRight = gaussFitYRight->GetParameter(1);
    double sigmaYRight = gaussFitYRight->GetParameter(2);
    double amplitudeYRight = gaussFitYRight->GetParameter(0);
    
    double meanYLeft = gaussFitYLeft->GetParameter(1);
    double sigmaYLeft = gaussFitYLeft->GetParameter(2);
    double amplitudeYLeft = gaussFitYLeft->GetParameter(0);
    
    // Compute weighted means - FIXED parameter names
    double x_pos = calculateWeightedMeanX(meanXTop, amplitudeXTop, sigmaXTop, meanXBottom, amplitudeXBottom, sigmaXBottom);
    double y_pos = calculateWeightedMeanY(meanYRight, amplitudeYRight, sigmaYRight, meanYLeft, amplitudeYLeft, sigmaYLeft);
    
    // Store the results in ntuple
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    
    // Fill 2D histogram with reconstructed position
    man->FillH2(0, x_pos, y_pos);
    
    // Store X and Y positions in separate ntuples
    man->FillNtupleDColumn(4, 0, x_pos);
    man->AddNtupleRow(4);
    
    man->FillNtupleDColumn(5, 0, y_pos);
    man->AddNtupleRow(5);
    
    // Write histograms to file
    rootFile->cd();
    histXTop->Write();
    histXBottom->Write();
    histYRight->Write();
    histYLeft->Write();
    
    // Cleanup
    delete histXTop; delete histXBottom;
    delete gaussFitXTop; delete gaussFitXBottom;
    delete histYRight; delete histYLeft;
    delete gaussFitYRight; delete gaussFitYLeft;
}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory*) {
    G4Track *track = aStep->GetTrack();
    G4ParticleDefinition* particle = track->GetDefinition();
    
    // Only process optical photons
    if (particle->GetParticleName() != "opticalphoton") {
        return false;
    }
    
    // Kill the photon after detection (optional)
    // track->SetTrackStatus(fStopAndKill);
    
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    
    G4ThreeVector posPhotons = postStepPoint->GetPosition();
    G4int copyNo = touchable->GetCopyNumber();
    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4ThreeVector posDetector = physVol->GetTranslation();
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4String volumeName = touchable->GetVolume()->GetName();
    // DEBUG: Print position and copy number
    G4cout << "Photon hit: volume=" << volumeName 
           << ", copyNo=" << copyNo
           << ", position=(" << posPhotons.x()/mm << ", " 
           << posPhotons.y()/mm << ", " << posPhotons.z()/mm << ") mm" 
           << G4endl;

    // Create a new hit
    SensorHit* aSensorHit = new SensorHit();
    aSensorHit->SetSensorPosition(posPhotons);
    aSensorHit->SetSensorEnergy(track->GetKineticEnergy());
    aSensorHit->SetVolumeName(volumeName);
    aSensorHit->SetSensorNumber(copyNo);  
    
    // Add hit to collection
    if (SensorCollection) {
        SensorCollection->insert(aSensorHit);
    } else {
        G4cerr << "SensorCollection not initialized!" << G4endl;
        delete aSensorHit;
        return false;
    }
    
    // Record data in ntuples based on which sensor was hit
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    
    if (volumeName == "sensor_Vol1") {
    // Right side sensors - measure Y position
    // Store: X = sensor center X (+23 mm), Y = photon Y position
    RecordSensorData(volumeName, 0, posDetector[0]/mm, posPhotons.y()/mm, evt, copyNo, man);
    SaveToCSV(posPhotons, copyNo);
    } 
    else if (volumeName == "sensor_Vol2") {
            // Left side sensors - measure Y position  
           // Store: X = sensor center X (-23 mm), Y = photon Y position
            RecordSensorData(volumeName, 1, posDetector[0]/mm, posPhotons.y()/mm, evt, copyNo, man);
            SaveToCSV(posPhotons, copyNo);
    } 
    else if (volumeName == "sensor_Vol3") {
    // Top side sensors - measure X position
    // Store: X = photon X position, Y = sensor center Y (23 mm)
            RecordSensorData(volumeName, 2, posPhotons.x()/mm, posDetector[1]/mm, evt, copyNo, man);
            SaveToCSV(posPhotons, copyNo);
    } 
    else if (volumeName == "sensor_Vol4") {
    // Bottom side sensors - measure X position
    // Store: X = photon X position, Y = sensor center Y (-23 mm)
            RecordSensorData(volumeName, 3, posPhotons.x()/mm, posDetector[1]/mm, evt, copyNo, man);
            SaveToCSV(posPhotons, copyNo);
    }
    else {
        G4cout << "Warning: Photon detected in unknown volume: " << volumeName << G4endl;
    }
    
    return true;
}

void MySensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {
    if (!SensorCollection) return;
    
    G4int nHits = SensorCollection->entries();
    G4int evtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    //G4cout << "Event " << evtID << ": Number of photon hits = " << nHits << G4endl;
    
    if (nHits > 0) {
        // CORRECTED: Update variable names to match actual geometry
        // yPositions1 = RIGHT sensors (sensor_Vol1) - measure Y at fixed X = +23 mm
        // yPositions2 = LEFT sensors (sensor_Vol2) - measure Y at fixed X = -23 mm  
        // xPositions1 = TOP sensors (sensor_Vol3) - measure X at fixed Y = +23 mm
        // xPositions2 = BOTTOM sensors (sensor_Vol4) - measure X at fixed Y = -23 mm
        
        std::vector<double> yPositionsRight, yPositionsLeft; // For RIGHT and LEFT sensors
        std::vector<double> xPositionsTop, xPositionsBottom; // For TOP and BOTTOM sensors
        
        // Process hits
        for (G4int i = 0; i < nHits; i++) {
            SensorHit* hit = (*SensorCollection)[i];
            if (!hit) continue;
            
            G4ThreeVector position = hit->GetSensorPosition();
            G4String volumeName = hit->GetVolumeName();
            
            // Separate data based on sensor type
            if (volumeName == "sensor_Vol1") {
                // RIGHT side - Y position at fixed X = +23 mm
                yPositionsRight.push_back(position.y());
            } 
            else if (volumeName == "sensor_Vol2") {
                // LEFT side - Y position at fixed X = -23 mm
                yPositionsLeft.push_back(position.y());
            } 
            else if (volumeName == "sensor_Vol3") {
                // TOP side - X position at fixed Y = +23 mm
                xPositionsTop.push_back(position.x());
            } 
            else if (volumeName == "sensor_Vol4") {
                // BOTTOM side - X position at fixed Y = -23 mm
                xPositionsBottom.push_back(position.x());
            }
        }
        
        // CORRECTED debug output
        G4cout << "Event " << evtID << " sensor hits distribution: " 
               << "Right(sensor_Vol1)=" << yPositionsRight.size() 
               << ", Left(sensor_Vol2)=" << yPositionsLeft.size()
               << ", Top(sensor_Vol3)=" << xPositionsTop.size()
               << ", Bottom(sensor_Vol4)=" << xPositionsBottom.size() << G4endl;
        
        // Perform Gaussian fits if we have enough data
        G4int minHitsForFit = 1;  // Minimum hits per sensor for fitting
        if (yPositionsRight.size() >= minHitsForFit && 
            yPositionsLeft.size() >= minHitsForFit &&
            xPositionsTop.size() >= minHitsForFit && 
            xPositionsBottom.size() >= minHitsForFit) {
            
            //G4cout << "Performing position reconstruction for event " << evtID << G4endl;
            FitHistogram(xPositionsTop, xPositionsBottom, yPositionsRight, yPositionsLeft, evtID);
            
        } else {
            //G4cout << "Event " << evtID << ": Not enough hits for fitting (need at least " 
            //       << minHitsForFit << " per sensor)" << G4endl;
        }
    }
    
    // Optional: Clear hits for next event (Geant4 usually handles this)
    // SensorCollection->clear();
}
