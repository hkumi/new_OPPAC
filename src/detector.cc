#include "detector.hh"
#include <numeric>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <TH1D.h>
#include <TTree.h>
#include "TF1.h"
#include "TFile.h"
#include "G4SDManager.hh"
#include "G4String.hh"
#include "G4HCofThisEvent.hh"
#include "SensorHit.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"

// Initialize static members 
TFile* MySensitiveDetector::analysisFile = nullptr;
bool MySensitiveDetector::analysisFileCreated = false;
std::mutex MySensitiveDetector::analysisFileMutex;

MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
    G4String HCname = "SensorCollection";
    collectionName.insert(HCname);
    
  
    
    SensorCollection = nullptr;  // Initialize to nullptr
}

MySensitiveDetector::~MySensitiveDetector()
{
   
}



void MySensitiveDetector::Initialize(G4HCofThisEvent* HCE) {
    if (!HCE) {
        G4cerr << "Error: HCE is null in Initialize!" << G4endl;
        return;
    }
    
    SensorCollection = new SensorHitsCollection(SensitiveDetectorName, collectionName[0]);
    
    static G4int HCID = -1;
    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    
    HCE->AddHitsCollection(HCID, SensorCollection);
    
    // Create analysis file ONCE with thread safety
    std::lock_guard<std::mutex> lock(analysisFileMutex);
    if (!analysisFileCreated && !analysisFile) {
        analysisFile = new TFile("sensor_response_fits.root", "RECREATE");
        if (analysisFile && !analysisFile->IsZombie()) {
            analysisFileCreated = true;
            G4cout << "Created analysis file: sensor_response_fits.root" << G4endl;
        } else {
            G4cerr << "Failed to create analysis file!" << G4endl;
            if (analysisFile) {
                delete analysisFile;
                analysisFile = nullptr;
            }
        }
    }
}

void MySensitiveDetector::SaveToCSV(G4ThreeVector posPhotons, G4int copyNo) {
    static std::ofstream outputFile("photon_hits.csv");
    if (outputFile.is_open()) {
        outputFile << posPhotons.x()/mm << "," 
                   << posPhotons.y()/mm << "," 
                   << posPhotons.z()/mm << "," 
                   << copyNo << "\n";
    } else {
        static bool warned = false;
        if (!warned) {
            G4cerr << "Failed to open CSV file!" << G4endl;
            warned = true;
        }
    }
}


void MySensitiveDetector::RecordSensorData(G4String volumeName, int ntupleIndex, 
                                          double posX, double posY, int event, int copyNo) {
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    if (!man) return;
    
    if (ntupleIndex >= 0 && ntupleIndex < 4) {
        man->FillNtupleDColumn(ntupleIndex, 0, posX);
        man->FillNtupleDColumn(ntupleIndex, 1, posY);
        man->FillNtupleIColumn(ntupleIndex, 2, event);
        man->FillNtupleIColumn(ntupleIndex, 3, copyNo);
        man->AddNtupleRow(ntupleIndex);
        
        // Fill histogram with sensor index
        if (volumeName == "sensor_Vol1") {
            man->FillH1(0, copyNo);
        }
        else if (volumeName == "sensor_Vol3") {
            man->FillH1(1, copyNo - 25);
        }
        else if (volumeName == "sensor_Vol2") {
            man->FillH1(2, copyNo - 50);
        }
        else if (volumeName == "sensor_Vol4") {
            man->FillH1(3, copyNo - 75);
        }
    }
}

std::pair<double, double> MySensitiveDetector::CalculateCenterOfGravity(
    const std::map<G4int, std::vector<G4double>>& rightXPos,
    const std::map<G4int, std::vector<G4double>>& rightYPos,
    const std::map<G4int, std::vector<G4double>>& leftXPos,
    const std::map<G4int, std::vector<G4double>>& leftYPos,
    const std::map<G4int, std::vector<G4double>>& topXPos,
    const std::map<G4int, std::vector<G4double>>& topYPos,
    const std::map<G4int, std::vector<G4double>>& bottomXPos,
    const std::map<G4int, std::vector<G4double>>& bottomYPos) {
    
    double x_sum = 0.0, y_sum = 0.0;
    double total_weight = 0.0;
    
    // Process right sensors - using ACTUAL hit positions
    for (const auto& [copyNo, xPosVec] : rightXPos) {
        for (const auto& x_pos : xPosVec) {
            x_sum += x_pos;
            total_weight += 1.0;
        }
    }
    for (const auto& [copyNo, yPosVec] : rightYPos) {
        for (const auto& y_pos : yPosVec) {
            y_sum += y_pos;
        }
    }
    
    // Process left sensors - using ACTUAL hit positions
    for (const auto& [copyNo, xPosVec] : leftXPos) {
        for (const auto& x_pos : xPosVec) {
            x_sum += x_pos;
            total_weight += 1.0;
        }
    }
    for (const auto& [copyNo, yPosVec] : leftYPos) {
        for (const auto& y_pos : yPosVec) {
            y_sum += y_pos;
        }
    }
    
    // Process top sensors - using ACTUAL hit positions
    for (const auto& [copyNo, xPosVec] : topXPos) {
        for (const auto& x_pos : xPosVec) {
            x_sum += x_pos;
            total_weight += 1.0;
        }
    }
    for (const auto& [copyNo, yPosVec] : topYPos) {
        for (const auto& y_pos : yPosVec) {
            y_sum += y_pos;
        }
    }
    
    // Process bottom sensors - using ACTUAL hit positions
    for (const auto& [copyNo, xPosVec] : bottomXPos) {
        for (const auto& x_pos : xPosVec) {
            x_sum += x_pos;
            total_weight += 1.0;
        }
    }
    for (const auto& [copyNo, yPosVec] : bottomYPos) {
        for (const auto& y_pos : yPosVec) {
            y_sum += y_pos;
        }
    }
    
    if (total_weight > 0) {
        return std::make_pair(x_sum / total_weight, y_sum / total_weight);
    }
    return std::make_pair(0.0, 0.0);
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

std::pair<double, double> MySensitiveDetector::PerformSensorGaussianFit(
    const std::map<G4int, std::vector<G4double>>& rightXPos,
    const std::map<G4int, std::vector<G4double>>& rightYPos,
    const std::map<G4int, std::vector<G4double>>& leftXPos,
    const std::map<G4int, std::vector<G4double>>& leftYPos,
    const std::map<G4int, std::vector<G4double>>& topXPos,
    const std::map<G4int, std::vector<G4double>>& topYPos,
    const std::map<G4int, std::vector<G4double>>& bottomXPos,
    const std::map<G4int, std::vector<G4double>>& bottomYPos,
    G4int eventID) {
    
    double reconstructedX = 0.0;
    double reconstructedY = 0.0;
    
    // Variables for equation 3.1 parameters
    double Px1 = 0.0, Nx1 = 0.0, sigmaX1 = 0.0;  // Top array
    double Px2 = 0.0, Nx2 = 0.0, sigmaX2 = 0.0;  // Bottom array
    double Py1 = 0.0, Ny1 = 0.0, sigmaY1 = 0.0;  // Right array  
    double Py2 = 0.0, Ny2 = 0.0, sigmaY2 = 0.0;  // Left array
    
    bool canFitX = false;
    bool canFitY = false;
    
    // === Fit TOP sensors (for X-coordinate reconstruction) ===
    // Fit ACTUAL X positions of photons hitting top sensors
    if (!topXPos.empty()) {
        // Create histogram of ACTUAL photon X positions
        TH1D histTop(Form("histTopX_evt%d", eventID), "Top Photon X Positions", 
                     50, -50, 50);
        
        // Fill with ACTUAL photon X positions from debug output
        for (const auto& [copyNo, xPosVec] : topXPos) {
            for (const auto& x_pos : xPosVec) {
                histTop.Fill(x_pos);
                Nx1++;
            }
        }
        
        // Fit Gaussian if we have enough hits
        if (histTop.GetEntries() >= 5) {
            TF1 gaussFitTop(Form("gaussFitTopX_%d", eventID), "gaus", -50.0, 50.0);
            
            // Initial parameters from histogram
            double meanInit = histTop.GetMean();
            double rmsInit = histTop.GetRMS();
            double maxVal = histTop.GetMaximum();
            
            // Use reasonable constraints
            double initialSigma = std::max(1.0, std::min(rmsInit, 10.0));
            
            gaussFitTop.SetParameters(maxVal, meanInit, initialSigma);
            gaussFitTop.SetParLimits(0, 0.1, maxVal * 2.0);
            gaussFitTop.SetParLimits(1, -50.0, 50.0);
            gaussFitTop.SetParLimits(2, 1.0, 15.0);  // Sigma 1-15 mm
            
            if (histTop.Fit(&gaussFitTop, "QN0") == 0) {
                double fitMean = gaussFitTop.GetParameter(1);
                double fitSigma = gaussFitTop.GetParameter(2);
                
                // Accept only realistic fits
                if (std::abs(fitMean) < 55.0 && fitSigma >= 0.5 && fitSigma <= 30.0) {
                    Px1 = fitMean;
                    sigmaX1 = fitSigma;
                    canFitX = true;
                } else {
                    // Use constrained histogram values
                    Px1 = std::min(std::max(histTop.GetMean(), -50.0), 50.0);
                    sigmaX1 = std::min(std::max(histTop.GetRMS(), 2.0), 10.0);
                    canFitX = true;
                }
            } else {
                // Fit failed, use constrained histogram values
                Px1 = std::min(std::max(histTop.GetMean(), -50.0), 50.0);
                sigmaX1 = std::min(std::max(histTop.GetRMS(), 2.0), 10.0);
                canFitX = true;
            }
            
            // Final sigma constraint
            sigmaX1 = std::min(std::max(sigmaX1, 1.0), 25.0);
        }
    }
    
    // === Fit BOTTOM sensors (for X-coordinate reconstruction) ===
    // Fit ACTUAL X positions of photons hitting bottom sensors
    if (!bottomXPos.empty()) {
        TH1D histBottom(Form("histBottomX_evt%d", eventID), "Bottom Photon X Positions", 
                        50, -50, 50);
        
        for (const auto& [copyNo, xPosVec] : bottomXPos) {
            for (const auto& x_pos : xPosVec) {
                histBottom.Fill(x_pos);
                Nx2++;
            }
        }
        
        if (histBottom.GetEntries() >= 5) {
            TF1 gaussFitBottom(Form("gaussFitBottomX_%d", eventID), "gaus", -50.0, 50.0);
            
            double meanInit = histBottom.GetMean();
            double rmsInit = histBottom.GetRMS();
            double maxVal = histBottom.GetMaximum();
            
            double initialSigma = std::max(1.0, std::min(rmsInit, 10.0));
            
            gaussFitBottom.SetParameters(maxVal, meanInit, initialSigma);
            gaussFitBottom.SetParLimits(0, 0.1, maxVal * 2.0);
            gaussFitBottom.SetParLimits(1, -50.0, 50.0);
            gaussFitBottom.SetParLimits(2, 1.0, 15.0);
            
            if (histBottom.Fit(&gaussFitBottom, "QN0") == 0) {
                double fitMean = gaussFitBottom.GetParameter(1);
                double fitSigma = gaussFitBottom.GetParameter(2);
                
                if (std::abs(fitMean) < 55.0 && fitSigma >= 0.5 && fitSigma <= 30.0) {
                    Px2 = fitMean;
                    sigmaX2 = fitSigma;
                    canFitX = true;
                } else {
                    Px2 = std::min(std::max(histBottom.GetMean(), -50.0), 50.0);
                    sigmaX2 = std::min(std::max(histBottom.GetRMS(), 2.0), 10.0);
                    canFitX = true;
                }
            } else {
                Px2 = std::min(std::max(histBottom.GetMean(), -50.0), 50.0);
                sigmaX2 = std::min(std::max(histBottom.GetRMS(), 2.0), 10.0);
                canFitX = true;
            }
            
            sigmaX2 = std::min(std::max(sigmaX2, 1.0), 25.0);
        }
    }
    
    // === Fit RIGHT sensors (for Y-coordinate reconstruction) ===
    // Fit ACTUAL Y positions of photons hitting right sensors
    if (!rightYPos.empty()) {
        TH1D histRight(Form("histRightY_evt%d", eventID), "Right Photon Y Positions", 
                       50, -50, 50);
        
        for (const auto& [copyNo, yPosVec] : rightYPos) {
            for (const auto& y_pos : yPosVec) {
                histRight.Fill(y_pos);
                Ny1++;
            }
        }
        
        if (histRight.GetEntries() >= 5) {
            TF1 gaussFitRight(Form("gaussFitRightY_%d", eventID), "gaus", -50.0, 50.0);
            
            double meanInit = histRight.GetMean();
            double rmsInit = histRight.GetRMS();
            double maxVal = histRight.GetMaximum();
            
            double initialSigma = std::max(1.0, std::min(rmsInit, 10.0));
            
            gaussFitRight.SetParameters(maxVal, meanInit, initialSigma);
            gaussFitRight.SetParLimits(0, 0.1, maxVal * 2.0);
            gaussFitRight.SetParLimits(1, -50.0, 50.0);
            gaussFitRight.SetParLimits(2, 1.0, 15.0);
            
            if (histRight.Fit(&gaussFitRight, "QN0") == 0) {
                double fitMean = gaussFitRight.GetParameter(1);
                double fitSigma = gaussFitRight.GetParameter(2);
                
                if (std::abs(fitMean) < 55.0 && fitSigma >= 0.5 && fitSigma <= 30.0) {
                    Py1 = fitMean;
                    sigmaY1 = fitSigma;
                    canFitY = true;
                } else {
                    Py1 = std::min(std::max(histRight.GetMean(), -50.0), 50.0);
                    sigmaY1 = std::min(std::max(histRight.GetRMS(), 2.0), 10.0);
                    canFitY = true;
                }
            } else {
                Py1 = std::min(std::max(histRight.GetMean(), -50.0), 50.0);
                sigmaY1 = std::min(std::max(histRight.GetRMS(), 2.0), 10.0);
                canFitY = true;
            }
            
            sigmaY1 = std::min(std::max(sigmaY1, 1.0), 25.0);
        }
    }
    
    // === Fit LEFT sensors (for Y-coordinate reconstruction) ===
    // Fit ACTUAL Y positions of photons hitting left sensors
    if (!leftYPos.empty()) {
        TH1D histLeft(Form("histLeftY_evt%d", eventID), "Left Photon Y Positions", 
                      50, -50, 50);
        
        for (const auto& [copyNo, yPosVec] : leftYPos) {
            for (const auto& y_pos : yPosVec) {
                histLeft.Fill(y_pos);
                Ny2++;
            }
        }
        
        if (histLeft.GetEntries() >= 5) {
            TF1 gaussFitLeft(Form("gaussFitLeftY_%d", eventID), "gaus", -50.0, 50.0);
            
            double meanInit = histLeft.GetMean();
            double rmsInit = histLeft.GetRMS();
            double maxVal = histLeft.GetMaximum();
            
            double initialSigma = std::max(1.0, std::min(rmsInit, 10.0));
            
            gaussFitLeft.SetParameters(maxVal, meanInit, initialSigma);
            gaussFitLeft.SetParLimits(0, 0.1, maxVal * 2.0);
            gaussFitLeft.SetParLimits(1, -50.0, 50.0);
            gaussFitLeft.SetParLimits(2, 1.0, 15.0);
            
            if (histLeft.Fit(&gaussFitLeft, "QN0") == 0) {
                double fitMean = gaussFitLeft.GetParameter(1);
                double fitSigma = gaussFitLeft.GetParameter(2);
                
                if (std::abs(fitMean) < 55.0 && fitSigma >= 0.5 && fitSigma <= 30.0) {
                    Py2 = fitMean;
                    sigmaY2 = fitSigma;
                    canFitY = true;
                } else {
                    Py2 = std::min(std::max(histLeft.GetMean(), -50.0), 50.0);
                    sigmaY2 = std::min(std::max(histLeft.GetRMS(), 2.0), 10.0);
                    canFitY = true;
                }
            } else {
                Py2 = std::min(std::max(histLeft.GetMean(), -50.0), 50.0);
                sigmaY2 = std::min(std::max(histLeft.GetRMS(), 2.0), 10.0);
                canFitY = true;
            }
            
            sigmaY2 = std::min(std::max(sigmaY2, 1.0), 25.0);
        }
    }
    
    // === APPLY THE LITERATURE FORMULA (Equation 3.1) ===
    
    // For X-coordinate: Weighted mean of Top and Bottom arrays
    if (canFitX && Nx1 > 0 && Nx2 > 0 && sigmaX1 > 0 && sigmaX2 > 0) {
        // Equation 3.1 for X:
        double numeratorX = (Px1 * Nx1 / sigmaX1) + (Px2 * Nx2 / sigmaX2);
        double denominatorX = (Nx1 / sigmaX1) + (Nx2 / sigmaX2);
        
        if (denominatorX != 0) {
            reconstructedX = numeratorX / denominatorX;
        }
    }
    
    // For Y-coordinate: Weighted mean of Right and Left arrays
    if (canFitY && Ny1 > 0 && Ny2 > 0 && sigmaY1 > 0 && sigmaY2 > 0) {
        // Equation 3.1 for Y:
        double numeratorY = (Py1 * Ny1 / sigmaY1) + (Py2 * Ny2 / sigmaY2);
        double denominatorY = (Ny1 / sigmaY1) + (Ny2 / sigmaY2);
        
        if (denominatorY != 0) {
            reconstructedY = numeratorY / denominatorY;
        }
    }
    
    // Check if reconstructed positions are physically reasonable
    if (std::abs(reconstructedX) > 50.0 || std::abs(reconstructedY) > 50.0) {
        reconstructedX = 0.0;
        reconstructedY = 0.0;
    }
    
    // Add detailed debug output for successful reconstructions
    if (canFitX && canFitY && reconstructedX != 0.0 && reconstructedY != 0.0) {
        if (eventID % 10 == 0) {  // Show every 10th event with successful reconstruction
            G4cout << "Event " << eventID << " Gaussian reconstruction:" << G4endl;
            G4cout << "  Top X: P=" << Px1 << " mm, N=" << Nx1 << ", σ=" << sigmaX1 << " mm" << G4endl;
            G4cout << "  Bottom X: P=" << Px2 << " mm, N=" << Nx2 << ", σ=" << sigmaX2 << " mm" << G4endl;
            G4cout << "  Right Y: P=" << Py1 << " mm, N=" << Ny1 << ", σ=" << sigmaY1 << " mm" << G4endl;
            G4cout << "  Left Y: P=" << Py2 << " mm, N=" << Ny2 << ", σ=" << sigmaY2 << " mm" << G4endl;
            G4cout << "  Reconstructed: (" << reconstructedX << ", " << reconstructedY << ") mm" << G4endl;
        }
    }
    
    
    return std::make_pair(reconstructedX, reconstructedY);
}
G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory*) {
    G4Track *track = aStep->GetTrack();
    if (!track) return false;
    
    G4ParticleDefinition* particle = track->GetDefinition();
    if (!particle) return false;
    
    if (particle->GetParticleName() != "opticalphoton") {
        return false;
    }
    
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    if (!postStepPoint) return false;
    
    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    if (!touchable) return false;
    
    G4ThreeVector posPhotons = postStepPoint->GetPosition();
    G4int copyNo = touchable->GetCopyNumber();
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4String volumeName = touchable->GetVolume()->GetName();
    
    // ADD DEBUG: Always print first few hits
    static int hitCount = 0;
    if (hitCount < 10) {
        hitCount++;
        G4cout << "DEBUG Hit " << hitCount << ": Event " << evt 
               << " Photon in " << volumeName 
               << ", copyNo=" << copyNo
               << ", ACTUAL pos=(" << posPhotons.x()/mm << ", " 
               << posPhotons.y()/mm << ", " << posPhotons.z()/mm << ") mm"
               << G4endl;
    }
    
    // Create hit
    if (!SensorCollection) {
        G4cerr << "ERROR: SensorCollection not initialized!" << G4endl;
        return false;
    }
    
    SensorHit* aSensorHit = new SensorHit();
    aSensorHit->SetSensorPosition(posPhotons);
    aSensorHit->SetSensorEnergy(track->GetKineticEnergy());
    aSensorHit->SetVolumeName(volumeName);
    aSensorHit->SetSensorNumber(copyNo);  
    
    SensorCollection->insert(aSensorHit);
    
    // Record data - using ACTUAL photon hit positions 
    double actualX = posPhotons.x()/mm;
    double actualY = posPhotons.y()/mm;
    
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    if (!man) {
        G4cerr << "ERROR: AnalysisManager not available!" << G4endl;
        return true;
    }
    
    // Determine correct ntuple index and record ACTUAL positions
    if (volumeName == "sensor_Vol1") {
        
        RecordSensorData(volumeName, 0, actualX, actualY, evt, copyNo);
    } 
    else if (volumeName == "sensor_Vol2") {
        
        RecordSensorData(volumeName, 1, actualX, actualY, evt, copyNo);
    } 
    else if (volumeName == "sensor_Vol3") {
       
        RecordSensorData(volumeName, 2, actualX, actualY, evt, copyNo);
    } 
    else if (volumeName == "sensor_Vol4") {
    
        RecordSensorData(volumeName, 3, actualX, actualY, evt, copyNo);
    }
    else {
        G4cout << "WARNING: Unknown volume: " << volumeName << G4endl;
    }
    
    // Optional: Save to CSV
    SaveToCSV(posPhotons, copyNo);
    
    return true;
}

void MySensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {
    if (!SensorCollection) {
        G4cerr << "WARNING: SensorCollection is null in EndOfEvent!" << G4endl;
        return;
    }
    
    G4int nHits = SensorCollection->entries();
    G4int evtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    // DEBUG: Print event info
    if (evtID < 10 || nHits > 0) {
        G4cout << "Event " << evtID << ": " << nHits << " photon hits" << G4endl;
    }
    
    if (nHits > 0) {
        // Collect ACTUAL hit positions
        std::map<G4int, std::vector<G4double>> rightX, rightY;
        std::map<G4int, std::vector<G4double>> leftX, leftY;
        std::map<G4int, std::vector<G4double>> topX, topY;
        std::map<G4int, std::vector<G4double>> bottomX, bottomY;
        
        // Process all hits to collect actual positions
        for (G4int i = 0; i < nHits; i++) {
            SensorHit* hit = (*SensorCollection)[i];
            if (!hit) continue;
            
            G4int copyNo = hit->GetSensorNumber();
            G4String volumeName = hit->GetVolumeName();
            G4ThreeVector pos = hit->GetSensorPosition();
            G4double x = pos.x()/mm;
            G4double y = pos.y()/mm;
            
            if (volumeName == "sensor_Vol1") { // Right
                rightX[copyNo].push_back(x);
                rightY[copyNo].push_back(y);
            }
            else if (volumeName == "sensor_Vol2") { // Left
                leftX[copyNo].push_back(x);
                leftY[copyNo].push_back(y);
            }
            else if (volumeName == "sensor_Vol3") { // Top
                topX[copyNo].push_back(x);
                topY[copyNo].push_back(y);
            }
            else if (volumeName == "sensor_Vol4") { // Bottom
                bottomX[copyNo].push_back(x);
                bottomY[copyNo].push_back(y);
            }
        }
        
        // Debug output
        if (evtID % 100 == 0) {
            G4cout << "Event " << evtID << " hits per sensor type: "
                   << "Right=" << rightX.size() << " sensors, "
                   << "Left=" << leftX.size() << " sensors, "
                   << "Top=" << topX.size() << " sensors, "
                   << "Bottom=" << bottomX.size() << " sensors"
                   << G4endl;
        }
        
        G4AnalysisManager* man = G4AnalysisManager::Instance();
        if (!man) return;
        
        // Method 1: Gaussian fit on ACTUAL hit positions
        bool canFitX = (!topX.empty() && !bottomX.empty());
        bool canFitY = (!rightY.empty() && !leftY.empty());
        
        if (canFitX && canFitY) {
            auto [xGauss, yGauss] = PerformSensorGaussianFit(
                rightX, rightY, leftX, leftY, 
                topX, topY, bottomX, bottomY, evtID);
            
            // Only store if fit produced reasonable values
            if (std::abs(xGauss) < 50.0 && std::abs(yGauss) < 50.0) {
                // Store Gaussian fit results
                man->FillH2(0, xGauss, yGauss);  // Gaussian fit positions
                man->FillNtupleDColumn(4, 0, xGauss);
                man->FillNtupleDColumn(4, 1, yGauss);
                man->AddNtupleRow(4);
                
                if (evtID % 100 == 0) {
                    G4cout << "Event " << evtID << " Gaussian fit: ("
                           << xGauss << ", " << yGauss << ") mm" << G4endl;
                }
            } else if (evtID % 100 == 0) {
                G4cout << "Event " << evtID << " Gaussian fit rejected: ("
                       << xGauss << ", " << yGauss << ") mm (outside detector)" << G4endl;
            }
        }
        
        // Method 2: Center of gravity using ACTUAL hit positions
        auto [xCOG, yCOG] = CalculateCenterOfGravity(
            rightX, rightY, leftX, leftY,
            topX, topY, bottomX, bottomY);
        
        // Store COG results
        man->FillH2(1, xCOG, yCOG);  // Center of gravity positions
        man->FillNtupleDColumn(5, 0, xCOG);
        man->FillNtupleDColumn(5, 1, yCOG);
        man->AddNtupleRow(5);
        
        if (evtID % 100 == 0) {
            G4cout << "Event " << evtID << " Center of gravity: ("
                   << xCOG << ", " << yCOG << ") mm" << G4endl;
        }
    }
    
    
}
