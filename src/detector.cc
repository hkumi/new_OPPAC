#include "detector.hh"
#include <numeric>
#include <map>
#include <vector>
#include <fstream>  // For writing to a file
#include <sstream>  // For string stream if you need formatting
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include "TF1.h"
#include "G4SDManager.hh"
#include "G4String.hh"
#include "G4HCofThisEvent.hh"
#include "SensorHit.hh"

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
 //   G4cout << "SensitiveDetectorName: " << SensitiveDetectorName << G4endl;
//G4cout << "Collection Name: " << collectionName[0] << G4endl;

    // Register the collection with the event
    static G4int HCID = -1; // Static to ensure it's set only once
    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    if (HCE) {
        HCE->AddHitsCollection(HCID, SensorCollection);
    } else {
        G4cerr << "Error: HCE is null in Initialize!" << G4endl;
    }
}



// Define the function to save positions and copy numbers to a CSV file
void MySensitiveDetector::SaveToCSV(G4ThreeVector posPhotons, G4int copyNo) {
    // Open the CSV file in append mode (to keep adding rows)
    std::ofstream outputFile("photon_positions_4.csv", std::ios::app);  

    // Check if the file is open
    if (outputFile.is_open()) {
        // Write the data to the CSV file
        // Each entry is written in a new line with positions (x, y, z) and copy number
        outputFile << posPhotons.x()/mm << "," 
                   << posPhotons.y()/mm << "," 
                   << posPhotons.z()/mm << "," 
                   << copyNo << "\n";
    } else {
        G4cerr << "Failed to open file!" << G4endl;
    }

    // Close the file
    outputFile.close();
}

void MySensitiveDetector::RecordSensorData(const std::string& volumeName, int ntupleIndex, double posX, double posY, int event, int copyNo, G4AnalysisManager* man) {
    // Fill ntuple with data
    man->FillNtupleDColumn(ntupleIndex, 0, posX );   // X position
    man->FillNtupleDColumn(ntupleIndex, 1, posY );   // Y position
    man->FillNtupleDColumn(ntupleIndex, 2, event);       // Event number
    man->FillNtupleIColumn(ntupleIndex, 3, copyNo);      // Sensor copy number
    man->AddNtupleRow(ntupleIndex);
}

// Calculate the weighted mean for the x-coordinate
double MySensitiveDetector::calculateWeightedMeanX(double P_x1, double N_x1, double sigma_x1, 
                                                   double P_x2, double N_x2, double sigma_x2) {
    // Check for invalid sigma values
    if (sigma_x1 == 0 || sigma_x2 == 0) {
        std::cout << "Error: Sigma values cannot be zero!" << std::endl;
        return -1; // Return a sentinel value to indicate error
    }

    // Calculate the weighted mean
    double numerator = (P_x1 * N_x1 / sigma_x1) + (P_x2 * N_x2 / sigma_x2);
    double denominator = (N_x1 / sigma_x1) + (N_x2 / sigma_x2);

    return numerator / denominator;
}

// Calculate the weighted mean for the y-coordinate
double MySensitiveDetector::calculateWeightedMeanY(double P_y1, double N_y1, double sigma_y1, 
                                                   double P_y2, double N_y2, double sigma_y2) {
    // Check for invalid sigma values
    if (sigma_y1 == 0 || sigma_y2 == 0) {
        std::cout << "Error: Sigma values cannot be zero!" << std::endl;
        return -1; // Return a sentinel value to indicate error
    }

    // Calculate the weighted mean
    double numerator = (P_y1 * N_y1 / sigma_y1) + (P_y2 * N_y2 / sigma_y2);
    double denominator = (N_y1 / sigma_y1) + (N_y2 / sigma_y2);

    return numerator / denominator;
}




void MySensitiveDetector::FitHistogram(const std::vector<double>& xPosition1, const std::vector<double>& xPosition2, 
                                       const std::vector<double>& yPosition1, const std::vector<double>& yPosition2) {
 
    if (xPosition1.empty() || xPosition2.empty() || yPosition1.empty() || yPosition2.empty()) return;

    int bins = 50;

    // Process X positions
    double minX1 = *std::min_element(xPosition1.begin(), xPosition1.end());
    double maxX1 = *std::max_element(xPosition1.begin(), xPosition1.end());
    double minX2 = *std::min_element(xPosition2.begin(), xPosition2.end());
    double maxX2 = *std::max_element(xPosition2.begin(), xPosition2.end());

    static int histCounterX1 = 0, histCounterX2 = 0;
    TH1D* histX1 = new TH1D(Form("histX1_%d", histCounterX1++), "X Sensor1 Distribution", bins, minX1, maxX1);
    TH1D* histX2 = new TH1D(Form("histX2_%d", histCounterX2++), "X Sensor2 Distribution", bins, minX2, maxX2);

    for (double x : xPosition1) histX1->Fill(x);
    for (double x : xPosition2) histX2->Fill(x);

    TF1* gaussFitX1 = new TF1("gaussFitX1", "gaus", minX1, maxX1);
    gaussFitX1->SetLineColor(kRed);
    histX1->Fit(gaussFitX1, "R");

    TF1* gaussFitX2 = new TF1("gaussFitX2", "gaus", minX2, maxX2);
    gaussFitX2->SetLineColor(kRed);
    histX2->Fit(gaussFitX2, "R");

    double meanX1 = gaussFitX1->GetParameter(1);
    double sigmaX1 = gaussFitX1->GetParameter(2);
    double amplitudeX1 = gaussFitX1->GetParameter(0);

    double meanX2 = gaussFitX2->GetParameter(1);
    double sigmaX2 = gaussFitX2->GetParameter(2);
    double amplitudeX2 = gaussFitX2->GetParameter(0);

    // Process Y positions
    double minY1 = *std::min_element(yPosition1.begin(), yPosition1.end());
    double maxY1 = *std::max_element(yPosition1.begin(), yPosition1.end());
    double minY2 = *std::min_element(yPosition2.begin(), yPosition2.end());
    double maxY2 = *std::max_element(yPosition2.begin(), yPosition2.end());

    static int histCounterY1 = 0, histCounterY2 = 0;
    TH1D* histY1 = new TH1D(Form("histY1_%d", histCounterY1++), "Y Sensor1 Distribution", bins, minY1, maxY1);
    TH1D* histY2 = new TH1D(Form("histY2_%d", histCounterY2++), "Y Sensor2 Distribution", bins, minY2, maxY2);

    for (double y : yPosition1) histY1->Fill(y);
    for (double y : yPosition2) histY2->Fill(y);

    TF1* gaussFitY1 = new TF1("gaussFitY1", "gaus", minY1, maxY1);
    gaussFitY1->SetLineColor(kBlue);
    histY1->Fit(gaussFitY1, "R");

    TF1* gaussFitY2 = new TF1("gaussFitY2", "gaus", minY2, maxY2);
    gaussFitY2->SetLineColor(kBlue);
    histY2->Fit(gaussFitY2, "R");

    double meanY1 = gaussFitY1->GetParameter(1);
    double sigmaY1 = gaussFitY1->GetParameter(2);
    double amplitudeY1 = gaussFitY1->GetParameter(0);

    double meanY2 = gaussFitY2->GetParameter(1);
    double sigmaY2 = gaussFitY2->GetParameter(2);
    double amplitudeY2 = gaussFitY2->GetParameter(0);

    // Compute weighted means
    double x_pos = calculateWeightedMeanX(meanX1, amplitudeX1, sigmaX1, meanX2, amplitudeX2, sigmaX2);
    double y_pos = calculateWeightedMeanY(meanY1, amplitudeY1, sigmaY1, meanY2, amplitudeY2, sigmaY2);

    // Debugging outputs
    //G4cout << "Weighted Mean X: " << x_pos << G4endl;
    //G4cout << "Weighted Mean Y: " << y_pos << G4endl;

    // Store the results in ntuple
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->FillH2(0, x_pos, y_pos);
    man->FillNtupleDColumn(4,0, x_pos);
    man->AddNtupleRow(4);

    man->FillNtupleDColumn(5,0, y_pos);
    man->AddNtupleRow(5);


    

    

    

    // Cleanup
    delete histX1; delete histX2;
    delete gaussFitX1; delete gaussFitX2;
    delete histY1; delete histY2;
    delete gaussFitY1; delete gaussFitY2;
}





G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory*)        
{

    SensorHit* aSensorHit = new SensorHit();

    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();    
    G4Track *track = aStep->GetTrack();

    //track->SetTrackStatus(fStopAndKill);
    G4ParticleDefinition*  particle =aStep->GetTrack()->GetDefinition();
    if (particle->GetParticleName() == "opticalphoton"){   
       G4StepPoint *preStepPoint = aStep->GetPreStepPoint();//used when the neutron enters the detector
       G4StepPoint *postStepPoint = aStep->GetPostStepPoint();//used when the neutron leaves the detector
       const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

       G4ThreeVector posPhotons = postStepPoint->GetPosition();//accessing the position       
       G4String particleName = particle->GetParticleName();
       //G4cout << "This is the particle in the sensor: " << particleName << G4endl;       

  
       G4int copyNo = touchable->GetCopyNumber();
       G4VPhysicalVolume *physVol = touchable->GetVolume();
       G4ThreeVector posDetector = physVol->GetTranslation();

       G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
       // Declare a vector to store the copy numbers
       std::vector<G4int> copyNumbers;

       G4String volumeName = touchable->GetVolume()->GetName(); // name of the volume
      // G4cout << "Photon detected in " << volumeName << "with Copy Number" << copyNo<< G4endl;
       G4AnalysisManager *man = G4AnalysisManager::Instance();

       G4double X1 = posPhotons.getX();
       G4double Y1 = posPhotons.getY();
       //std::vector<int> sensor1CopyNumbers;
      
        
       if (volumeName == "sensor_Vol1") {
          RecordSensorData(volumeName, 0, posDetector[0], posDetector[1], evt, copyNo, man);

          aSensorHit->SetSensorPosition(postStepPoint->GetPosition());
          aSensorHit->SetSensorEnergy(track->GetKineticEnergy());
          aSensorHit->SetVolumeName(volumeName);

          if (SensorCollection) {
             HitID = SensorCollection->insert(aSensorHit);
            
          } else {
            G4cerr << "SensorCollection not initialized!" << G4endl;
            delete aSensorHit; // Avoid memory leaks
            }
       } else if (volumeName == "sensor_Vol2") {
         // Create a new hit object for sensor_Vol2
         SensorHit* anotherSensorHit = new SensorHit();

         RecordSensorData(volumeName, 1, posDetector[0], posDetector[1], evt, copyNo, man);

         anotherSensorHit->SetSensorPosition(postStepPoint->GetPosition());
         anotherSensorHit->SetSensorEnergy(track->GetKineticEnergy());
         anotherSensorHit->SetVolumeName(volumeName);

         if (SensorCollection) {
            HitID = SensorCollection->insert(anotherSensorHit);

         } else {
           G4cerr << "SensorCollection not initialized!" << G4endl;
           delete anotherSensorHit; // Avoid memory leaks
           }

         } else if (volumeName == "sensor_Vol3") {
         // Create a new hit object for sensor_Vol3
         SensorHit* SensorHit3 = new SensorHit();

         RecordSensorData(volumeName, 2, posDetector[0], posDetector[1], evt, copyNo, man);

         SensorHit3->SetSensorPosition(postStepPoint->GetPosition());
         SensorHit3->SetSensorEnergy(track->GetKineticEnergy());
         SensorHit3->SetVolumeName(volumeName);

         if (SensorCollection) {
            HitID = SensorCollection->insert(SensorHit3);

         } else {
           G4cerr << "SensorCollection not initialized!" << G4endl;
           delete SensorHit3; // Avoid memory leaks
           }


           } else if (volumeName == "sensor_Vol4") {
         // Create a new hit object for sensor_Vol4
         SensorHit* SensorHit4 = new SensorHit();

         RecordSensorData(volumeName, 3, posDetector[0], posDetector[1], evt, copyNo, man);

         SensorHit4->SetSensorPosition(postStepPoint->GetPosition());
         SensorHit4->SetSensorEnergy(track->GetKineticEnergy());
         SensorHit4->SetVolumeName(volumeName);

         if (SensorCollection) {
            HitID = SensorCollection->insert(SensorHit4);

         } else {
           G4cerr << "SensorCollection not initialized!" << G4endl;
           delete SensorHit4; // Avoid memory leaks
           }


        
  
       }

         /*
            } else if (volumeName == "sensor_Vol3") {
              RecordSensorData(volumeName, 2, posDetector[0], posDetector[1], evt, copyNo, man);
              } else if (volumeName == "sensor_Vol4") {
                RecordSensorData(volumeName, 3, posDetector[0], posDetector[1], evt, copyNo, man);
                SaveToCSV(posDetector, copyNo);
                }*/

         //  FitHistogram(sensor1CopyNumbers);


       
    }     


    return true;

}


void MySensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {
    // Log the total number of hits for debugging
    G4int nHits = SensorCollection->entries();
    //G4cout << "End of Event: Number of hits in SensorCollection: " << nHits << G4endl;

    // Check if there are any hits
    if (nHits > 0) {
        std::vector<double> yPositions1, yPositions2; // For sensors 1 and 2
        std::vector<double> xPositions1, xPositions2; // For sensor 3 and 4
        std::vector<G4int> sensor1CopyNumbers, sensor2CopyNumbers;

        // Process hits
        for (G4int i = 0; i < nHits; i++) {
            SensorHit* hit = (*SensorCollection)[i];
            G4ThreeVector position = hit->GetSensorPosition(); // Get the hit position
            G4String VolumeName = hit->GetVolumeName();       // Copy number (optional)

            // Separate data based on the sensor volume name
            if (VolumeName == "sensor_Vol1") { // For sensor_Vol1
               yPositions1.push_back(position.y()/mm);
               //sensor1CopyNumbers.push_back(sensorNumber);
            } else if (VolumeName == "sensor_Vol2") { // For sensor_Vol2
              yPositions2.push_back(position.y()/mm);
              //sensor2CopyNumbers.push_back(sensorNumber);
            } else if (VolumeName == "sensor_Vol3") { // For sensor_Vol3
              xPositions1.push_back(position.x()/mm);
            } else if (VolumeName == "sensor_Vol4") { // For sensor_Vol4
              xPositions2.push_back(position.x()/mm);
            }
        }

        // Output for debugging
        G4cout << "Sensor_Vol1 Hits: " << yPositions1.size() 
               << ", Sensor_Vol2 Hits: " << yPositions2.size() << G4endl;

        // Perform Gaussian fits if we have enough data
        if (!yPositions1.empty() && !yPositions2.empty() && 
            !xPositions1.empty() && !xPositions2.empty()) {
            FitHistogram(xPositions1, xPositions2, yPositions1, yPositions2);

        }

    }
}

