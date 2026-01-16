#include "DC.hh"
#include "detector.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialPropertiesTable.hh"

#include <cmath>

// Constructor - simplified
DC::DC(G4double density, G4double collimatorLength)
    : G4VUserDetectorConstruction(),
      fDensity(density),
      fCollimatorLength(collimatorLength),
      fTeflon(nullptr),           // Only keep Teflon
      fSipmSurf(nullptr),
      fMylarSurf(nullptr),
      fWorldPhys(nullptr),
      fWorldLog(nullptr),
      fScoringVolume(nullptr),
      fCheckOverlaps(true)
{
}

// Destructor - simplified
DC::~DC()
{
    // Clean up optical surfaces
    if (fSipmSurf) delete fSipmSurf;
    if (fMylarSurf) delete fMylarSurf;
}

void DC::DefineMaterials()
{
    G4int ncomponents, natoms;
    G4double massfraction;

    // Get elements
    auto nistManager = G4NistManager::Instance();
    auto C = nistManager->FindOrBuildElement("C");
    auto F = nistManager->FindOrBuildElement("F");
    auto Ar = nistManager->FindOrBuildElement("Ar");
    auto H = nistManager->FindOrBuildElement("H");
    auto O = nistManager->FindOrBuildElement("O");
    auto Si = nistManager->FindOrBuildElement("Si");
    auto Al = nistManager->FindOrBuildElement("Al");

   
     
    // CF4
    G4Material* CF4 = new G4Material("CF4", 3.72 * mg / cm3, 2, kStateGas, 
                                     293.15 * kelvin, 1.0 * atmosphere);
    CF4->AddElement(C, 1);
    CF4->AddElement(F, 4);

    // Ar gas
    G4Material* Ar_gas = new G4Material("Ar_gas", 1.782 * mg / cm3, 1, 
                                        kStateGas, 293.15 * kelvin, 1.0 * atmosphere);
    Ar_gas->AddElement(Ar, 1);

    // Ar:CF4 (90:10) - Scintillating gas
    G4Material* Ar_CF4 = new G4Material("Ar_CF4", 0.061 * mg / cm3, 2, 
                                        kStateGas, 293.15 * kelvin, 30.e-3 * bar);
    Ar_CF4->AddMaterial(Ar_gas, 0.9);
    Ar_CF4->AddMaterial(CF4, 0.1);

    // Mylar
    G4Material* mylar = new G4Material("mylar", 1.39 * g / cm3, 3);
    mylar->AddElement(C, 5);
    mylar->AddElement(H, 4);
    mylar->AddElement(O, 2);

    // SiPM sensor
    G4Material* sensor = new G4Material("SiPM", 2.33 * g / cm3, 2);
    sensor->AddElement(Si, 1);
    sensor->AddElement(O, 2);

    // Aluminium
    G4Material* alum = new G4Material("alum", 2.7 * g / cm3, 1);
    alum->AddElement(Al, 1);

    // HDPE
    G4Material* HDPE = new G4Material("HDPE", 0.93 * g / cm3, 2);
    HDPE->AddElement(C, 2);
    HDPE->AddElement(H, 4);

    // Optical properties --------------------------------------------------------

    // Scintillating gas properties
    const G4int nScint = 18;
    G4double gasEnergy[nScint] = {
        1.55 * eV, 1.61 * eV, 1.70 * eV, 1.78 * eV, 1.88 * eV,
        2.00 * eV, 2.11 * eV, 2.22 * eV, 2.54 * eV, 2.87 * eV,
        3.42 * eV, 3.89 * eV, 4.35 * eV, 4.72 * eV, 4.88 * eV, 
        5.40 * eV, 5.83 * eV, 6.16 * eV
    };
    
    G4double gasScintSp[nScint] = {
        0.01, 0.01, 0.04, 0.10, 0.20,
        0.26, 0.18, 0.04, 0.00, 0.02,
        0.04, 0.13, 0.30, 0.19, 0.25,
        0.03, 0.00, 0.00
    };
    
    G4double gasAbsLength[nScint];
    for (int i = 0; i < nScint; i++) {
        gasAbsLength[i] = 4 * m;
    }

    G4double gasScatLength[nScint];
    for (int i = 0; i < nScint; i++) {
        gasScatLength[i] = 5.0 * cm;
    }

    auto ArCF4_mpt = new G4MaterialPropertiesTable();
    ArCF4_mpt->AddProperty("SCINTILLATIONCOMPONENT1", gasEnergy, gasScintSp, nScint);
    ArCF4_mpt->AddProperty("RAYLEIGH", gasEnergy, gasScatLength, nScint);
    ArCF4_mpt->AddProperty("RINDEX", "Air");
    ArCF4_mpt->AddConstProperty("SCINTILLATIONYIELD", 5000 / MeV);
    ArCF4_mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    ArCF4_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 15 * ns);
    ArCF4_mpt->AddProperty("ABSLENGTH", gasEnergy, gasAbsLength, nScint);
    ArCF4_mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    Ar_CF4->SetMaterialPropertiesTable(ArCF4_mpt);

    // SiPM sensor properties
    const G4int nSiPM = 6;
    G4double SiPMEnergy[nSiPM] = {
        2.07 * eV, 2.34 * eV, 2.64 * eV, 2.75 * eV, 2.95 * eV, 3.10 * eV
    };
    G4double SiPMQE[nSiPM] = {
        0.25, 0.35, 0.45, 0.48, 0.50, 0.45
    };
    G4double SiPMrIndex[nSiPM] = {
        3.4, 3.4, 3.4, 3.4, 3.4, 3.4
    };
    
    auto SiPM_mpt = new G4MaterialPropertiesTable();
    SiPM_mpt->AddProperty("RINDEX", SiPMEnergy, SiPMrIndex, nSiPM);
    SiPM_mpt->AddProperty("EFFICIENCY", SiPMEnergy, SiPMQE, nSiPM);
    sensor->SetMaterialPropertiesTable(SiPM_mpt);

    // Mylar & Al properties
    const G4int nAl = 6;
    G4double AlEnergy[nAl] = {
        2.07 * eV, 2.34 * eV, 2.64 * eV, 2.75 * eV, 2.95 * eV, 3.10 * eV
    };
    G4double AlrIndex[nAl] = {
        1.373, 1.373, 1.373, 1.373, 1.373, 1.373
    };
    G4double AlAbsLength[nAl] = {
        4 * cm, 4 * cm, 4 * cm, 4 * cm, 4 * cm, 4 * cm
    };
    
    auto Al_mtp = new G4MaterialPropertiesTable();
    Al_mtp->AddProperty("RINDEX", AlEnergy, AlrIndex, nAl);
    Al_mtp->AddProperty("ABSLENGTH", AlEnergy, AlAbsLength, nAl);
    alum->SetMaterialPropertiesTable(Al_mtp);
    mylar->SetMaterialPropertiesTable(Al_mtp);

    // Teflon - collimator material
    const G4int TNbEntries = 4;
    G4Material* Teflon = new G4Material("Teflon", 2.2 * g / cm3, 2, kStateSolid);
    Teflon->AddElement(C, 0.240183);
    Teflon->AddElement(F, 0.759817);

    // Teflon optical properties (highly reflective)
    G4double pdTeflonPhotonMomentum[TNbEntries] = {
        2.07 * eV, 2.34 * eV, 2.64 * eV, 3.10 * eV
    };
    G4double pdTeflonRefractiveIndex[TNbEntries] = {1.34, 1.34, 1.34, 1.34};
    G4double pdTeflonReflectivity[TNbEntries] = {0.95, 0.95, 0.95, 0.95};
    G4double pdTeflonAbsLength[TNbEntries] = {
        1.0 * m, 1.0 * m, 1.0 * m, 1.0 * m
    };

    G4MaterialPropertiesTable* pTeflonPropertiesTable = new G4MaterialPropertiesTable();
    pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, 
                                        pdTeflonRefractiveIndex, TNbEntries);
    pTeflonPropertiesTable->AddProperty("ABSLENGTH", pdTeflonPhotonMomentum, 
                                        pdTeflonAbsLength, TNbEntries);
    pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, 
                                        pdTeflonReflectivity, TNbEntries);
    Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);

    // Store Teflon in class variable
    fTeflon = Teflon;
}

G4VPhysicalVolume* DC::Construct()
{
    DefineMaterials();
    
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    
    // Get materials
    auto atmosphere = G4Material::GetMaterial("Ar_CF4");
    auto sensor = G4Material::GetMaterial("SiPM");
    auto HDPE = G4Material::GetMaterial("HDPE");
    auto mylar = G4Material::GetMaterial("mylar");
    auto Teflon = fTeflon;  // Use stored Teflon material
    
    if (!Teflon) {
        G4cerr << "ERROR: Teflon material not found!" << G4endl;
        return nullptr;
    }
    
    // -------------------------- //
    // World volume - ArCF4 (90/10)
    // -------------------------- //
    G4double worldSize = 100 * mm;
    G4double worldDepth = 10 * mm;
    
    G4Box* worldBox = new G4Box("World", worldSize/2, worldSize/2, worldDepth/2);
    fWorldLog = new G4LogicalVolume(worldBox, atmosphere, "World");
    fWorldPhys = new G4PVPlacement(0, G4ThreeVector(), fWorldLog, "World", 
                                   0, false, 0, fCheckOverlaps);
    
    // Geometry parameters 
    G4double sipmSize = 1.0 * mm;  
    G4double cellSize = sipmSize + 0.84 * mm;     // = 1.84 mm
    G4double cellLength = fCollimatorLength;      // Use constructor parameter
    G4int nCells = 25;
    G4double collSize = nCells * cellSize;        // Total collimator size: 46 mm
    G4double sipmWidth = 1.0 * mm;
    
   
    
    // Create collimator cell - TEFLON
    G4Box* sBlock = new G4Box("cellBlock", cellLength/2, cellSize/2, cellSize/2);
    G4Box* sHole = new G4Box("cellHole", cellLength/2, sipmSize/2, sipmSize/2);
    G4SubtractionSolid* sCell = new G4SubtractionSolid("cell", sBlock, sHole, 
                                                       0, G4ThreeVector(0, 0, 0));
    G4LogicalVolume* fLCell = new G4LogicalVolume(sCell, Teflon, "Cell");
    
    // Place collimator cells (4 sides) 
    G4int copyNo = 0;
    for (G4int i = 0; i < 4; i++) {
        G4double angle = i * 90. * deg;
        auto cellRot = new G4RotationMatrix();
        cellRot->rotateZ(angle); 

        for (G4int j = 0; j < nCells; j++) {
            G4double xPos = 0, yPos = 0;
            G4double cellCenterOffset = (-nCells/2.0 + j + 0.5) * cellSize;
            
            if (i == 0) { // RIGHT side (+X)
                xPos = collSize/2 + cellLength/2;
                yPos = cellCenterOffset;
            }
            else if (i == 1) { // TOP side (+Y)
                xPos = cellCenterOffset;
                yPos = collSize/2 + cellLength/2;
            }
            else if (i == 2) { // LEFT side (-X)
                xPos = -(collSize/2 + cellLength/2);
                yPos = cellCenterOffset;
            }
            else if (i == 3) { // BOTTOM side (-Y)
                xPos = cellCenterOffset;
                yPos = -(collSize/2 + cellLength/2);
            }

            new G4PVPlacement(cellRot, G4ThreeVector(xPos, yPos, 0), fLCell, "Cell", 
                             fWorldLog, false, copyNo++, fCheckOverlaps);
        }
    }
    
    // Mylar electrodes
    G4double mylThickness = 0.012 * mm;
    G4double mylPos = cellSize;

    G4Box* sMyl = new G4Box("myl", collSize/2, collSize/2, mylThickness/2);
    G4LogicalVolume* fLMyl = new G4LogicalVolume(sMyl, mylar, "Myl");

    G4VPhysicalVolume* fPMylA = new G4PVPlacement(0, G4ThreeVector(0, 0, mylPos/2), 
                                                 fLMyl, "MylA", fWorldLog, false, 0, fCheckOverlaps);
    G4VPhysicalVolume* fPMylK = new G4PVPlacement(0, G4ThreeVector(0, 0, -mylPos/2), 
                                                 fLMyl, "MylK", fWorldLog, false, 1, fCheckOverlaps);
    
    // HDPE converter
    G4double convThickness = 0.01 * mm;
    G4double convPos = mylPos/2 + mylThickness/2 + convThickness/2;

    G4Box* sConv = new G4Box("conv", collSize/2, collSize/2, convThickness/2);
    G4LogicalVolume* fLConv = new G4LogicalVolume(sConv, HDPE, "Conv");

    new G4PVPlacement(0, G4ThreeVector(0, 0, convPos), fLConv, "Conv", 
                     fWorldLog, false, 0, fCheckOverlaps);
    
    // SiPM sensors
    G4Box* sSiPM = new G4Box("sipm", sipmWidth/2, sipmSize/2, sipmSize/2);
    
    // Create 4 logical volumes for sensors (one for each side)
    G4LogicalVolume* sipmLog1 = new G4LogicalVolume(sSiPM, sensor, "sensor_log1"); // RIGHT
    G4LogicalVolume* sipmLog2 = new G4LogicalVolume(sSiPM, sensor, "sensor_log2"); // LEFT
    G4LogicalVolume* sipmLog3 = new G4LogicalVolume(sSiPM, sensor, "sensor_log3"); // TOP
    G4LogicalVolume* sipmLog4 = new G4LogicalVolume(sSiPM, sensor, "sensor_log4"); // BOTTOM
    
    // Place sensors - align with collimator holes
    copyNo = 0;
    for (G4int i = 0; i < 4; i++) {
        G4LogicalVolume* currentSipmLog = nullptr;
        G4String volumeName;
        
        if (i == 0) { // Right side
            currentSipmLog = sipmLog1;
            volumeName = "sensor_Vol1";
        } else if (i == 1) { // Top side
            currentSipmLog = sipmLog3;
            volumeName = "sensor_Vol3";
        } else if (i == 2) { // Left side
            currentSipmLog = sipmLog2;
            volumeName = "sensor_Vol2";
        } else { // Bottom side
            currentSipmLog = sipmLog4;
            volumeName = "sensor_Vol4";
        }
        
        for (G4int j = 0; j < nCells; j++) {
            G4double xPos = 0, yPos = 0;
            G4double angle = i * 90. * deg;
            auto sipmRot = new G4RotationMatrix();
            sipmRot->rotateZ(angle); 

            G4double sensorCenterOffset = (-nCells/2.0 + j + 0.5) * cellSize;
            
            if (i == 0) { // RIGHT side
                xPos = collSize/2 + cellLength + sipmWidth/2;
                yPos = sensorCenterOffset;
            }
            else if (i == 1) { // TOP side
                xPos = sensorCenterOffset;
                yPos = collSize/2 + cellLength + sipmWidth/2;
            }
            else if (i == 2) { // LEFT side
                xPos = -(collSize/2 + cellLength + sipmWidth/2);
                yPos = sensorCenterOffset;
            }
            else if (i == 3) { // BOTTOM side
                xPos = sensorCenterOffset;
                yPos = -(collSize/2 + cellLength + sipmWidth/2);
            }
            
            new G4PVPlacement(sipmRot, G4ThreeVector(xPos, yPos, 0), 
                currentSipmLog, volumeName, fWorldLog, false, copyNo++, fCheckOverlaps);
        }
    }
    
    // Set up optical surfaces 
    auto sipmSurf = new G4OpticalSurface("sipmSurf");
    sipmSurf->SetType(dielectric_dielectric);
    sipmSurf->SetFinish(polished);
    sipmSurf->SetModel(unified);
    
    auto mylarSurf = new G4OpticalSurface("mylSurf");
    mylarSurf->SetType(dielectric_metal);
    mylarSurf->SetFinish(polished);
    mylarSurf->SetModel(unified);
    
    // Optical properties for surfaces
    const G4int nOptic = 6;
    G4double photonEnergy[nOptic] = {
        2.07 * eV, 2.34 * eV, 2.64 * eV, 2.75 * eV, 2.95 * eV, 3.10 * eV
    };
    G4double mylarRefl[nOptic] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.9};
    G4double mylarTrans[nOptic] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double SiPMeff[nOptic] = {0.25, 0.35, 0.45, 0.48, 0.50, 0.45};
    
    auto sipm_mtp = new G4MaterialPropertiesTable();
    sipm_mtp->AddProperty("EFFICIENCY", photonEnergy, SiPMeff, nOptic);
    sipmSurf->SetMaterialPropertiesTable(sipm_mtp);
    
    auto mylarS_mtp = new G4MaterialPropertiesTable();
    mylarS_mtp->AddProperty("REFLECTIVITY", photonEnergy, mylarRefl, nOptic);
    mylarS_mtp->AddProperty("TRANSMITTANCE", photonEnergy, mylarTrans, nOptic);
    mylarSurf->SetMaterialPropertiesTable(mylarS_mtp);
    
    // Apply optical surfaces
    // Teflon uses material properties, no need for separate surface
    
    // Mylar surfaces
    new G4LogicalBorderSurface("MylarA-Air_in", fWorldPhys, fPMylA, mylarSurf);
    new G4LogicalBorderSurface("MylarA-Air_out", fPMylA, fWorldPhys, mylarSurf);
    new G4LogicalBorderSurface("MylarK-Air_in", fWorldPhys, fPMylK, mylarSurf);
    new G4LogicalBorderSurface("MylarK-Air_out", fPMylK, fWorldPhys, mylarSurf);
    
    // SiPM surfaces
    new G4LogicalSkinSurface("SiPMSurface1", sipmLog1, sipmSurf);
    new G4LogicalSkinSurface("SiPMSurface2", sipmLog2, sipmSurf);
    new G4LogicalSkinSurface("SiPMSurface3", sipmLog3, sipmSurf);
    new G4LogicalSkinSurface("SiPMSurface4", sipmLog4, sipmSurf);
    
    // Visualization attributes
    fLCell->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.3))); // Yellow, transparent
    fLMyl->SetVisAttributes(G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.5))); // Grey, transparent
    fLConv->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.5))); // Green, transparent
    
    sipmLog1->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.0, 0.0))); // Red
    sipmLog2->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 0.0, 1.0))); // Blue
    sipmLog3->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.5, 0.0))); // Orange
    sipmLog4->SetVisAttributes(G4VisAttributes(G4Colour(0.5, 0.0, 0.5))); // Purple
    
    fWorldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
    
   
    
    return fWorldPhys;
}

void DC::ConstructSDandField()
{
    // Set up sensitive detectors
    MySensitiveDetector* sensDet = new MySensitiveDetector("SensitiveDetector");
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    SDman->AddNewDetector(sensDet);
    
    // Set SiPM logical volumes as sensitive
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
    
    for (auto lv : *lvStore) {
        G4String lvName = lv->GetName();
        if (lvName == "sensor_log1" || lvName == "sensor_log2" || 
            lvName == "sensor_log3" || lvName == "sensor_log4") {
            lv->SetSensitiveDetector(sensDet);
        }
    }
}
