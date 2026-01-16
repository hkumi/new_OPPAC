#ifndef DC_h
#define DC_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4OpticalSurface;

class DC : public G4VUserDetectorConstruction
{
  public:
    DC(G4double density, G4double collimatorLength);
    virtual ~DC();

    virtual G4VPhysicalVolume* Construct() override;
    virtual void ConstructSDandField() override;
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  private:
    void DefineMaterials();
    
    // Geometry parameters
    G4double fDensity;
    G4double fCollimatorLength;
    
    // Materials (only what we actually use)
    G4Material* fTeflon;           // Collimator material
    
    // Optical surfaces
    G4OpticalSurface* fSipmSurf;
    G4OpticalSurface* fMylarSurf;
    
    // World volume
    G4VPhysicalVolume* fWorldPhys;
    G4LogicalVolume* fWorldLog;
    
    // Scoring volume
    G4LogicalVolume* fScoringVolume;
    
    // Check overlaps
    G4bool fCheckOverlaps;
};

#endif
