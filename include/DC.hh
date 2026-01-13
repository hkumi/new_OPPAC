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
    
    // Materials
    G4Material* fGasMaterial;
    G4Material* fTeflon;
    G4Material* fHDPE;
    G4Material* fSilicon;
    G4Material* fAluminium;
    G4Material* fMylar;
    G4Material* fVacuum;
    G4Material* fLead;
    G4Material* fPolyethylene;
    G4Material* fBoratedPolyethylene;
    G4Material* CF4;    
    // Optical surfaces
    G4OpticalSurface* fSipmSurf;
    G4OpticalSurface* fMylarSurf;
    G4OpticalSurface* fTeflonSurf;
    G4OpticalSurface* fCf4SiSurface;
    
    // World volume
    G4VPhysicalVolume* fWorldPhys;
    G4LogicalVolume* fWorldLog;
    
    // Scoring volume
    G4LogicalVolume* fScoringVolume;
    
    // Check overlaps
    G4bool fCheckOverlaps;
};

#endif
