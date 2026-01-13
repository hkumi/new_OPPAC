#ifndef SensorHit_h
#define SensorHit_h 1

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include <G4VProcess.hh>

class SensorHit : public G4VHit {
  public:
    SensorHit();
    ~SensorHit();
    SensorHit(const SensorHit&);

    const SensorHit& operator=(const SensorHit&);
    int operator==(const SensorHit&) const;

    // Declare operator new and delete (defined inline in the header file)
    static void* operator new(size_t);
    static void operator delete(void*);

  private:
    G4ThreeVector positionSensor;
    G4double energySensor;
    G4int SensorNumberID;
    G4String VolumeNameID; 

  public:
    inline void SetSensorPosition(G4ThreeVector position) { positionSensor = position; };
    inline G4ThreeVector GetSensorPosition() const { return positionSensor; };

    inline void SetSensorNumber(G4int SensorNumber) { SensorNumberID = SensorNumber; }; 
    inline G4int GetSensorNumber() const { return SensorNumberID; };

    inline void SetVolumeName(G4String VolumeName) { VolumeNameID = VolumeName; }; 
    inline G4String GetVolumeName() const { return VolumeNameID; }; 

 

    inline void SetSensorEnergy(G4double energy) { energySensor = energy; };
    inline G4double GetSensorEnergy() const { return energySensor; };
};

// Vector collection of one type of hits
typedef G4THitsCollection<SensorHit> SensorHitsCollection;

// Allocator declaration
extern G4Allocator<SensorHit> SensorHitsAllocator;

// Define operator new inline (already done in the header)
inline void* SensorHit::operator new(size_t) {
  return (void*)SensorHitsAllocator.MallocSingle();
}

// Define operator delete inline (already done in the header)
inline void SensorHit::operator delete(void* aSensorHit) {
  SensorHitsAllocator.FreeSingle((SensorHit*)aSensorHit);
}

#endif
