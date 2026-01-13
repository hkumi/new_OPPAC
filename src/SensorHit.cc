#include "SensorHit.hh"
#include "G4UnitsTable.hh"

#include <iomanip>

// Define the allocator with the correct name
G4Allocator<SensorHit> SensorHitsAllocator;

SensorHit::SensorHit()
    : positionSensor(G4ThreeVector(0., 0., 0.)) { }

SensorHit::~SensorHit() { }

SensorHit::SensorHit(const SensorHit& right) 
    : positionSensor(right.positionSensor) { }

const SensorHit& SensorHit::operator=(const SensorHit& right) {
    positionSensor = right.positionSensor;
    return *this;
}

int SensorHit::operator==(const SensorHit& right) const {
    return 0;
}

