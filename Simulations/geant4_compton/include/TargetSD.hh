#ifndef TargetSD_h
#define TargetSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class TargetSD : public G4VSensitiveDetector {
public:
    TargetSD(const G4String& name);
    virtual ~TargetSD();

    // Questo metodo verrà chiamato ad ogni hit nel volume target
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
};

#endif
