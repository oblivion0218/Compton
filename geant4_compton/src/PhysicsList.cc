#include "PhysicsList.hh"

// Includi la fisica elettromagnetica standard
#include "G4EmStandardPhysics.hh"

#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList() {
    // Imposta il valore di cut predefinito (ad esempio, 1 mm)
    defaultCutValue = 1.0 * mm;
    SetVerboseLevel(1);

    // Registra la fisica elettromagnetica: include scattering, fotoelettrico, Compton, etc.
    RegisterPhysics(new G4EmStandardPhysics());
}

PhysicsList::~PhysicsList() {}
