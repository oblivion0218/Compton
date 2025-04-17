#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
    // Inizializza il generatore con un solo fotone per evento
    fParticleGun = new G4ParticleGun(1);

    // Definizione del tipo di particella: fotone
    G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);

    // Imposta l'energia a 511 keV
    fParticleGun->SetParticleEnergy(511*keV);

    // Imposta la posizione di emissione (origine)
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));

    // Imposta la direzione (verso il target, lungo l'asse z positivo)
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
