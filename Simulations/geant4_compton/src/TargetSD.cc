#include "TargetSD.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4SystemOfUnits.hh"

TargetSD::TargetSD(const G4String& name) 
    : G4VSensitiveDetector(name) {}

TargetSD::~TargetSD() {}

G4bool TargetSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    // Ottenere il processo che ha generato questo step
    const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();

    // Verifica se il processo è uno scattering Compton ("compt")
    if (process && process->GetProcessName() == "compt") {
        G4double edep = aStep->GetTotalEnergyDeposit();

        if (edep > 0.) {  // Solo se c'è energia depositata
            // Recupera l'EventAction per registrare l'evento Compton
            EventAction* eventAction = (EventAction*)G4EventManager::GetEventManager()->GetUserEventAction();
            if (eventAction) {
                eventAction->AddComptonEvent(edep);
            }
        }
    }
    return true;
}

