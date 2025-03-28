#include "EventAction.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction() : G4UserEventAction(), fComptonCount(0), fTotalEdep(0.) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->CreateNtuple("ComptonData", "Scattering Compton");
    analysisManager->CreateNtupleIColumn("NumCompton");  // Numero di eventi Compton
    analysisManager->CreateNtupleDColumn("Edep");        // Energia depositata totale
    analysisManager->FinishNtuple();
}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*) {
    fComptonCount = 0;
    fTotalEdep = 0.;
}

void EventAction::EndOfEventAction(const G4Event*) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(0, fComptonCount);
    analysisManager->FillNtupleDColumn(1, fTotalEdep / keV);
    analysisManager->AddNtupleRow();
}

void EventAction::AddComptonEvent(G4double edep) {
    fComptonCount++;
    fTotalEdep += edep;
}

