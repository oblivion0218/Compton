#include "RunAction.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

RunAction::RunAction() : G4UserRunAction() {
    // Creazione del gestore di analisi ROOT
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFileName("output"); // Salverà i dati in output.root
}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*) {
    // Apertura del file ROOT
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run*) {
    // Scrittura e chiusura del file ROOT
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    
    G4cout << "Run terminata, dati salvati in output.root" << G4endl;
}
