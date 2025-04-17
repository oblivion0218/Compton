#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

int main(int argc, char** argv) {
    // Inizializzazione dell'interfaccia utente
    G4UIExecutive* ui = nullptr;
    if (argc == 1) {
        ui = new G4UIExecutive(argc, argv);
    }

    // Creazione del RunManager
    auto* runManager = new G4RunManager();

    // Imposta le classi utente personalizzate
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserAction(new PrimaryGeneratorAction());
    runManager->SetUserAction(new RunAction());
    runManager->SetUserAction(new EventAction());

    // Inizializzazione Geant4
    runManager->Initialize();

    // Avvio della visualizzazione
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();

    // Avvio interfaccia utente o esecuzione batch
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    if (!ui) {
        UImanager->ApplyCommand("/control/execute init.mac");
    } else {
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    }

    // Cleanup
    delete visManager;
    delete runManager;

    return 0;
}
