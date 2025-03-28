#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "TargetSD.hh"

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction() {}
DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    // 1. Definizione del materiale
    G4NistManager* nist = G4NistManager::Instance();
    // Materiale del target: rame
    G4Material* copper = nist->FindOrBuildMaterial("G4_Cu");

    // 2. Definizione del volume World
    G4double world_sizeXY = 1.0*m;
    G4double world_sizeZ  = 1.0*m;
    G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");

    G4Box* solidWorld =
        new G4Box("World",                       // nome
                  0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);

    G4LogicalVolume* logicWorld =
        new G4LogicalVolume(solidWorld,            // solid
                            worldMat,              // materiale
                            "World");

    G4VPhysicalVolume* physWorld =
        new G4PVPlacement(0,                       // nessuna rotazione
                          G4ThreeVector(),         // posizione (origine)
                          logicWorld,              // volume logico
                          "World",                 // nome
                          0,                       // padre
                          false,                   // non duplicato
                          0,                       // copy number
                          true);                   // check overlaps

    // 3. Definizione del target: un cilindro di rame
    // Specifiche: spessore 1 cm (altezza) e raggio 2.5 cm
    G4double innerRadius = 0.*cm;
    G4double outerRadius = 2.5*cm;
    G4double height = 1.0*cm;
    G4double startAngle = 0.*deg;
    G4double spanningAngle = 360.*deg;

    G4Tubs* solidTarget =
        new G4Tubs("Target", innerRadius, outerRadius, 0.5*height, startAngle, spanningAngle);

    G4LogicalVolume* logicTarget =
        new G4LogicalVolume(solidTarget, copper, "Target");
    
    
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    TargetSD* targetSD = new TargetSD("TargetSD");
    sdManager->AddNewDetector(targetSD);
    logicTarget->SetSensitiveDetector(targetSD);
    
    // 4. Posizionamento del target: a 4 cm dalla sorgente lungo l'asse z
    // Supponiamo che la sorgente sia all'origine, posizioniamo il target a z = 4 cm.
    G4ThreeVector targetPosition = G4ThreeVector(0, 0, 4*cm);

    new G4PVPlacement(0,                      // nessuna rotazione
                      targetPosition,         // posizione
                      logicTarget,            // volume logico
                      "Target",               // nome
                      logicWorld,             // volume padre
                      false,                  // non duplicato
                      0,                      // copy number
                      true);                  // check overlaps

    // 5. (Opzionale) Definizione degli attributi visivi
    logicWorld->SetVisAttributes(new G4VisAttributes(false)); // rendi invisibile il world
    G4VisAttributes* visTarget = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // target in rosso
    logicTarget->SetVisAttributes(visTarget);

    return physWorld;
}
