# Imposta la versione minima di CMake richiesta
cmake_minimum_required(VERSION 3.10)

# Imposta il nome del progetto
project(geant4_compton)

# Trova i pacchetti necessari di Geant4
find_package(Geant4 REQUIRED)

# Aggiungi le cartelle di inclusione di Geant4
include_directories(${Geant4_INCLUDE_DIRS})

# Definisci le cartelle dei file di sorgente e di inclusione
include_directories(${PROJECT_SOURCE_DIR}/include)

# Aggiungi i file di sorgente
set(SOURCE_FILES
    src/DetectorConstruction.cc
    src/PhysicsList.cc
    src/PrimaryGeneratorAction.cc
    src/RunAction.cc
    src/EventAction.cc
    src/TargetSD.cc
    main.cc  # Aggiungi il file main.cc
)

# Crea l'eseguibile
add_executable(geant4_compton ${SOURCE_FILES})

# Linka il progetto con le librerie Geant4
target_link_libraries(geant4_compton ${Geant4_LIBRARIES})

# Aggiungi il supporto per la visualizzazione di Geant4 (se necessario)
# Trova i pacchetti di visualizzazione Geant4
find_package(Geant4 REQUIRED ui_all vis_all)

# Linka anche le librerie di visualizzazione (se usate)
target_link_libraries(geant4_compton ${Geant4_LIBRARIES})



