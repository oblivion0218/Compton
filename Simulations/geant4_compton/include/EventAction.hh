#ifndef EventAction_hh
#define EventAction_hh

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class EventAction : public G4UserEventAction {
public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void AddComptonEvent(G4double edep);

private:
    G4int fComptonCount;
    G4double fTotalEdep;
};

#endif

