#ifndef G4DMEventAction_h
#define G4DMEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4DMRunAction.hh"
#include "globals.hh"


class G4DMEventAction : public G4UserEventAction
{
   public:
      G4DMEventAction(G4DMRunAction* runAction);
      virtual ~G4DMEventAction();
 
      virtual void BeginOfEventAction(const G4Event* event);
      virtual void EndOfEventAction(const G4Event* event);

      void AddEdep(G4double edep) { fEdep += edep; }
      G4DMRunAction* GetRunAction() const;
      G4int GetnBrem() { return nBrem; }
      void incnBrem() { nBrem++;}
      void rstnBrem() { nBrem = 0;}

   private:
      G4DMRunAction* fRunAction;
      G4double fEdep;
      G4int nBrem;
};

#endif

