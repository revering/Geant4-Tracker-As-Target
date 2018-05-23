#ifndef TTargRun_h
#define TTargRun_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"

class TTargRun : public G4Run {
   public:
      TTargRun();
      virtual ~TTargRun();
      virtual void RecordEvent(const G4Event*);
   
   private:
      G4int nEvent;
      G4THitsMap<G4double> f_pt;
      G4THitsMap<G4double> i_pt;

   public:
      G4int get_nEvent(){return nEvent;}
      G4THitsMap<G4double> get_f_pt(){return f_pt;}
      G4THitsMap<G4double> get_i_pt(){return i_pt;}
};

#endif
