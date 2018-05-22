#ifndef TTargPrimaryGeneratorAction_h
#define TTargPrimaryGeneratorAction_h 

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class TTargPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
   public:
      TTargPrimaryGeneratorAction();
      virtual ~TTargPrimaryGeneratorAction();

      virtual void GeneratePrimaries(G4Event*);

      const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

   private:
      G4ParticleGun* fParticleGun;
};

#endif
