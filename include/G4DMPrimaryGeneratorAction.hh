#ifndef G4DMPrimaryGeneratorAction_h
#define G4DMPrimaryGeneratorAction_h 

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class G4DMPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
   public:
      G4DMPrimaryGeneratorAction();
      virtual ~G4DMPrimaryGeneratorAction();

      virtual void GeneratePrimaries(G4Event*);

      const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

   private:
      G4ParticleGun* fParticleGun;
};

#endif
