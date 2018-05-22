#ifndef G4eDarkBremsstrahlung_h
#define G4eDarkBremsstrahlung_h

#include "G4VEnergyLossProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4Material;

class G4eDarkBremsstrahlung : public G4VEnergyLossProcess
{

   public:

      G4eDarkBremsstrahlung(const G4String& name = "eDBrem");

      virtual ~G4eDarkBremsstrahlung();

      virtual G4bool IsApplicable(const G4ParticleDefinition& p);

      virtual void PrintInfo();

   protected:

      virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                               const G4ParticleDefinition*);
      G4bool isInitialised;

   private:

      G4eDarkBremsstrahlung & operator=(const G4eDarkBremsstrahlung &right);
      G4eDarkBremsstrahlung(const G4eDarkBremsstrahlung&);

};


#endif
