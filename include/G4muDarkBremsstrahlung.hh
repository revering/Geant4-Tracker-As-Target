#ifndef G4muDarkBremsstrahlung_h
#define G4muDarkBremsstrahlung_h

#include "G4VEnergyLossProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4Material;

class G4muDarkBremsstrahlung : public G4VEnergyLossProcess
{

   public:

      G4muDarkBremsstrahlung(const G4String& name = "eDBrem");

      virtual ~G4muDarkBremsstrahlung();

      virtual G4bool IsApplicable(const G4ParticleDefinition& p);

      virtual void PrintInfo();

   protected:

      virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                               const G4ParticleDefinition*);
      G4bool isInitialised;

   private:

      G4muDarkBremsstrahlung & operator=(const G4muDarkBremsstrahlung &right);
      G4muDarkBremsstrahlung(const G4muDarkBremsstrahlung&);

};


#endif
