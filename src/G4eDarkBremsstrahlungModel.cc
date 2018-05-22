#include "G4eDarkBremsstrahlungModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4DarkPhoton.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ModifiedTsai.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

G4eDarkBremsstrahlungModel::G4eDarkBremsstrahlungModel(const G4ParticleDefinition* p,
                                                       const G4String& nam)
   : G4VEmModel(nam),
     particle(0),
     isElectron(true),
     probsup(1.0),
     isInitialised(false)
{
   if(p) { SetParticle(p); }
   theDarkPhoton = G4DarkPhoton::DarkPhoton();
   MA = G4DarkPhoton::DarkPhoton()->GetPDGMass()/CLHEP::GeV;
   SetAngularDistribution(new G4ModifiedTsai());
   highKinEnergy = HighEnergyLimit();
   lowKinEnergy = LowEnergyLimit();
   fParticleChange = 0;
}

G4eDarkBremsstrahlungModel::~G4eDarkBremsstrahlungModel()
{
   size_t n = partialSumSigma.size();
   if(n>0)
   {
      for(size_t i=0; i<n; i++)
      {
         delete partialSumSigma[i];
      }
   }
}

void G4eDarkBremsstrahlungModel::SetParticle(const G4ParticleDefinition* p)
{
   particle = p;
   if(p == G4Electron::Electron())
   {
      isElectron = true; 
   }
   else
   {
      isElectron = false;
   }
}

void G4eDarkBremsstrahlungModel::Initialise(const G4ParticleDefinition* p,
                                            const G4DataVector& cuts)
{
   if(p) 
   { 
      SetParticle(p);
   }

   highKinEnergy = HighEnergyLimit();
   lowKinEnergy = LowEnergyLimit();
   const G4ProductionCutsTable* theCoupleTable=G4ProductionCutsTable::GetProductionCutsTable();
   
   if(theCoupleTable)
   {
      G4int numOfCouples = theCoupleTable->GetTableSize();
      G4int nn = partialSumSigma.size();
      G4int nc = cuts.size();
      if(nn > 0)
      {
      for (G4int ii=0; ii<nn; ii++)
         {
            G4DataVector* a=partialSumSigma[ii];
            if ( a )  delete a;
         }
         partialSumSigma.clear();
      }
       if(numOfCouples>0) 
       {
          for (G4int i=0; i<numOfCouples; i++) 
	  {
             G4double cute   = DBL_MAX;
             if(i < nc) cute = cuts[i];
             const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
             const G4Material* material = couple->GetMaterial();
             G4DataVector* dv = ComputePartialSumSigma(material, 0.5*highKinEnergy,
             std::min(cute, 0.25*highKinEnergy));
             partialSumSigma.push_back(dv);
          }
       }
   }

   if(isInitialised) return;
   fParticleChange = GetParticleChangeForLoss();
   isInitialised = true;
}

/*
G4double G4eDarkBremsstrahlungModel::ComputeDEDXPerVolume(const G4Material*,
                                                      const G4ParticleDefinition*,
						      G4double kineticEnergy,
						      G4double cutEnergy)
{
   return 0.0;
}

G4double G4eDarkBremsstrahlungModel::CrossSectionPerVolume(const G4Material* material,
                                                           const G4ParticleDefinition* p,
							   G4double ekin,
							   G4double emin,
							   G4double emax)
{
   SetupForMaterial(p, material, ekin);
   G4double cross = 0.0;
   const G4ElementVector* theElementVector = material->GetElementVector();
   const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
   G4int nelm = material->GetNumberOfElements();
   if(nelm > nsec)
   {
      xsec.resize(nelm);
      nsec = nelm;
   }
   for (G4int i=0; i<nelm; i++)
   {
      cross += theAtomNumDensityVector[i] * 
                ComputeCrossSectionPerAtom(p,(*theElementVector)[i],ekin,emin,emax);
      xsec[i] = cross;
   }
   return cross;
}
*/

G4double G4eDarkBremsstrahlungModel::DsigmaDx (double x, void * pp)
{
   ParamsForChi* params = (ParamsForChi*)pp;

   G4double Mel = 5.1E-04;

   G4double beta = sqrt(1 - (params->MMA)*(params->MMA)/(params->EE0)/(params->EE0));
   G4double num = 1.-x+x*x/3.;
   G4double denom = (params->MMA)*(params->MMA)*(1.-x)/x+Mel*Mel*x;
   G4double DsDx = beta*num/denom;

   return DsDx;
}

G4double G4eDarkBremsstrahlungModel::chi (double t, void * pp) 
{
   ParamsForChi* params = (ParamsForChi*)pp;
/* Reminder II:
 * params->AA;
 * params->ZZ;
 * params->MMA;
 * params->EE0;
 */
   G4double Mel = 5.1E-04;
   G4double MUp = 2.79;
   G4double Mpr = 0.938;

   G4double d = 0.164/pow((params->AA),2./3.);
   G4double ap = 773.0/Mel/pow((params->ZZ),2./3.);
   G4double a = 111.0/Mel/pow((params->ZZ),1./3.);
   G4double G2el = (params->ZZ)*(params->ZZ)*a*a*a*a*t*t/(1.0+a*a*t)/(1.0+a*a*t)/(1.0+t/d)/(1.0+t/d);
   G4double G2in = (params->ZZ)*ap*ap*ap*ap*t*t/(1.0+ap*ap*t)/(1.0+ap*ap*t)/(1.0+t/0.71)/(1.0+t/0.71)
    /(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)
    *(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr)*(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr);
   G4double G2 = G2el+G2in;
   G4double ttmin = (params->MMA)*(params->MMA)*(params->MMA)*(params->MMA)/4.0/(params->EE0)/(params->EE0);
//   G4double ttmin = lowerLimit(x,theta,p);
   G4double Under = G2*(t-ttmin)/t/t;
//   G4cout << "Under: " << Under << " MMA: " << params->MMA << " EEO: " << params->EE0 << G4endl;

   return Under;
}

G4double G4eDarkBremsstrahlungModel::ComputeCrossSectionPerAtom(
                                            const G4ParticleDefinition*,
                                            G4double E0, 
                                            G4double Z,   G4double A,
                                            G4double cut, G4double)
// Calculates the cross section per atom in GeV, then changes back to GEANT4 internal units.
{
   G4double cross = 0.0 ;
   if ( E0 < keV || E0 < cut) 
   {
      return cross; 
   }

   E0 = E0 / CLHEP::GeV;
   G4double Mel = 5.1E-04;
   if(E0 < 2.*MA) return 0.;

   //begin: chi-formfactor calculation
   gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
   
   G4double result, error;
   G4double tmin = MA*MA*MA*MA/(4.*E0*E0);
   G4double tmax = MA*MA;
    
   gsl_function F;
   ParamsForChi alpha = {1.0, 1.0, 1.0, 1.0};
   F.function = &G4eDarkBremsstrahlungModel::chi;
   F.params = &alpha;
   alpha.AA = A;
   alpha.ZZ = Z;
   alpha.MMA = MA;
   alpha.EE0 = E0;
    
   gsl_integration_qags (&F, tmin, tmax, 0, 1e-7, 1000,
   w, &result, &error);
   //printf ("chi/Z^2 = % .18f\n", result/(ZNucl*ZNucl));
   //printf ("result    = % .18f\n", result);
   //printf ("estimated error = % .18f\n", error);
   //printf ("intervals =  %d\n", w->size);
   
   G4double ChiRes = result;
   gsl_integration_workspace_free (w);
   
   gsl_integration_workspace * dxspace
      = gsl_integration_workspace_alloc (1000);
   gsl_function G;
   G.function = &DsigmaDx;
   G.params = &alpha;
   G4double xmin = 0;
   G4double xmax = 1;
   if((Mel/E0)>(MA/E0)) xmax = 1-Mel/E0;
   else xmax = 1-MA/E0;
   G4double res, err;

   gsl_integration_qags (&G, xmin, xmax, 0, 1e-7, 1000,
                         dxspace, &res, &err);

   G4double DsDx = res;
   gsl_integration_workspace_free(dxspace);

   G4double GeVtoPb = 3.894E08;   
   G4double alphaEW = 1.0/137.0;
   G4double epsilBench = 0.0001;
 
   cross= GeVtoPb*4.*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*ChiRes*DsDx*CLHEP::picobarn;
   if(cross < 0.) 
   { 
      cross = 0.; 
   }
//   G4cout << "Kinetic Energy: " << E0 << " Cross section: " << cross << G4endl;
   E0 = E0*CLHEP::GeV;
   return cross*1E14;
}

G4DataVector* G4eDarkBremsstrahlungModel::ComputePartialSumSigma(
                                        const G4Material* material,
                                        G4double kineticEnergy,
                                        G4double cut)

// Build the table of cross section per element. 
//The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material.
{
   G4int nElements = material->GetNumberOfElements();
   const G4ElementVector* theElementVector = material->GetElementVector();
   const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();
   G4DataVector* dv = new G4DataVector();
   G4double cross = 0.0;

   for (G4int i=0; i<nElements; i++ ) 
   {
      cross += theAtomNumDensityVector[i] * ComputeCrossSectionPerAtom( particle,
      kineticEnergy, (*theElementVector)[i]->GetZ(), (*theElementVector)[i]->GetA(), cut);
      dv->push_back(cross);
   }
   return dv;
}


void G4eDarkBremsstrahlungModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp, 
                                                const G4MaterialCutsCouple* couple,
                                                const G4DynamicParticle* dp,
                                                G4double tmin,
                                                G4double maxEnergy)
{
   G4double E0 = dp->GetKineticEnergy();
   G4double tmax = min(maxEnergy, E0);
   if(tmin >= tmax) { return; }
   // limits of the energy sampling
   E0  = E0 + electron_mass_c2;
   G4double Mel = 5.1E-04;
   E0 = E0 / 1000.;

   G4double Xmin = MA/E0;
   G4double Xmax = 0.998;
   if(Xmin>Xmax) { return; }
   double ThetaMaxA = 0.06;
   double integratedX = -log(1-Xmax);
   double XAcc, ThetaAcc, PhiAcc;
   double PhiEv = G4UniformRand() * 2. * 3.1415962;
   double ThetaConst = 100.;
   for(int iii=1;iii<10000;iii++)
   {
      double A = G4UniformRand()*integratedX;
      double XEv = -exp(-A)+1;
      double intTheta = 1./ThetaConst*log(1./ThetaConst/20.+ThetaMaxA);
      double intTzero = 1./ThetaConst*log(1./ThetaConst/20.);
      double B = G4UniformRand()*(intTzero-intTheta)+intTheta;
      double ThetaEv = -1./ThetaConst/20. + exp(ThetaConst*B);
      double smax = 16./(MA*MA*(1.-XEv))/(ThetaConst*ThetaEv+MA/E0)/E0/E0;
      double UU = G4UniformRand() * smax;

      double Uxtheta = E0*E0*ThetaEv*ThetaEv*XEv + MA*MA*(1.0-XEv)/XEv + Mel*Mel*XEv;
      double AA = (1. - XEv + XEv*XEv/2.) / (Uxtheta*Uxtheta);
      double BB = (1. - XEv)*(1. - XEv)*MA*MA/(Uxtheta*Uxtheta*Uxtheta*Uxtheta);
      double CC = MA*MA - Uxtheta*XEv/(1. - XEv);
      double sigma = ThetaEv * XEv * (AA + BB*CC);
      if(sigma > smax)  printf ("Maximum violated: ratio = % .18f, X: %e, Theta: %e\n", sigma/smax, XEv, ThetaEv);

      if(sigma >= UU) 
      {
         XAcc = XEv;
         ThetaAcc = ThetaEv;
         PhiAcc = PhiEv;
/*         fParticle.E0 = XAcc;
         fParticle.Theta = ThetaAcc;
         fParticle.Phi = PhiAcc;
*/	 break;
//         printf ("Accepted at iteration %d, X: %e, Theta: %e\n", iii, XEv, ThetaEv);
//       printf( "Ee = %e XAcc = %e ThetaAcc = %e\n ", E0, XAcc, ThetaAcc);
      }
   }
   
   if((ThetaAcc==0)&&(XAcc==0)&&(PhiAcc==0))
   {
      printf ("Did not manage to simulate.\n");
      return;
   }
   E0 = E0*1000.;
   G4double dmomentum = sqrt((XAcc*E0*XAcc*E0)-MA*MA*1000*1000);
   G4ThreeVector newDirection;
   newDirection.set(std::sin(ThetaAcc)*std::cos(PhiAcc),std::sin(ThetaAcc)*std::sin(PhiAcc), std::cos(ThetaAcc));
   newDirection.rotateUz(dp->GetMomentumDirection());
   newDirection.setMag(dmomentum);
   // create G4DynamicParticle object for the Gamma
   G4DynamicParticle* dphoton = new G4DynamicParticle(theDarkPhoton,newDirection);
   vdp->push_back(dphoton);
   G4double totMomentum = dp->GetMomentum().mag(); 
   G4double dpMomentum = dphoton->GetMomentum().mag(); 
  
   G4ThreeVector direction = (totMomentum*dp->GetMomentumDirection()
                              - dpMomentum*dphoton->GetMomentumDirection()).unit();
 
   // energy of primary
   G4double finalmomentum = (totMomentum*dp->GetMomentumDirection()-dpMomentum*dphoton->GetMomentumDirection()).mag();
   G4double finalKE = sqrt(finalmomentum*finalmomentum+electron_mass_c2*electron_mass_c2) - electron_mass_c2;
  
   // stop tracking and create new secondary instead of primary
   if(finalKE > SecondaryThreshold()) {
     fParticleChange->ProposeTrackStatus(fStopAndKill);
     fParticleChange->SetProposedKineticEnergy(0.0);
     G4DynamicParticle* el = 
       new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
                             direction, finalKE);
     vdp->push_back(el);
     // continue tracking
   } else {
     fParticleChange->SetProposedMomentumDirection(direction);
     fParticleChange->SetProposedKineticEnergy(finalKE);
   }
    
 } 
   
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
   
 const G4Element* G4eDarkBremsstrahlungModel::SelectRandomAtom(
            const G4MaterialCutsCouple* couple) 
 {
   // select randomly 1 element within the material
 
   const G4Material* material = couple->GetMaterial();
   G4int nElements = material->GetNumberOfElements();
   const G4ElementVector* theElementVector = material->GetElementVector();
 
   const G4Element* elm = 0;
 
   if(1 < nElements) {
 
     --nElements; 
     G4DataVector* dv = partialSumSigma[couple->GetIndex()];
     G4double rval = G4UniformRand()*((*dv)[nElements]);
 
     elm = (*theElementVector)[nElements];
     for (G4int i=0; i<nElements; ++i) {
       if (rval <= (*dv)[i]) {
         elm = (*theElementVector)[i];
         break;
       }
     }
   } else { elm = (*theElementVector)[0]; }
  
   SetCurrentElement(elm);
   return elm;
 }
 
 

