#include "G4muDarkBremsstrahlungModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
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

G4muDarkBremsstrahlungModel::G4muDarkBremsstrahlungModel(const G4ParticleDefinition* p,
                                                       const G4String& nam)
   : G4VEmModel(nam),
     particle(0),
     isMuon(true),
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
   MuonMass = 105.658;
}

G4muDarkBremsstrahlungModel::~G4muDarkBremsstrahlungModel()
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

void G4muDarkBremsstrahlungModel::SetParticle(const G4ParticleDefinition* p)
{
   particle = p;
   if((p==G4MuonPlus::MuonPlus())||(p==G4MuonMinus::MuonMinus()))
   {
      isMuon = true; 
   }
   else
   {
      isMuon = false;
   }
}

void G4muDarkBremsstrahlungModel::Initialise(const G4ParticleDefinition* p,
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


/*G4double G4muDarkBremsstrahlungModel::ComputeDEDXPerVolume(
                                 const G4Material* material,
				 const G4ParticleDefinition* p,
				       G4double kineticEnergy,
				       G4double cutEnergy)
{
   return 0.0;
}
*/
/*
   if(!particle)
   {
      SetParticle(p);
   }
   if(kineticEnergy < lowKinEnergy)
   {
      return 0.0;
   }
   const G4double thigh = 100.*GeV;
 
   G4double cut = std::min(cutEnergy, kineticEnergy);
 
   G4double rate, loss;
   const G4double factorHigh = 36./(1450.*GeV);
   const G4double coef1 = -0.5;
   const G4double coef2 = 2./9.;
 
   const G4ElementVector* theElementVector = material->GetElementVector();
   const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();
 
   G4double totalEnergy = kineticEnergy + MuonMass;
   G4double dedx = 0.0;
 
   //  loop for elements in the material
   for (size_t i=0; i<material->GetNumberOfElements(); i++) 
   {
 
      G4double Z     = (*theElementVector)[i]->GetZ();
      G4double natom = theAtomicNumDensityVector[i];
 
     // loss for MinKinEnergy<KineticEnergy<=100 GeV
      if (kineticEnergy <= thigh) 
      { 
         // x = log(totalEnergy/MuonMass);
         loss = ComputeBremLoss(Z, kineticEnergy, cut) ;
         if (!isMuon) loss *= PositronCorrFactorLoss(Z, kineticEnergy, cut);
         // extrapolation for KineticEnergy>100 GeV
      } 
      else 
      {
         // G4double xhigh = log(thigh/MuonMass);
         G4double cuthigh = thigh*0.5;
 
         if (cut < thigh) 
	 {
            loss = ComputeBremLoss(Z, thigh, cut) ;
            if (!isMuon) loss *= PositronCorrFactorLoss(Z, thigh, cut) ;
            rate = cut/totalEnergy;
            loss *= (1. + coef1*rate + coef2*rate*rate);
            rate = cut/thigh;
            loss /= (1.+coef1*rate+coef2*rate*rate);
         } 
         else 
	 {
            loss = ComputeBremLoss(Z, thigh, cuthigh) ;
            if (!isMuon) loss *= PositronCorrFactorLoss(Z, thigh, cuthigh) ;
            rate = cut/totalEnergy;
            loss *= (1. + coef1*rate + coef2*rate*rate);
            loss *= cut*factorHigh;
         }
      }
      loss *= natom;
 
      G4double kp2 = MigdalConstant*totalEnergy*totalEnergy
                  * (material->GetElectronDensity()) ;
 
      // now compute the correction due to the supression(s)
      G4double kmin = 1.*eV;
      G4double kmax = cut;
 
      if (kmax > kmin) 
      {
         G4double floss = 0.;
         G4int nmax = 100;
 
         G4double vmin=log(kmin);
         G4double vmax=log(kmax) ;
         G4int nn = (G4int)(nmax*(vmax-vmin)/(log(highKinEnergy)-vmin)) ;
         G4double u,fac,c,v,dv ;
         if(nn > 0) 
         {
            dv = (vmax-vmin)/nn ;
            v  = vmin-dv ;
 
            for(G4int n=0; n<=nn; n++)
            {
               v += dv;
               u = exp(v);
               fac = u*SupressionFunction(material,kineticEnergy,u);
               fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;
               if ((n==0)||(n==nn)) c=0.5;
               else                 c=1. ;
               fac   *= c ;
               floss += fac ;
            }
         floss *=dv/(kmax-kmin);
         } 
	 else 
         {
         floss = 1.;
         }
         if(floss > 1.) floss = 1.;
         // correct the loss
         loss *= floss;
      } 
      dedx += loss;
   }
   if(dedx < 0.) { dedx = 0.; }
   return dedx;
}
*/

/*
   static const G4double beta=1.0, ksi=2.0;
   static const G4double clossh = 0.254 , closslow = 1./3. , alosslow = 1. ;
   static const G4double Tlim= 10.*MeV ;
 
   static const G4double xlim = 1.2 ;
   static const G4int NZ = 8 ;
   static const G4int Nloss = 11 ;
   static const G4double ZZ[NZ] =
          {2.,4.,6.,14.,26.,50.,82.,92.};
   static const G4double coefloss[NZ][Nloss] = {
    // Z=2
   { 0.98916,        0.47564,        -0.2505,       -0.45186,        0.14462,
     0.21307,      -0.013738,      -0.045689,     -0.0042914,      0.0034429,
     0.00064189},
 
   // Z=4
   { 1.0626,        0.37662,       -0.23646,       -0.45188,        0.14295,
     0.22906,      -0.011041,      -0.051398,     -0.0055123,      0.0039919,
     0.00078003},
    // Z=6
   { 1.0954,          0.315,       -0.24011,       -0.43849,        0.15017,
     0.23001,      -0.012846,      -0.052555,     -0.0055114,      0.0041283,
     0.00080318},
 
    // Z=14
   { 1.1649,        0.18976,       -0.24972,       -0.30124,         0.1555,
     0.13565,      -0.024765,      -0.027047,    -0.00059821,      0.0019373,
     0.00027647},
 
    // Z=26
   { 1.2261,        0.14272,       -0.25672,       -0.28407,        0.13874,
     0.13586,      -0.020562,      -0.026722,    -0.00089557,      0.0018665,
     0.00026981},
 
    // Z=50
   { 1.3147,       0.020049,       -0.35543,       -0.13927,        0.17666,
     0.073746,      -0.036076,      -0.013407,      0.0025727,     0.00084005,
    -1.4082e-05},
 
    // Z=82
   { 1.3986,       -0.10586,       -0.49187,     -0.0048846,        0.23621,
     0.031652,      -0.052938,     -0.0076639,      0.0048181,     0.00056486,
    -0.00011995},
 
    // Z=92
   { 1.4217,         -0.116,       -0.55497,      -0.044075,        0.27506,
     0.081364,      -0.058143,      -0.023402,      0.0031322,      0.0020201,
     0.00017519}
 
     } ;
   static G4double aaa = 0.414;
   static G4double bbb = 0.345;
   static G4double ccc = 0.460;
 
   G4int iz = 0;
   G4double delz = 1.e6;
   for (G4int ii=0; ii<NZ; ii++)
   {
      G4double dz = std::abs(Z-ZZ[ii]);
      if(dz < delz)
      {
         iz = ii;
         delz = dz;
      }
   }
 
   G4double xx = log10(T/MeV);
   G4double fl = 1.;
 
   if (xx <= xlim)
   {
      xx /= xlim;
      G4double yy = 1.0;
      fl = 0.0;
      for (G4int j=0; j<Nloss; j++) 
      {
         fl += yy+coefloss[iz][j];
         yy *= xx;
      }
      if (fl < 0.00001) fl = 0.00001;
      else if (fl > 1.0) fl = 1.0;
   }
 
   G4double loss;
   G4double E = T+electron_mass_c2 ;
 
   loss = Z*(Z+ksi)*E*E/(T+E)*exp(beta*log(Cut/T))*(2.-clossh*exp(log(Z)/4.));
   if (T <= Tlim) loss /= exp(closslow*log(Tlim/T));
   if( T <= Cut)  loss *= exp(alosslow*log(T/Cut));
   //  correction
   loss *= (aaa+bbb*T/Tlim)/(1.+ccc*T/Tlim);
   loss *= fl;
   loss /= Avogadro;
 
   return loss;
}
*/
/*
G4double G4muDarkBremsstrahlungModel::PositronCorrFactorLoss(G4double Z, G4double kineticEnergy, G4double cut)
{
   static const G4double K = 132.9416*eV ;
   static const G4double a1=4.15e-1, a3=2.10e-3, a5=54.0e-5 ;

   G4double x   = log(kineticEnergy/(K*Z*Z)), x2 = x*x, x3 = x2*x;
   G4double eta = 0.5+atan(a1*x+a3*x3+a5*x3*x2)/pi;
   G4double e0  = cut/kineticEnergy;

   G4double factor = 0.0;
   if (e0 < 1.0) 
   {
      factor=log(1.-e0)/eta;
      factor=exp(factor);
   }
   factor = eta*(1.-factor)/e0;

return factor;
}
*/
/*
G4double G4muDarkBremsstrahlungModel::CrossSectionPerVolume(
                                    const G4Material* material,
                                    const G4ParticleDefinition* p,
                                    G4double kineticEnergy,
                                    G4double cutEnergy,
                                    G4double maxEnergy)
{
   if(!particle) 
   {
      SetParticle(p);
   }
   G4double cross = 0.0;
   G4double tmax = min(maxEnergy, kineticEnergy);
   G4double cut  = min(cutEnergy, kineticEnergy);
   if(cut >= tmax) 
   {
      return cross; 
   }
   const G4ElementVector* theElementVector = material->GetElementVector();
   const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();
   for (size_t i=0; i<material->GetNumberOfElements(); i++) 
   {
       cross += theAtomNumDensityVector[i] * ComputeCrossSectionPerAtom(p,
       kineticEnergy, (*theElementVector)[i]->GetZ(), (*theElementVector)[i]->GetA(), cut);
       if (tmax < kineticEnergy) 
       {
          cross -= theAtomNumDensityVector[i] * ComputeCrossSectionPerAtom(p,
          kineticEnergy, (*theElementVector)[i]->GetZ(), (*theElementVector)[i]->GetA(), tmax);
       }
   }
   
   G4double kmax = tmax;
   G4double kmin = cut;
  
   G4double totalEnergy = kineticEnergy+MuonMass ;
   G4double kp2 = MigdalConstant*totalEnergy*totalEnergy
   *(material->GetElectronDensity());
  
   G4double fsig = 0.;
   G4int nmax = 100;
   G4double vmin=log(kmin);
   G4double vmax=log(kmax) ;
   G4int nn = (G4int)(nmax*(vmax-vmin)/(log(highKinEnergy)-vmin));
   G4double u,fac,c,v,dv,y ;
   if(nn > 0) 
   {
      dv = (vmax-vmin)/nn ;
      v  = vmin-dv ;
      for(G4int n=0; n<=nn; n++) 
      { 
         v += dv;  
         u = exp(v);              
         fac = SupressionFunction(material, kineticEnergy, u);
         y = u/kmax;
         fac *= (4.-4.*y+3.*y*y)/3.;
         fac *= probsup*(u*u/(u*u+kp2))+1.-probsup; 
         if ((n==0)||(n==nn)) c=0.5;
         else                 c=1. ;  
         fac  *= c;
         fsig += fac;
      }
      y = kmin/kmax ;
      fsig *=dv/(-4.*log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));
  
   }
   else 
   {  
      fsig = 1.;
   }
   if (fsig > 1.) fsig = 1.;
  
 // correct the cross section
   cross *= fsig;
  
   return cross;
}   
*/

G4double G4muDarkBremsstrahlungModel::DsigmaDx (double x, void * pp)
{
   ParamsForChi* params = (ParamsForChi*)pp;

   G4double MMu = 105.658/1000.;
   G4double beta = sqrt(1- (params->MMA)*(params->MMA)/(params->EE0)/(params->EE0));
   G4double num = 1.-x+x*x/3.;
   G4double denom = (params->MMA)*(params->MMA)*(1.-x)/x+MMu*MMu*x;
   G4double DsDx = beta*num/denom;

   return DsDx;
}

G4double G4muDarkBremsstrahlungModel::chi(double t, void * pp) 
{
   ParamsForChi* params = (ParamsForChi*)pp;
/* Reminder II:
 * params->AA;
 * params->ZZ;
 * params->MMA;
 * params->EE0;
 */
   G4double Me = 0.511/1000.;
   G4double MUp = 2.79;
   G4double Mpr = 0.938;

   G4double d = 0.164/pow((params->AA),2./3.);
   G4double ap = 773.0/Me/pow((params->ZZ),2./3.);
   G4double a = 111.0/Me/pow((params->ZZ),1./3.);
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

G4double G4muDarkBremsstrahlungModel::ComputeCrossSectionPerAtom(
                                            const G4ParticleDefinition*,
                                            G4double E0, 
                                            G4double Z,   G4double A,
                                            G4double cut, G4double)
// Calculates the cross section per atom in GEANT4 internal units.
{
   G4double cross = 0.0 ;
   if ( E0 < keV || E0 < cut) 
   {
      return cross; 
   }

   E0 = E0 / CLHEP::GeV;
   G4double Mmu = MuonMass/1000.;
   if(E0 < 2.*MA) return 0.;

   //begin: chi-formfactor calculation
   gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
   
   G4double result, error;
   G4double tmin = MA*MA*MA*MA/(4.*E0*E0);
   G4double tmax = MA*MA;
    
   gsl_function F;
   ParamsForChi alpha = {1.0, 1.0, 1.0, 1.0};
   F.function = &G4muDarkBremsstrahlungModel::chi;
   F.params = &alpha;
//   G4cout << "MA: " << MA << G4endl; 
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
   if((Mmu/E0)>(MA/E0)) xmax = 1-Mmu/E0;
   else xmax = 1-MA/E0;
   G4double res, err;

   gsl_integration_qags (&G, xmin, xmax, 0, 1e-7, 1000,
                         dxspace, &res, &err);

   G4double DsDx = res;
   gsl_integration_workspace_free(dxspace);      
      
   
   G4double GeVtoPb = 3.894E08;   
   G4double alphaEW = 1.0/137.0;
//   G4double epsilBench = 0.0001;
   G4double epsilBench = 1;

   cross= GeVtoPb*4.*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*ChiRes*DsDx*CLHEP::picobarn;
   if(cross < 0.) 
   { 
      cross = 0.; 
   }
//   G4cout << "A mass: " << MA << G4endl;
//   G4cout << "Kinetic Energy: " << E0 << " Cross section: " << cross << G4endl;
   E0 = E0*CLHEP::GeV;
//   G4cout << "Muon Energy is: " << E0/CLHEP::GeV << ". Cross section is: " << cross/CLHEP::picobarn << " picobarns.\n";
   return cross*1E5;
}

/*   G4double G4muDarkBremsstrahlungModel::PositronCorrFactorSigma( G4double Z,
   G4double kineticEnergy, G4double cut)

//Calculates the correction factor for the total cross section of the positron
//  bremsstrahlung.
// Eta is the ratio of positron to electron energy loss by bremstrahlung. 
// A parametrized formula from L. Urban is used to estimate eta. It is a fit to
// the results of L. Kim & al: Phys Rev. A33,3002 (1986)

{   
   static const G4double K = 132.9416*eV;
   static const G4double a1 = 4.15e-1, a3 = 2.10e-3, a5 = 54.0e-5;

   G4double x    = log(kineticEnergy/(K*Z*Z));
   G4double x2 = x*x;
   G4double x3 = x2*x;
   G4double eta  = 0.5 + atan(a1*x + a3*x3 + a5*x3*x2)/pi ;
   G4double alfa = (1. - eta)/eta;
   return eta*pow((1. - cut/kineticEnergy), alfa);
}
*/
G4DataVector* G4muDarkBremsstrahlungModel::ComputePartialSumSigma(
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


void G4muDarkBremsstrahlungModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp, 
                                                const G4MaterialCutsCouple* couple,
                                                const G4DynamicParticle* dp,
                                                G4double tmin,
                                                G4double maxEnergy)
{
   G4double E0 = dp->GetKineticEnergy();
   G4double tmax = min(maxEnergy, E0);
   if(tmin >= tmax) { return; }
   // limits of the energy sampling
   E0  = E0 + MuonMass;
   G4double Mmu = MuonMass/1000.;
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

      double Uxtheta = E0*E0*ThetaEv*ThetaEv*XEv + MA*MA*(1.0-XEv)/XEv + Mmu*Mmu*XEv;
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
   G4double finalKE = sqrt(finalmomentum*finalmomentum+MuonMass*MuonMass) - MuonMass;
  
   // stop tracking and create new secondary instead of primary
   if(finalKE > SecondaryThreshold()) {
     fParticleChange->ProposeTrackStatus(fStopAndKill);
     fParticleChange->SetProposedKineticEnergy(0.0);
     G4DynamicParticle* mu = 
       new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
                             direction, finalKE);
     vdp->push_back(mu);
     // continue tracking
   } else {
     fParticleChange->SetProposedMomentumDirection(direction);
     fParticleChange->SetProposedKineticEnergy(finalKE);
   }
    
 } 
   
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
   
 const G4Element* G4muDarkBremsstrahlungModel::SelectRandomAtom(
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
 
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
 G4double G4muDarkBremsstrahlungModel::SupressionFunction(const G4Material* material,
                                  G4double kineticEnergy, G4double gammaEnergy)
 {
   // supression due to the LPM effect+polarisation of the medium/
   // supression due to the polarisation alone
 
 
   G4double totEnergy = kineticEnergy+MuonMass ;
   G4double totEnergySquare = totEnergy*totEnergy ;
 
   G4double LPMEnergy = LPMconstant*(material->GetRadlen()) ;
 
   G4double gammaEnergySquare = gammaEnergy*gammaEnergy ;
 
   G4double muonDensity = 0;//material->GetmuonDensity();
 
   G4double sp = gammaEnergySquare/
    (gammaEnergySquare+MigdalConstant*totEnergySquare*muonDensity);
 
   G4double supr = 1.0;
 
   if (false) {
 
     G4double s2lpm = LPMEnergy*gammaEnergy/totEnergySquare;
 
     if (s2lpm < 1.) {
 
       G4double LPMgEnergyLimit = totEnergySquare/LPMEnergy ;
       G4double LPMgEnergyLimit2 = LPMgEnergyLimit*LPMgEnergyLimit;
       G4double splim = LPMgEnergyLimit2/
         (LPMgEnergyLimit2+MigdalConstant*totEnergySquare*muonDensity);
       G4double w = 1.+1./splim ;
 
       if ((1.-sp) < 1.e-6) w = s2lpm*(3.-sp);
       else                 w = s2lpm*(1.+1./sp);
 
       supr = (sqrt(w*w+4.*s2lpm)-w)/(sqrt(w*w+4.)-w) ;
       supr /= sp;    
     } 
   
   } 
   return supr;
}
*/
