#include "FAZIAPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "TMath.h"
#include <iostream>
//#define restM 105 //MeV

using namespace CLHEP;
using namespace std;

  FAZIAPrimaryGeneratorAction::FAZIAPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
}

FAZIAPrimaryGeneratorAction::~FAZIAPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void FAZIAPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  costheta = RandGauss::shoot(0,0.05); //mean , std 
  random = G4UniformRand();  
  phirand = G4UniformRand();
  Parti_ene = G4UniformRand();
  //cout <<"rand " <<random << endl;
  //cout <<"partene " <<Parti_ene << endl;
  //cout << "" << endl;
  ///////////////////////////////////////////////////////////////////////////
  G4double Prototflux = 3339.4;
  G4double Deutotflux = 3339.4;
  G4double He3totflux = 3339.4;
  G4double Alptotflux = 3339.4;

  G4double E_avg[84];
  G4double Efor[84];
  G4double Elat[84];

  G4double fluxPro[84];
  G4double normfluxPro[84];
  G4double fluxDeu[84];
  G4double normfluxDeu[84];
  G4double fluxHe3[84];
  G4double normfluxHe3[84];
  G4double fluxAlp[84];
  G4double normfluxAlp[84];

  for(G4int i = 0 ; i < 84 ; i++)
  {	
	//file >> p0 >> p1 >> mup  >> mum;

	if(i==0)
	{
	  Efor[i] = PreEnergy[i];
	  Elat[i] = LatEnergy[i];

	  E_avg[i] = (Efor[i] + Elat[i])/2;

	  fluxPro[i] = probPro[i]*E_avg[i];
	  fluxDeu[i] = probDeu[i]*E_avg[i];
	  fluxHe3[i] = probHe3[i]*E_avg[i];
	  fluxAlp[i] = probAlp[i]*E_avg[i];
	}
	else
	{
	  Efor[i] = PreEnergy[i];
	  Elat[i] = LatEnergy[i];

	  E_avg[i] = (Efor[i] + Elat[i])/2;

	  fluxPro[i] = fluxPro[i-1] + probPro[i]*E_avg[i];
	  fluxDeu[i] = fluxDeu[i-1] + probDeu[i]*E_avg[i];
	  fluxHe3[i] = fluxHe3[i-1] + probHe3[i]*E_avg[i];
	  fluxAlp[i] = fluxAlp[i-1] + probAlp[i]*E_avg[i];
	}

	//E_avg[i] = (Efor[i] + Elat[i])/2;

	normfluxPro[i] = fluxPro[i]/148653.0;
	normfluxDeu[i] = fluxDeu[i]/148653.0;	
	normfluxHe3[i] = fluxHe3[i]/148653.0;
	normfluxAlp[i] = fluxAlp[i]/148653.0;	

	//cout << normfluxPro[i] << endl;
	//cout << normfluxDeu[i] << endl;
	//cout << normfluxHe3[i] << endl;
	//cout << normfluxAlp[i] << endl;
  }
  //////////////////////////////////////////////////////////////////////////////
  G4int n_particle = 1;

  G4double Phi_rand = 2*pi*phirand;
  //G4double Theta_rand = (60*CLHEP::degree)*((TMath::Cos(90*CLHEP::degree*costheta))*(TMath::Cos(90*CLHEP::degree*costheta)) );
  G4double Theta_rand = (pi/2)*costheta;
  G4int r_value = 20;

  G4double x = r_value*(TMath::Sin(Theta_rand))*(TMath::Cos(Phi_rand))*m;
  G4double y = r_value*(TMath::Sin(Theta_rand))*(TMath::Sin(Phi_rand))*m;
  G4double z = r_value*(TMath::Cos(Theta_rand))*m;

  fParticleGun = new G4ParticleGun(n_particle);
  G4ThreeVector gunposition = G4ThreeVector(0,0,0);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  G4ParticleDefinition* particle; 

  for(G4int j = 0 ; j < 84 ; j++)
  {
	//cout << random << " and " << fluxM[j]/fluxP[j] << endl;
	if(0. <= random && random < Prototflux/13357.600)//select particle sort
	{
	  if(normfluxPro[j] <= Parti_ene && Parti_ene < normfluxPro[j+1])//select particle Energy
	  {
		particle = particleTable -> FindParticle(particleName = "proton");		

		fParticleGun -> SetParticlePosition(gunposition);  
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
		fParticleGun -> SetParticleEnergy(E_avg[j]*MeV);
		fParticleGun -> GeneratePrimaryVertex(anEvent);
		std::cout <<"event" << eventNum++ << G4endl;
		std::cout <<"particle" << particleName << G4endl;
		std::cout <<"Energy" << E_avg[j] << G4endl;
	  }
	  else {
		continue;
	  }
	}

	if(Prototflux/13357.600 <= random && random <= (Prototflux + Deutotflux)/13357.600 )
	{
	  if(normfluxDeu[j] <= Parti_ene && Parti_ene < normfluxDeu[j+1])
	  {
		particle = particleTable -> FindParticle(particleName = "deuteron");		

		fParticleGun -> SetParticlePosition(gunposition);  
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
		fParticleGun -> SetParticleEnergy(E_avg[j]*MeV);
		fParticleGun -> GeneratePrimaryVertex(anEvent);
		std::cout <<"event" << eventNum++ << G4endl;
		std::cout <<"particle" << particleName << G4endl;
		std::cout <<"Energy" << E_avg[j] << G4endl;
	  }
	  else {
		continue;
	  }
	}

	if((Prototflux + Deutotflux)/13357.600  <= random && random < (Prototflux + Deutotflux + He3totflux)/13357.600)
	{
	  if(normfluxHe3[j] <= Parti_ene && Parti_ene < normfluxHe3[j+1])
	  {
		particle = particleTable -> FindParticle(particleName = "He3");		

		fParticleGun -> SetParticlePosition(gunposition);  
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
		fParticleGun -> SetParticleEnergy(E_avg[j]*MeV);
		fParticleGun -> GeneratePrimaryVertex(anEvent);
		std::cout <<"event" << eventNum++ << G4endl;
		std::cout <<"particle" << particleName << G4endl;
		std::cout <<"Energy" << E_avg[j] << G4endl;
	  }
	  else {
		continue;
	  }
	}

	if( (Prototflux + Deutotflux + He3totflux)/13357.600 <= random && random <= (Prototflux + Deutotflux + He3totflux + Alptotflux)/13357.600 )
	{
	  if(normfluxAlp[j] <= Parti_ene && Parti_ene < normfluxAlp[j+1])
	  {
		particle = particleTable -> FindParticle(particleName = "alpha");		

		fParticleGun -> SetParticlePosition(gunposition);  
		fParticleGun -> SetParticleDefinition(particle);
		fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(x,y,z));
		fParticleGun -> SetParticleEnergy(E_avg[j]*MeV);
		fParticleGun -> GeneratePrimaryVertex(anEvent);
		std::cout <<"event" << eventNum++ << G4endl;
		std::cout <<"particle" << particleName << G4endl;
		std::cout <<"Energy" << E_avg[j] << G4endl;
	  }
	  else {
		continue;
	  }
	}

  }
  //this function is called at the begining of each event

  fParticleGun -> GeneratePrimaryVertex(anEvent);
  //std::cout << eventNum++ << G4endl;
}

