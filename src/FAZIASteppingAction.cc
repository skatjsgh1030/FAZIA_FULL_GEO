#include "FAZIASteppingAction.hh"
#include "FAZIAEventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4Step.hh"

#include "G4ThreeVector.hh"

#include "G4EventManager.hh"
//#include "G4cvs.hh"
//#include "G4root.hh"

  FAZIASteppingAction::FAZIASteppingAction()
: G4UserSteppingAction()
{
  //fScintillationCounter = 0;
  //fCerenkovCounter      = 0;
  //fEventNumber = 0;
}

FAZIASteppingAction::~FAZIASteppingAction()
{
}

void FAZIASteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();

  /*G4int eventNumber = G4RunManager::GetRunManager()->
	GetCurrentEvent()->GetEventID();

	if (eventNumber != fEventNumber) {
	fEventNumber = eventNumber;
	fScintillationCounter = 0;
  //fCerenkovCounter = 0;
  }

  G4Track* track = step->GetTrack();

  G4String ParticleName = track->GetDynamicParticle()->
  GetParticleDefinition()->GetParticleName();

  if (ParticleName == "opticalphoton") return;

  const std::vector<const G4Track*>* secondaries =
  step->GetSecondaryInCurrentStep();

  if (secondaries -> size() > 0) 
  {
  for(unsigned int i = 0 ; i < secondaries->size() ; ++i) 
  {
  if (secondaries -> at(i) -> GetParentID() > 0) 
  {
  if(secondaries -> at(i) -> GetDynamicParticle() -> GetParticleDefinition() 
  == G4OpticalPhoton::OpticalPhotonDefinition())
  {

  if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Scintillation")
  fScintillationCounter++;

  //if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
  //== "Cerenkov")fCerenkovCounter++;
  }
  }
  }
  }*/
  G4int volumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();
  G4double edep = step -> GetTotalEnergyDeposit();
  G4double x = step -> GetDeltaPosition().x();  	
  G4double y = step -> GetDeltaPosition().y();  	
  G4double z = step -> GetDeltaPosition().z();

  //G4double moment = step -> GetDeltaMomentum();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> FillNtupleIColumn(0, eventID);
  analysisManager -> FillNtupleIColumn(1, volumeID);
  analysisManager -> FillNtupleDColumn(2, edep);
  analysisManager -> FillNtupleDColumn(3, x);
  analysisManager -> FillNtupleDColumn(4, y);
  analysisManager -> FillNtupleDColumn(5, z);
  //analysisManager -> FillNtupleDColumn(6, moment);
  analysisManager -> AddNtupleRow();

  FAZIAEventAction *eventAction[5145];
  for(int volI = 0 ; volI < 5145 ; volI++)
  {
	eventAction[volI] = (FAZIAEventAction *) G4EventManager::GetEventManager() -> GetUserEventAction();
	if (volumeID == volI)
	{
	  eventAction[volI] -> AddEnergyDeposit1(edep);
	}

  }
}
