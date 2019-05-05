#include "FAZIAEventAction.hh"
#include "G4RunManager.hh"

  FAZIAEventAction::FAZIAEventAction()
: G4UserEventAction()
{
}

FAZIAEventAction::~FAZIAEventAction()
{
}

void FAZIAEventAction::BeginOfEventAction(const G4Event*)
{
  eventID = 0;
  //volumeID = 0;
  edepTot[5145] = 0.;
}

void FAZIAEventAction::EndOfEventAction(const G4Event*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4int eventID = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
  //G4int volumeID = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo(); 


  analysisManager -> FillNtupleDColumn(1, 0, eventID);
  //analysisManager -> FillNtupleDColumn(1, 1, volumeID);
  analysisManager -> FillNtupleDColumn(1, 1, edepTot[5145]);
  analysisManager -> AddNtupleRow(1);
}
