#include "FAZIARunAction.hh"
#include "g4root.hh"
//#include "g4cvs.hh"

FAZIARunAction::FAZIARunAction()
: G4UserRunAction()
{
}

FAZIARunAction::~FAZIARunAction()
{
  delete G4AnalysisManager::Instance();

}

void FAZIARunAction::BeginOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> OpenFile("FAZIASimdata.root");

  analysisManager -> CreateNtuple("step", "step");
  analysisManager -> CreateNtupleIColumn("eventID");
  analysisManager -> CreateNtupleIColumn("volumeID");
  analysisManager -> CreateNtupleDColumn("edep");
  analysisManager -> CreateNtupleDColumn("x");
  analysisManager -> CreateNtupleDColumn("y");
  analysisManager -> CreateNtupleDColumn("z");
  analysisManager -> FinishNtuple();

  analysisManager -> CreateNtuple("edep","edep");
  analysisManager -> CreateNtupleDColumn("eventID");
  //analysisManager -> CreateNtupleDColumn("volumeID");
  analysisManager -> CreateNtupleDColumn("edepTot");
  analysisManager -> FinishNtuple();
  
}

void FAZIARunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager -> Write();
  analysisManager -> CloseFile();
}
