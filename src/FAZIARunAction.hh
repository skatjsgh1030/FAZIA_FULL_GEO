#ifndef FAZIARUNACTION_HH
#define FAZIARUNACTION_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "globals.hh"
#include "g4root.hh"

class FAZIARunAction : public G4UserRunAction
{
  public:
    FAZIARunAction();
    virtual ~FAZIARunAction();

    // method from the base class
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
};

#endif
