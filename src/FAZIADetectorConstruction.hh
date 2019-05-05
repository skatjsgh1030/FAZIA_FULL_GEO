#ifndef FAZIADETECTORCONSTRUCTION_HH
#define FAZIADETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class FAZIADetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    FAZIADetectorConstruction();
    virtual ~FAZIADetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
};

#endif
