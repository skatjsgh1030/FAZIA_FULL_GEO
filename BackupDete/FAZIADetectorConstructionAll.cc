#include "FAZIADetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformElectricField.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iostream>
#include "TMath.h"
//COR
#define Blktheta1 0.10728646 //0.10728646
#define Blktheta2 0.10062866 //0.10162866
#define Blktheta3 0.08820 //0.08920
#define Blktheta4 0.0812759 //0.0822759
#define Blkphi1 1.1902899
#define Blkphi2 1.3734008
#define Blkphi3 1.1071487
#define Blkphi4 1.3258177
//MID
#define Blktheta41 0.15495 //0.15495
#define Blktheta42 0.140489 //
#define Blktheta43 0.143232 //
#define Blktheta44 0.1273692 //
#define Blkphi41 0.6947382
#define Blkphi42 0.78539816
#define Blkphi43 0.58800260
#define Blkphi44 0.67474094
//OUTER
#define Blktheta51 0.18234275 //
#define Blktheta52 0.17916567 //
#define Blktheta53 0.16345286 //
#define Blktheta54 0.15986910 //
#define Blkphi51 1.3521274
#define Blkphi52 1.4601391
#define Blkphi53 1.3258177
#define Blkphi54 1.4464413

using std::stringstream;
using namespace std;
using namespace CLHEP;
  FAZIADetectorConstruction::FAZIADetectorConstruction()
: G4VUserDetectorConstruction()
{
}

FAZIADetectorConstruction::~FAZIADetectorConstruction()
{
}

G4VPhysicalVolume* FAZIADetectorConstruction::Construct()
{
  //----------------Material --------------------
  G4NistManager* nist = G4NistManager::Instance();

  const G4double labTemp = CLHEP::STP_Temperature + 20.*kelvin;

  G4Element*  elCs  = new G4Element("Caesium","Cs",  55,   132.90545*g/mole);
  G4Element*  elI   = new G4Element("Iodine", "I" , 53,   126.90447*g/mole);

  G4Material* CsI = new G4Material("CsI", 4.51*g/CLHEP::cm3, 2, kStateSolid, labTemp);
  CsI -> AddElement(elCs, 1);
  CsI -> AddElement(elI,  1);

  // ------------ Generate & Add Material Properties Table ------------
  //
  G4double photonEnergy[] =
  {	1.90738*eV, 1.98368*eV, 2.06633*eV,
	2.15617*eV, 2.25418*eV, 2.36152*eV,
	2.47960*eV, 2.61011*eV, 2.75511*eV,
	2.91718*eV, 3.0995*eV};
  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  // CsI
  //fixed
  G4double refractiveIndex1[] =
  { 1.7789, 1.7820, 1.7856, 1.7897, 1.7945,
	1.8001, 1.8067, 1.8146, 1.8242, 1.8359,
	1.8506};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] =
  {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
	30.000*m, 28.500*m, 27.000*m,17.500*m, 14.500*m };

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] =
  { 1.00, 1.00, 1.00, 
	1.00, 1.00, 1.00, 
	1.00,1.00, 1.00,
	1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4double scintilSlow[] =
  { 0.01, 1.00, 2.00,
	3.00, 4.00, 5.00,
	6.00,7.00, 6.00,
	5.00, 4.00 };

  assert(sizeof(scintilSlow) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
	->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
	->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
	->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
	->SetSpline(true);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

  CsI->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator

  CsI->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // -----------------------------------------------------
  // World

  G4Material* world_mat = nist -> FindOrBuildMaterial("G4_Galactic");
  //G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");		
  G4double world_size = 10*m;

  G4Box* solidWorld =
	new G4Box("World",                       // its name
		3*world_size,                // half x
		3*world_size,                // half y
		3*world_size);               // half zNDDetectorConstruction::LANDDetectorConstruction()


  G4LogicalVolume* logicWorld =
	new G4LogicalVolume(solidWorld,          //its solid
		world_mat,           //its material
		"World");            //its name

  G4VPhysicalVolume* physWorld =
	new G4PVPlacement(0,                     //no rotation
		G4ThreeVector(),       //at (0,0,0)
		logicWorld,            //its logical volume
		"World",               //its name
		0,                     //its mother  volume
		false,                 //no boolean operation
		0,                     //copy number
		true);                 //overlaps checking


  // -----------------------------------------------------
  // Detector
  const int blkcor = 4;
  const int tel = 4;

  G4Material* scintillator_mat = nist -> FindOrBuildMaterial("CsI");
  G4Material* Silicon_mat = nist -> FindOrBuildMaterial("G4_Si");
  G4Material* randomizer_mat = nist -> FindOrBuildMaterial("G4_Al");

  G4Box*randomizer = new G4Box("target",(20*mm)/2,(20*mm)/2,(0.1*mm)/2);
  G4LogicalVolume* logicalDetector0 
	= new G4LogicalVolume(randomizer,randomizer_mat,"target");

  G4VPhysicalVolume* PhysicalDetector0 
	= new G4PVPlacement(0,
		G4ThreeVector(0*mm,0*mm,0*mm),
		logicalDetector0,
		"target",
		logicWorld,
		false,
		1,
		true);

  /*new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_1)*copyN*degree, Angle1, 0*degree), G4ThreeVector(0*mm, 0*mm, 400*mm).rotateX(Angle1).rotateZ((360/CsIPhiN_1)*copyN*degree),"Si_PV", Si_LV_1, labPV, false, 111+copyN);
	1000*(TMath::Sin(0.00142*(Blk+1)))*(TMath::Cos((pi/4.0)+Tel*(pi/2.0))*mm) , 1000*(TMath::Sin(0.00142*(Blk+1)))*(TMath::Sin((pi/4.0)+Tel*(pi/2.0)))*mm , 1000*(TMath::Cos(0.00142*(Blk+1)))*mm)*/

//////////////COR////////////////////////////////////////////////////////////////////detecor kind -> BlockNum -> Quad -> Tele
  G4int silicon_11Num[blkcor][tel] = {{11,22,33,44},
	{111,122,133,144},
	{211,222,233,244},
	{311,322,333,344}};

  G4int silicon_12Num[blkcor][tel] = {{2011,2022,2033,2044},
	{2111,2122,2133,2144},
	{2211,2222,2233,2244},
	{2311,2322,2333,2344}};

  G4int csicrystal1Num[blkcor][tel] = {{3011,3022,3033,3044},
	{3111,3122,3133,3144},
	{3211,3222,3233,3244},
	{3311,3322,3333,3344}};

  G4int silicon_21Num[blkcor][tel] = {{12,23,34,41},
	{112,123,134,141},
	{212,223,234,241},
	{312,323,334,341}};

  G4int silicon_22Num[blkcor][tel] = {{2012,2023,2034,2041},
	{2112,2123,2134,2141},
	{2212,2223,2234,2241},
	{2312,2323,2334,2341}};

  G4int csicrystal2Num[blkcor][tel] = {{3012,3023,3034,3041},
	{3112,3123,3134,3141},
	{3212,3223,3234,3241},
	{3312,3323,3334,3341}};

  G4int silicon_31Num[blkcor][tel] = {{14,21,32,43},
	{114,121,132,143},
	{214,221,232,243},
	{314,321,332,343}};

  G4int silicon_32Num[blkcor][tel] = {{2014,2021,2032,2043},
	{2114,2121,2132,2143},
	{2214,2221,2232,2243},
	{2314,2321,2332,2343}};

  G4int csicrystal3Num[blkcor][tel] = {{3014,3021,3032,3043},
	{3114,3121,3132,3143},
	{3214,3221,3232,3243},
	{3314,3321,3332,3343}};

  G4int silicon_41Num[blkcor][tel] = {{13,24,31,42},
	{113,124,131,142},
	{213,224,231,242},
	{313,324,331,342}};

  G4int silicon_42Num[blkcor][tel] = {{2013,2024,2031,2042},
	{2113,2124,2131,2142},
	{2213,2224,2231,2242},
	{2313,2324,2331,2342}};

  G4int csicrystal4Num[blkcor][tel] = {{3013,3024,3031,3042},
	{3113,3124,3131,3142},
	{3213,3224,3231,3242},
	{3313,3324,3331,3342}};
//////////////////////////////////////////////////////////////////////////////
  string silicon_11[blkcor][tel] ={ {"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"},
	{"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"},
	{"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"},
	{"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"}};

  string silicon_12[blkcor][tel] ={ {"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"},
	{"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"},
	{"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"},
	{"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"}};

  string csicrystal1[blkcor][tel] ={ {"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"},
	{"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"},
	{"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"},
	{"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"}};

  string silicon_21[blkcor][tel] ={ {"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"},
	{"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"},
	{"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"},
	{"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"}};

  string silicon_22[blkcor][tel] ={ {"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"},
	{"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"},
	{"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"},
	{"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"}};

  string csicrystal2[blkcor][tel] = {{"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"},
	{"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"},
	{"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"},
	{"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"}};

  string silicon_31[blkcor][tel] ={ {"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"},
	{"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"},
	{"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"},
	{"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"}};

  string silicon_32[blkcor][tel] ={ {"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"},
	{"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"},
	{"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"},
	{"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"}};

  string csicrystal3[blkcor][tel] ={ {"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"},
	{"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"},
	{"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"},
	{"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"}};

  string silicon_41[blkcor][tel] ={ {"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"},
	{"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"},
	{"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"},
	{"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"}};

  string silicon_42[blkcor][tel] ={ {"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"},
	{"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"},
	{"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"},
	{"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"}};

  string csicrystal4[blkcor][tel] ={ {"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"},
	{"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"},
	{"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"},
	{"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"}};

  string TubeSurface[blkcor][tel] ={ {"CsI11","CsI12","CsI13","CsI14"},
	{"CsI11","CsI12","CsI13","CsI14"},
	{"CsI11","CsI12","CsI13","CsI14"},
	{"CsI11","CsI12","CsI13","CsI14"}};

  G4Box* Silicon_11[blkcor][tel];
  G4Box* Silicon_12[blkcor][tel];
  G4Trd* CsIcrystal1[blkcor][tel];

  G4Box* Silicon_21[blkcor][tel];
  G4Box* Silicon_22[blkcor][tel];
  G4Trd* CsIcrystal2[blkcor][tel];

  G4Box* Silicon_31[blkcor][tel];
  G4Box* Silicon_32[blkcor][tel];
  G4Trd* CsIcrystal3[blkcor][tel];

  G4Box* Silicon_41[blkcor][tel];
  G4Box* Silicon_42[blkcor][tel];
  G4Trd* CsIcrystal4[blkcor][tel];
  //G4Trd* CsIcrystal[1][qua][blkcor][tel];

  G4LogicalVolume* logicalDetector11[blkcor][tel];
  G4LogicalVolume* logicalDetector12[blkcor][tel];
  G4LogicalVolume* logicalDetector13[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector11[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector12[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector13[blkcor][tel];

  G4LogicalVolume* logicalDetector21[blkcor][tel];
  G4LogicalVolume* logicalDetector22[blkcor][tel];
  G4LogicalVolume* logicalDetector23[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector21[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector22[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector23[blkcor][tel];

  G4LogicalVolume* logicalDetector31[blkcor][tel];
  G4LogicalVolume* logicalDetector32[blkcor][tel];
  G4LogicalVolume* logicalDetector33[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector31[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector32[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector33[blkcor][tel];

  G4LogicalVolume* logicalDetector41[blkcor][tel];
  G4LogicalVolume* logicalDetector42[blkcor][tel];
  G4LogicalVolume* logicalDetector43[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector41[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector42[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector43[blkcor][tel];


  G4VisAttributes* Silicon_11VisAtt[4][4];
  G4VisAttributes* Silicon_12VisAtt[4][4];
  G4VisAttributes* CsIcrystal1VisAtt[4][4];
  G4VisAttributes* Silicon_21VisAtt[4][4];
  G4VisAttributes* Silicon_22VisAtt[4][4];
  G4VisAttributes* CsIcrystal2VisAtt[4][4];
  G4VisAttributes* Silicon_31VisAtt[4][4];
  G4VisAttributes* Silicon_32VisAtt[4][4];
  G4VisAttributes* CsIcrystal3VisAtt[4][4];
  G4VisAttributes* Silicon_41VisAtt[4][4];
  G4VisAttributes* Silicon_42VisAtt[4][4];
  G4VisAttributes* CsIcrystal4VisAtt[4][4];
  /////////////////////////////////////////////////////////////////////////////////////
  
//////////////OUTER////////////////////////////////////////////////////////////////////5,7,9,11
  G4int silicon_511Num[blkcor][tel] = {{511,522,533,544},
	{711,722,733,744},
	{911,922,933,944},
	{1111,1122,1133,1144}};

  G4int silicon_512Num[blkcor][tel] = {{2511,2522,2533,2544},
	{2711,2722,2733,2744},
	{2911,2922,2933,2944},
	{21111,21122,21133,21144}};

  G4int csicrystal51Num[blkcor][tel] = {{3511,3522,3533,3544},
	{3711,3722,3733,3744},
	{3911,3922,3933,3944},
	{31111,31122,31133,31144}};

  G4int silicon_521Num[blkcor][tel] = {{512,523,534,541},
	{712,723,734,741},
	{912,923,934,941},
	{1112,1123,1134,1141}};

  G4int silicon_522Num[blkcor][tel] = {{2512,2523,2534,2541},
	{2712,2723,2734,2741},
	{2912,2923,2934,2941},
	{21112,21123,21134,21141}};

  G4int csicrystal52Num[blkcor][tel] = {{3512,3523,3534,3541},
	{3712,3723,3734,3741},
	{3912,3923,3934,3941},
	{31112,31123,31134,31141}};

  G4int silicon_531Num[blkcor][tel] = {{514,521,532,543},
	{714,721,732,743},
	{914,921,932,943},
	{1114,1121,1132,1143}};

  G4int silicon_532Num[blkcor][tel] = {{2514,2521,2532,2543},
	{2714,2721,2732,2743},
	{2914,2921,2932,2943},
	{21114,21121,21132,21143}};

  G4int csicrystal53Num[blkcor][tel] = {{3514,3521,3532,3543},
	{3714,3721,3732,3743},
	{3914,3921,3932,3943},
	{31114,31121,31132,31143}};

  G4int silicon_541Num[blkcor][tel] = {{513,524,531,542},
	{713,724,731,742},
	{913,924,931,942},
	{1113,1124,1131,1142}};

  G4int silicon_542Num[blkcor][tel] = {{2513,2524,2531,2542},
	{2713,2724,2731,2742},
	{2913,2924,2931,2942},
	{21113,21124,21131,21142}};

  G4int csicrystal54Num[blkcor][tel] = {{3513,3524,3531,3542},
	{3713,3724,3731,3742},
	{3913,3924,3931,3942},
	{31113,31124,31131,31142}};
//////////////////////////////////////////////////////////////////////////////
  string silicon_511[blkcor][tel] ={ {"S1B5Q1T1","S1B5Q2T2","S1B5Q3T3","S1B5Q4T4"},
	{"S1B7Q1T1","S1B7Q2T2","S1B7Q3T3","S1B7Q4T4"},
	{"S1B9Q1T1","S1B9Q2T2","S1B9Q3T3","S1B9Q4T4"},
	{"S1B11Q1T1","S1B11Q2T2","S1B11Q3T3","S1B11Q4T4"}};

  string silicon_512[blkcor][tel] ={ {"S2B5Q1T1","S2B5Q2T2","S2B5Q3T3","S2B5Q4T4"},
	{"S2B7Q1T1","S2B7Q2T2","S2B7Q3T3","S2B7Q4T4"},
	{"S2B9Q1T1","S2B9Q2T2","S2B9Q3T3","S2B9Q4T4"},
	{"S2B11Q1T1","S2B11Q2T2","S2B11Q3T3","S2B11Q4T4"}};

  string csicrystal51[blkcor][tel] ={ {"CsIB5Q1T1","CsIB5Q2T2","CsIB5Q3T3","CsIB5Q4T4"},
	{"CsIB7Q1T1","CsIB7Q2T2","CsIB7Q3T3","CsIB7Q4T4"},
	{"CsIB9Q1T1","CsIB9Q2T2","CsIB9Q3T3","CsIB9Q4T4"},
	{"CsIB11Q1T1","CsIB11Q2T2","CsIB11Q3T3","CsIB11Q4T4"}};

  string silicon_521[blkcor][tel] ={ {"S1B5Q1T2","S1B5Q2T3","S1B5Q3T4","S1B5Q4T1"},
	{"S1B7Q1T2","S1B7Q2T3","S1B7Q3T4","S1B7Q4T1"},
	{"S1B9Q1T2","S1B9Q2T3","S1B9Q3T4","S1B9Q4T1"},
	{"S1B11Q1T2","S1B11Q2T3","S1B11Q3T4","S1B11Q4T1"}};

  string silicon_522[blkcor][tel] ={ {"S2B5Q1T2","S2B5Q2T3","S2B5Q3T4","S2B5Q4T1"},
	{"S2B7Q1T2","S2B7Q2T3","S2B7Q3T4","S2B7Q4T1"},
	{"S2B9Q1T2","S2B9Q2T3","S2B9Q3T4","S2B9Q4T1"},
	{"S2B11Q1T2","S2B11Q2T3","S2B11Q3T4","S2B11Q4T1"}};

  string csicrystal52[blkcor][tel] = {{"CsIB5Q1T2","CsIB5Q2T3","CsIB5Q3T4","CsIB5Q4T1"},
	{"CsIB7Q1T2","CsIB7Q2T3","CsIB7Q3T4","CsIB7Q4T1"},
	{"CsIB9Q1T2","CsIB9Q2T3","CsIB9Q3T4","CsIB9Q4T1"},
	{"CsIB11Q1T2","CsIB11Q2T3","CsIB11Q3T4","CsIB11Q4T1"}};

  string silicon_531[blkcor][tel] ={ {"S1B5Q1T4","S1B5Q2T1","S1B5Q3T2","S1B5Q4T3"},
	{"S1B7Q1T4","S1B7Q2T1","S1B7Q3T2","S1B7Q4T3"},
	{"S1B9Q1T4","S1B9Q2T1","S1B9Q3T2","S1B9Q4T3"},
	{"S1B11Q1T4","S1B11Q2T1","S1B11Q3T2","S1B11Q4T3"}};

  string silicon_532[blkcor][tel] ={ {"S2B5Q1T4","S2B5Q2T1","S2B5Q3T2","S2B5Q4T3"},
	{"S2B7Q1T4","S2B7Q2T1","S2B7Q3T2","S2B7Q4T3"},
	{"S2B9Q1T4","S2B9Q2T1","S2B9Q3T2","S2B9Q4T3"},
	{"S2B11Q1T4","S2B11Q2T1","S2B11Q3T2","S2B11Q4T3"}};

  string csicrystal53[blkcor][tel] ={ {"CsIB5Q1T4","CsIB5Q2T1","CsIB5Q3T2","CsIB5Q4T3"},
	{"CsIB7Q1T4","CsIB7Q2T1","CsIB7Q3T2","CsIB7Q4T3"},
	{"CsIB9Q1T4","CsIB9Q2T1","CsIB9Q3T2","CsIB9Q4T3"},
	{"CsIB11Q1T4","CsIB11Q2T1","CsIB11Q3T2","CsIB11Q4T3"}};

  string silicon_541[blkcor][tel] ={ {"S1B5Q1T3","S1B5Q2T4","S1B5Q3T1","S1B5Q4T2"},
	{"S1B7Q1T3","S1B7Q2T4","S1B7Q3T1","S1B7Q4T2"},
	{"S1B9Q1T3","S1B9Q2T4","S1B9Q3T1","S1B9Q4T2"},
	{"S1B11Q1T3","S1B11Q2T4","S1B11Q3T1","S1B11Q4T2"}};

  string silicon_542[blkcor][tel] ={ {"S2B5Q1T3","S2B5Q2T4","S2B5Q3T1","S2B5Q4T2"},
	{"S2B7Q1T3","S2B7Q2T4","S2B7Q3T1","S2B7Q4T2"},
	{"S2B9Q1T3","S2B9Q2T4","S2B9Q3T1","S2B9Q4T2"},
	{"S2B11Q1T3","S2B11Q2T4","S2B11Q3T1","S2B11Q4T2"}};

  string csicrystal54[blkcor][tel] ={ {"CsIB5Q1T3","CsIB5Q2T4","CsIB5Q3T1","CsIB5Q4T2"},
	{"CsIB7Q1T3","CsIB7Q2T4","CsIB7Q3T1","CsIB7Q4T2"},
	{"CsIB9Q1T3","CsIB9Q2T4","CsIB9Q3T1","CsIB9Q4T2"},
	{"CsIB11Q1T3","CsIB11Q2T4","CsIB11Q3T1","CsIB11Q4T2"}};

  string TubeSurface5[blkcor][tel] ={ {"CsI511","CsI512","CsI513","CsI514"},
	{"CsI711","CsI712","CsI713","CsI714"},
	{"CsI911","CsI912","CsI913","CsI914"},
	{"CsI1111","CsI1112","CsI1113","CsI1114"}};

  G4Box* Silicon_511[blkcor][tel];
  G4Box* Silicon_512[blkcor][tel];
  G4Trd* CsIcrystal51[blkcor][tel];

  G4Box* Silicon_521[blkcor][tel];
  G4Box* Silicon_522[blkcor][tel];
  G4Trd* CsIcrystal52[blkcor][tel];

  G4Box* Silicon_531[blkcor][tel];
  G4Box* Silicon_532[blkcor][tel];
  G4Trd* CsIcrystal53[blkcor][tel];

  G4Box* Silicon_541[blkcor][tel];
  G4Box* Silicon_542[blkcor][tel];
  G4Trd* CsIcrystal54[blkcor][tel];
  //G4Trd* CsIcrystal[1][qua][blkcor][tel];

  G4LogicalVolume* logicalDetector511[blkcor][tel];
  G4LogicalVolume* logicalDetector512[blkcor][tel];
  G4LogicalVolume* logicalDetector513[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector511[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector512[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector513[blkcor][tel];

  G4LogicalVolume* logicalDetector521[blkcor][tel];
  G4LogicalVolume* logicalDetector522[blkcor][tel];
  G4LogicalVolume* logicalDetector523[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector521[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector522[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector523[blkcor][tel];

  G4LogicalVolume* logicalDetector531[blkcor][tel];
  G4LogicalVolume* logicalDetector532[blkcor][tel];
  G4LogicalVolume* logicalDetector533[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector531[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector532[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector533[blkcor][tel];

  G4LogicalVolume* logicalDetector541[blkcor][tel];
  G4LogicalVolume* logicalDetector542[blkcor][tel];
  G4LogicalVolume* logicalDetector543[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector541[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector542[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector543[blkcor][tel];


  G4VisAttributes* Silicon_511VisAtt[4][4];
  G4VisAttributes* Silicon_512VisAtt[4][4];
  G4VisAttributes* CsIcrystal51VisAtt[4][4];
  G4VisAttributes* Silicon_521VisAtt[4][4];
  G4VisAttributes* Silicon_522VisAtt[4][4];
  G4VisAttributes* CsIcrystal52VisAtt[4][4];
  G4VisAttributes* Silicon_531VisAtt[4][4];
  G4VisAttributes* Silicon_532VisAtt[4][4];
  G4VisAttributes* CsIcrystal53VisAtt[4][4];
  G4VisAttributes* Silicon_541VisAtt[4][4];
  G4VisAttributes* Silicon_542VisAtt[4][4];
  G4VisAttributes* CsIcrystal54VisAtt[4][4];
  /////////////////////////////////////////////////////////////////////////////////////

//////////////----MID----////////////////////////////////////////////////////////////////////4,6,8,10
  G4int silicon_411Num[blkcor][tel] = {{411,422,433,444},
	{611,622,633,644},
	{811,822,833,844},
	{1011,1022,1033,1044}};

  G4int silicon_412Num[blkcor][tel] = {{2411,2422,2433,2444},
	{2611,2622,2633,2644},
	{2811,2822,2833,2844},
	{21011,21022,21033,21044}};

  G4int csicrystal41Num[blkcor][tel] = {{3411,3422,3433,3444},
	{3611,3622,3633,3644},
	{3811,3822,3833,3844},
	{31011,31022,31033,31044}};

  G4int silicon_421Num[blkcor][tel] = {{412,423,434,441},
	{612,623,634,641},
	{812,823,834,841},
	{1012,1023,1034,1041}};

  G4int silicon_422Num[blkcor][tel] = {{2412,2423,2434,2441},
	{2612,2623,2634,2641},
	{2812,2823,2834,2841},
	{21012,21023,21034,21041}};

  G4int csicrystal42Num[blkcor][tel] = {{3412,3423,3434,3441},
	{3612,3623,3634,3641},
	{3812,3823,3834,3841},
	{31012,31023,31034,31041}};

  G4int silicon_431Num[blkcor][tel] = {{414,421,432,443},
	{614,621,632,643},
	{814,821,832,843},
	{1014,1021,1032,1043}};

  G4int silicon_432Num[blkcor][tel] = {{2414,2421,2432,2443},
	{2614,2621,2632,2643},
	{2814,2821,2832,2843},
	{21014,21021,21032,21043}};

  G4int csicrystal43Num[blkcor][tel] = {{3414,3421,3432,3443},
	{3614,3621,3632,3643},
	{3814,3821,3832,3843},
	{31014,31021,31032,31043}};

  G4int silicon_441Num[blkcor][tel] = {{413,424,431,442},
	{613,624,631,642},
	{813,824,831,842},
	{1013,1024,1031,1042}};

  G4int silicon_442Num[blkcor][tel] = {{2413,2424,2431,2442},
	{2613,2624,2631,2642},
	{2813,2824,2831,2842},
	{21013,21024,21031,21042}};

  G4int csicrystal44Num[blkcor][tel] = {{3413,3424,3431,3442},
	{3613,3624,3631,3642},
	{3813,3824,3831,3842},
	{31013,31024,31031,31042}};
//////////////////////////////////////////////////////////////////////////////
  string silicon_411[blkcor][tel] ={ {"S1B4Q1T1","S1B4Q2T2","S1B4Q3T3","S1B4Q4T4"},
	{"S1B6Q1T1","S1B6Q2T2","S1B6Q3T3","S1B6Q4T4"},
	{"S1B8Q1T1","S1B8Q2T2","S1B8Q3T3","S1B8Q4T4"},
	{"S1B10Q1T1","S1B10Q2T2","S1B10Q3T3","S1B10Q4T4"}};

  string silicon_412[blkcor][tel] ={ {"S2B4Q1T1","S2B4Q2T2","S2B4Q3T3","S2B4Q4T4"},
	{"S2B6Q1T1","S2B6Q2T2","S2B6Q3T3","S2B6Q4T4"},
	{"S2B8Q1T1","S2B8Q2T2","S2B8Q3T3","S2B8Q4T4"},
	{"S2B10Q1T1","S2B10Q2T2","S2B10Q3T3","S2B10Q4T4"}};

  string csicrystal41[blkcor][tel] ={ {"CsIB4Q1T1","CsIB4Q2T2","CsIB4Q3T3","CsIB4Q4T4"},
	{"CsIB6Q1T1","CsIB6Q2T2","CsIB6Q3T3","CsIB6Q4T4"},
	{"CsIB8Q1T1","CsIB8Q2T2","CsIB8Q3T3","CsIB8Q4T4"},
	{"CsIB10Q1T1","CsIB10Q2T2","CsIB10Q3T3","CsIB10Q4T4"}};

  string silicon_421[blkcor][tel] ={ {"S1B4Q1T2","S1B4Q2T3","S1B4Q3T4","S1B4Q4T1"},
	{"S1B6Q1T2","S1B6Q2T3","S1B6Q3T4","S1B6Q4T1"},
	{"S1B8Q1T2","S1B8Q2T3","S1B8Q3T4","S1B8Q4T1"},
	{"S1B10Q1T2","S1B10Q2T3","S1B10Q3T4","S1B10Q4T1"}};

  string silicon_422[blkcor][tel] ={ {"S2B4Q1T2","S2B4Q2T3","S2B4Q3T4","S2B4Q4T1"},
	{"S2B6Q1T2","S2B6Q2T3","S2B6Q3T4","S2B6Q4T1"},
	{"S2B8Q1T2","S2B8Q2T3","S2B8Q3T4","S2B8Q4T1"},
	{"S2B10Q1T2","S2B10Q2T3","S2B10Q3T4","S2B10Q4T1"}};

  string csicrystal42[blkcor][tel] = {{"CsIB4Q1T2","CsIB4Q2T3","CsIB4Q3T4","CsIB4Q4T1"},
	{"CsIB6Q1T2","CsIB6Q2T3","CsIB6Q3T4","CsIB6Q4T1"},
	{"CsIB8Q1T2","CsIB8Q2T3","CsIB8Q3T4","CsIB8Q4T1"},
	{"CsIB10Q1T2","CsIB10Q2T3","CsIB10Q3T4","CsIB10Q4T1"}};

  string silicon_431[blkcor][tel] ={ {"S1B4Q1T4","S1B4Q2T1","S1B4Q3T2","S1B4Q4T3"},
	{"S1B6Q1T4","S1B6Q2T1","S1B6Q3T2","S1B6Q4T3"},
	{"S1B8Q1T4","S1B8Q2T1","S1B8Q3T2","S1B8Q4T3"},
	{"S1B10Q1T4","S1B10Q2T1","S1B10Q3T2","S1B10Q4T3"}};

  string silicon_432[blkcor][tel] ={ {"S2B4Q1T4","S2B4Q2T1","S2B4Q3T2","S2B4Q4T3"},
	{"S2B6Q1T4","S2B6Q2T1","S2B6Q3T2","S2B6Q4T3"},
	{"S2B8Q1T4","S2B8Q2T1","S2B8Q3T2","S2B8Q4T3"},
	{"S2B10Q1T4","S2B10Q2T1","S2B10Q3T2","S2B10Q4T3"}};

  string csicrystal43[blkcor][tel] ={ {"CsIB4Q1T4","CsIB4Q2T1","CsIB4Q3T2","CsIB4Q4T3"},
	{"CsIB6Q1T4","CsIB6Q2T1","CsIB6Q3T2","CsIB6Q4T3"},
	{"CsIB8Q1T4","CsIB8Q2T1","CsIB8Q3T2","CsIB8Q4T3"},
	{"CsIB10Q1T4","CsIB10Q2T1","CsIB10Q3T2","CsIB10Q4T3"}};

  string silicon_441[blkcor][tel] ={ {"S1B4Q1T3","S1B4Q2T4","S1B4Q3T1","S1B4Q4T2"},
	{"S1B6Q1T3","S1B6Q2T4","S1B6Q3T1","S1B6Q4T2"},
	{"S1B8Q1T3","S1B8Q2T4","S1B8Q3T1","S1B8Q4T2"},
	{"S1B10Q1T3","S1B10Q2T4","S1B10Q3T1","S1B10Q4T2"}};

  string silicon_442[blkcor][tel] ={ {"S2B4Q1T3","S2B4Q2T4","S2B4Q3T1","S2B4Q4T2"},
	{"S2B6Q1T3","S2B6Q2T4","S2B6Q3T1","S2B6Q4T2"},
	{"S2B8Q1T3","S2B8Q2T4","S2B8Q3T1","S2B8Q4T2"},
	{"S2B10Q1T3","S2B10Q2T4","S2B10Q3T1","S2B10Q4T2"}};

  string csicrystal44[blkcor][tel] ={ {"CsIB4Q1T3","CsIB4Q2T4","CsIB4Q3T1","CsIB4Q4T2"},
	{"CsIB6Q1T3","CsIB6Q2T4","CsIB6Q3T1","CsIB6Q4T2"},
	{"CsIB8Q1T3","CsIB8Q2T4","CsIB8Q3T1","CsIB8Q4T2"},
	{"CsIB10Q1T3","CsIB10Q2T4","CsIB10Q3T1","CsIB10Q4T2"}};

  string TubeSurface4[blkcor][tel] ={ {"CsI411","CsI412","CsI413","CsI414"},
	{"CsI611","CsI612","CsI613","CsI614"},
	{"CsI811","CsI812","CsI813","CsI814"},
	{"CsI1011","CsI1012","CsI1013","CsI1014"}};

  G4Box* Silicon_411[blkcor][tel];
  G4Box* Silicon_412[blkcor][tel];
  G4Trd* CsIcrystal41[blkcor][tel];

  G4Box* Silicon_421[blkcor][tel];
  G4Box* Silicon_422[blkcor][tel];
  G4Trd* CsIcrystal42[blkcor][tel];

  G4Box* Silicon_431[blkcor][tel];
  G4Box* Silicon_432[blkcor][tel];
  G4Trd* CsIcrystal43[blkcor][tel];

  G4Box* Silicon_441[blkcor][tel];
  G4Box* Silicon_442[blkcor][tel];
  G4Trd* CsIcrystal44[blkcor][tel];
  //G4Trd* CsIcrystal[1][qua][blkcor][tel];

  G4LogicalVolume* logicalDetector411[blkcor][tel];
  G4LogicalVolume* logicalDetector412[blkcor][tel];
  G4LogicalVolume* logicalDetector413[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector411[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector412[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector413[blkcor][tel];

  G4LogicalVolume* logicalDetector421[blkcor][tel];
  G4LogicalVolume* logicalDetector422[blkcor][tel];
  G4LogicalVolume* logicalDetector423[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector421[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector422[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector423[blkcor][tel];

  G4LogicalVolume* logicalDetector431[blkcor][tel];
  G4LogicalVolume* logicalDetector432[blkcor][tel];
  G4LogicalVolume* logicalDetector433[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector431[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector432[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector433[blkcor][tel];

  G4LogicalVolume* logicalDetector441[blkcor][tel];
  G4LogicalVolume* logicalDetector442[blkcor][tel];
  G4LogicalVolume* logicalDetector443[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector441[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector442[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector443[blkcor][tel];


  G4VisAttributes* Silicon_411VisAtt[4][4];
  G4VisAttributes* Silicon_412VisAtt[4][4];
  G4VisAttributes* CsIcrystal41VisAtt[4][4];
  G4VisAttributes* Silicon_421VisAtt[4][4];
  G4VisAttributes* Silicon_422VisAtt[4][4];
  G4VisAttributes* CsIcrystal42VisAtt[4][4];
  G4VisAttributes* Silicon_431VisAtt[4][4];
  G4VisAttributes* Silicon_432VisAtt[4][4];
  G4VisAttributes* CsIcrystal43VisAtt[4][4];
  G4VisAttributes* Silicon_441VisAtt[4][4];
  G4VisAttributes* Silicon_442VisAtt[4][4];
  G4VisAttributes* CsIcrystal44VisAtt[4][4];
  /////////////////////////////////////////////////////////////////////////////////////
  G4OpticalSurface* opTubeSurface[blkcor][tel];

  logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());
  // for(G4int Blk = 0 ; Blk < 1 ; Blk++)
  //{
  //for(G4int Qua = 0 ; Qua < 1 ; Qua++)
  //{
  G4double AngleN = -(0.810)*CLHEP::degree;
  G4double AngleNN = -(0.80)*CLHEP::degree;

  G4RotationMatrix* Rot_1Si[blkcor][tel];
  G4RotationMatrix* Rot_2Si[blkcor][tel];
  G4RotationMatrix* Rot_3Si[blkcor][tel];
  G4RotationMatrix* Rot_4Si[blkcor][tel];

  //----------------------------------------------------------------------------------------------------------
  for(G4int Blkcor = 0 ; Blkcor < 4 ; Blkcor++)
  {
	for(G4int Tel = 0 ; Tel < 4 ; Tel++)
	{
	  Rot_1Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_2Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_3Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_4Si[Blkcor][Tel] = new G4RotationMatrix;

	  if(Tel == 0){
		Rot_1Si[Blkcor][Tel] -> rotateY(AngleNN*2);
		Rot_1Si[Blkcor][Tel] -> rotateX(AngleNN*-2);}
	  if(Tel == 1){
		Rot_1Si[Blkcor][Tel] -> rotateY(AngleNN*-2);
		Rot_1Si[Blkcor][Tel] -> rotateX(AngleNN*-2);}
	  if(Tel == 2){
		Rot_1Si[Blkcor][Tel] -> rotateY(AngleNN*-2);
		Rot_1Si[Blkcor][Tel] -> rotateX(AngleNN*2);}
	  if(Tel == 3){
		Rot_1Si[Blkcor][Tel] -> rotateY(AngleNN*2);
		Rot_1Si[Blkcor][Tel] -> rotateX(AngleNN*2);}

	  if(Tel == 0){
		Rot_2Si[Blkcor][Tel] -> rotateY(AngleN*1);
		Rot_2Si[Blkcor][Tel] -> rotateX(AngleNN*-2);}
	  if(Tel == 1){
		Rot_2Si[Blkcor][Tel] -> rotateY(AngleNN*-2);
		Rot_2Si[Blkcor][Tel] -> rotateX(AngleN*-1);}
	  if(Tel == 2){
		Rot_2Si[Blkcor][Tel] -> rotateY(AngleN*-1);
		Rot_2Si[Blkcor][Tel] -> rotateX(AngleNN*2);}
	  if(Tel == 3){
		Rot_2Si[Blkcor][Tel] -> rotateY(AngleNN*2);
		Rot_2Si[Blkcor][Tel] -> rotateX(AngleN*1);}

	  if(Tel == 0){
		Rot_3Si[Blkcor][Tel] -> rotateY(AngleNN*2);
		Rot_3Si[Blkcor][Tel] -> rotateX(AngleN*-1);}
	  if(Tel == 1){
		Rot_3Si[Blkcor][Tel] -> rotateY(AngleN*-1);
		Rot_3Si[Blkcor][Tel] -> rotateX(AngleNN*-2);}
	  if(Tel == 2){
		Rot_3Si[Blkcor][Tel] -> rotateY(AngleNN*-2);
		Rot_3Si[Blkcor][Tel] -> rotateX(AngleN*1);}
	  if(Tel == 3){
		Rot_3Si[Blkcor][Tel] -> rotateY(AngleN*1);
		Rot_3Si[Blkcor][Tel] -> rotateX(AngleNN*2);}

	  if(Tel == 0){
		Rot_4Si[Blkcor][Tel] -> rotateY(AngleN*1);
		Rot_4Si[Blkcor][Tel] -> rotateX(AngleN*-1);}
	  if(Tel == 1){
		Rot_4Si[Blkcor][Tel] -> rotateY(AngleN*-1);
		Rot_4Si[Blkcor][Tel] -> rotateX(AngleN*-1);}
	  if(Tel == 2){
		Rot_4Si[Blkcor][Tel] -> rotateY(AngleN*-1);
		Rot_4Si[Blkcor][Tel] -> rotateX(AngleN*1);}
	  if(Tel == 3){
		Rot_4Si[Blkcor][Tel] -> rotateY(AngleN*1);
		Rot_4Si[Blkcor][Tel] -> rotateX(AngleN*1);}

	  //one block's geometry-COR---------------------------------------------------------------------------------------------------------		
	  {//1
		//Si1
		Silicon_11[Blkcor][Tel] = new G4Box(silicon_11[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector11[Blkcor][Tel] = new G4LogicalVolume(Silicon_11[Blkcor][Tel],Silicon_mat,silicon_11[Blkcor][Tel]);
		PhysicalDetector11[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Cos(Blktheta1))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector11[Blkcor][Tel],//theta == 0.0144
			silicon_11[Blkcor][Tel],
			logicWorld,
			false,
			silicon_11Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_12[Blkcor][Tel] = new G4Box(silicon_12[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector12[Blkcor][Tel] = new G4LogicalVolume(Silicon_12[Blkcor][Tel] ,Silicon_mat ,silicon_12[Blkcor][Tel]);

		PhysicalDetector12[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Cos(Blktheta1))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector12[Blkcor][Tel],
			silicon_12[Blkcor][Tel],
			logicWorld,
			false,
			silicon_12Num[Blkcor][Tel],
			true);

		//CsI crystal		  
		CsIcrystal1[Blkcor][Tel] = new G4Trd(csicrystal1[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector13[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal1[Blkcor][Tel] ,scintillator_mat ,csicrystal1[Blkcor][Tel]);

		PhysicalDetector13[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Cos(Blktheta1))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector13[Blkcor][Tel],
			csicrystal1[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal1Num[Blkcor][Tel],
			true);
	  }
	  {//2
		Silicon_21[Blkcor][Tel] = new G4Box(silicon_21[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector21[Blkcor][Tel] = new G4LogicalVolume(Silicon_21[Blkcor][Tel],Silicon_mat,silicon_21[Blkcor][Tel]);
		PhysicalDetector21[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1000.*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(pi*Blkcor/2))*mm, 1000*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(pi*Blkcor/2))*mm, 1000*(TMath::Cos(Blktheta2))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector21[Blkcor][Tel],//theta == 0.0316
			silicon_21[Blkcor][Tel],
			logicWorld,
			false,
			silicon_21Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_22[Blkcor][Tel] = new G4Box(silicon_22[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector22[Blkcor][Tel] = new G4LogicalVolume(Silicon_22[Blkcor][Tel] ,Silicon_mat ,silicon_22[Blkcor][Tel]);

		PhysicalDetector22[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(pi*Blkcor/2))*mm, 1000.625*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(pi*Blkcor/2))*mm, 1000.625*(TMath::Cos(Blktheta2))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector22[Blkcor][Tel],
			silicon_22[Blkcor][Tel],
			logicWorld,
			false,
			silicon_22Num[Blkcor][Tel],
			true);

		//CsI crystal		  
		CsIcrystal2[Blkcor][Tel] = new G4Trd(csicrystal2[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector23[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal2[Blkcor][Tel] ,scintillator_mat ,csicrystal2[Blkcor][Tel]);

		PhysicalDetector23[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(pi*Blkcor/2))*mm, 1051.058*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(pi*Blkcor/2))*mm, 1051.058*(TMath::Cos(Blktheta2))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector23[Blkcor][Tel],
			csicrystal2[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal2Num[Blkcor][Tel],
			true);
	  }
	  {//3
		Silicon_31[Blkcor][Tel] = new G4Box(silicon_31[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector31[Blkcor][Tel] = new G4LogicalVolume(Silicon_31[Blkcor][Tel],Silicon_mat,silicon_31[Blkcor][Tel]);
		PhysicalDetector31[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(pi*Blkcor/2))*mm, 1000*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(pi*Blkcor/2))*mm, 1000*(TMath::Cos(Blktheta3))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector31[Blkcor][Tel],//theta == 0.0147
			silicon_31[Blkcor][Tel],
			logicWorld,
			false,
			silicon_31Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_32[Blkcor][Tel] = new G4Box(silicon_32[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector32[Blkcor][Tel] = new G4LogicalVolume(Silicon_32[Blkcor][Tel] ,Silicon_mat ,silicon_32[Blkcor][Tel]);

		PhysicalDetector32[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(pi*Blkcor/2))*mm, 1000.625*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(pi*Blkcor/2))*mm, 1000.625*(TMath::Cos(Blktheta3))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector32[Blkcor][Tel],
			silicon_32[Blkcor][Tel],
			logicWorld,
			false,
			silicon_32Num[Blkcor][Tel],
			true);

		//CsI crystal		  
		CsIcrystal3[Blkcor][Tel] = new G4Trd(csicrystal3[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector33[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal3[Blkcor][Tel] ,scintillator_mat ,csicrystal3[Blkcor][Tel]);

		PhysicalDetector33[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(pi*Blkcor/2))*mm, 1051.058*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(pi*Blkcor/2))*mm, 1051.058*(TMath::Cos(Blktheta3))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector33[Blkcor][Tel],
			csicrystal3[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal3Num[Blkcor][Tel],
			true);
	  }
	  {//4
		Silicon_41[Blkcor][Tel] = new G4Box(silicon_41[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector41[Blkcor][Tel] = new G4LogicalVolume(Silicon_41[Blkcor][Tel],Silicon_mat,silicon_41[Blkcor][Tel]);
		PhysicalDetector41[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Cos(Blktheta4))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector41[Blkcor][Tel],//theta == 0.0147
			silicon_41[Blkcor][Tel],
			logicWorld,
			false,
			silicon_41Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_42[Blkcor][Tel] = new G4Box(silicon_42[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector42[Blkcor][Tel] = new G4LogicalVolume(Silicon_42[Blkcor][Tel] ,Silicon_mat ,silicon_42[Blkcor][Tel]);

		PhysicalDetector42[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Cos(Blktheta4))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector42[Blkcor][Tel],
			silicon_42[Blkcor][Tel],
			logicWorld,
			false,
			silicon_42Num[Blkcor][Tel],
			true);

		//CsI crystal		  
		CsIcrystal4[Blkcor][Tel] = new G4Trd(csicrystal4[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector43[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal4[Blkcor][Tel] ,scintillator_mat ,csicrystal4[Blkcor][Tel]);

		PhysicalDetector43[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Cos(Blktheta4))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector43[Blkcor][Tel],
			csicrystal4[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal4Num[Blkcor][Tel],
			true);
	  }
	  //-------------------------------------------------------------------------------------------------------------------------------		
	  //one block's geometry-MID---------------------------------------------------------------------------------------------------------		
	  {//1
		//Si1
		Silicon_411[Blkcor][Tel] = new G4Box(silicon_411[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector411[Blkcor][Tel] = new G4LogicalVolume(Silicon_411[Blkcor][Tel],Silicon_mat,silicon_411[Blkcor][Tel]);
		PhysicalDetector411[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta41))*TMath::Cos((Blkphi41+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Sin(Blktheta41))*TMath::Sin((Blkphi41+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Cos(Blktheta41))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector411[Blkcor][Tel],//theta4 == 0.0144
			silicon_411[Blkcor][Tel],
			logicWorld,
			false,
			silicon_411Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_412[Blkcor][Tel] = new G4Box(silicon_412[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector412[Blkcor][Tel] = new G4LogicalVolume(Silicon_412[Blkcor][Tel] ,Silicon_mat ,silicon_412[Blkcor][Tel]);

		PhysicalDetector412[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta41))*TMath::Cos((Blkphi41+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Sin(Blktheta41))*TMath::Sin((Blkphi41+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Cos(Blktheta41))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector412[Blkcor][Tel],
			silicon_412[Blkcor][Tel],
			logicWorld,
			false,
			silicon_412Num[Blkcor][Tel],
			true);

		//CsI crystal4		  
		CsIcrystal41[Blkcor][Tel] = new G4Trd(csicrystal41[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector413[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal41[Blkcor][Tel] ,scintillator_mat ,csicrystal41[Blkcor][Tel]);

		PhysicalDetector413[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta41))*TMath::Cos((Blkphi41+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Sin(Blktheta41))*TMath::Sin((Blkphi41+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Cos(Blktheta41))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector413[Blkcor][Tel],
			csicrystal41[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal41Num[Blkcor][Tel],
			true);
	  }
	  {//2
		Silicon_421[Blkcor][Tel] = new G4Box(silicon_421[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector421[Blkcor][Tel] = new G4LogicalVolume(Silicon_421[Blkcor][Tel],Silicon_mat,silicon_421[Blkcor][Tel]);
		PhysicalDetector421[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1000.*(TMath::Sin(Blktheta42))*TMath::Cos(Blkphi42+(pi*Blkcor/2))*mm, 1000*(TMath::Sin(Blktheta42))*TMath::Sin(Blkphi42+(pi*Blkcor/2))*mm, 1000*(TMath::Cos(Blktheta42))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector421[Blkcor][Tel],//theta4 == 0.0316
			silicon_421[Blkcor][Tel],
			logicWorld,
			false,
			silicon_421Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_422[Blkcor][Tel] = new G4Box(silicon_422[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector422[Blkcor][Tel] = new G4LogicalVolume(Silicon_422[Blkcor][Tel] ,Silicon_mat ,silicon_422[Blkcor][Tel]);

		PhysicalDetector422[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta42))*TMath::Cos(Blkphi42+(pi*Blkcor/2))*mm, 1000.625*(TMath::Sin(Blktheta42))*TMath::Sin(Blkphi42+(pi*Blkcor/2))*mm, 1000.625*(TMath::Cos(Blktheta42))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector422[Blkcor][Tel],
			silicon_422[Blkcor][Tel],
			logicWorld,
			false,
			silicon_422Num[Blkcor][Tel],
			true);

		//CsI crystal4		  
		CsIcrystal42[Blkcor][Tel] = new G4Trd(csicrystal42[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector423[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal42[Blkcor][Tel] ,scintillator_mat ,csicrystal42[Blkcor][Tel]);

		PhysicalDetector423[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta42))*TMath::Cos(Blkphi42+(pi*Blkcor/2))*mm, 1051.058*(TMath::Sin(Blktheta42))*TMath::Sin(Blkphi42+(pi*Blkcor/2))*mm, 1051.058*(TMath::Cos(Blktheta42))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector423[Blkcor][Tel],
			csicrystal42[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal42Num[Blkcor][Tel],
			true);
	  }
	  {//3
		Silicon_431[Blkcor][Tel] = new G4Box(silicon_431[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector431[Blkcor][Tel] = new G4LogicalVolume(Silicon_431[Blkcor][Tel],Silicon_mat,silicon_431[Blkcor][Tel]);
		PhysicalDetector431[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta43))*TMath::Cos(Blkphi43+(pi*Blkcor/2))*mm, 1000*(TMath::Sin(Blktheta43))*TMath::Sin(Blkphi43+(pi*Blkcor/2))*mm, 1000*(TMath::Cos(Blktheta43))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector431[Blkcor][Tel],//theta4 == 0.0147
			silicon_431[Blkcor][Tel],
			logicWorld,
			false,
			silicon_431Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_432[Blkcor][Tel] = new G4Box(silicon_432[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector432[Blkcor][Tel] = new G4LogicalVolume(Silicon_432[Blkcor][Tel] ,Silicon_mat ,silicon_432[Blkcor][Tel]);

		PhysicalDetector432[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta43))*TMath::Cos(Blkphi43+(pi*Blkcor/2))*mm, 1000.625*(TMath::Sin(Blktheta43))*TMath::Sin(Blkphi43+(pi*Blkcor/2))*mm, 1000.625*(TMath::Cos(Blktheta43))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector432[Blkcor][Tel],
			silicon_432[Blkcor][Tel],
			logicWorld,
			false,
			silicon_432Num[Blkcor][Tel],
			true);

		//CsI crystal4		  
		CsIcrystal43[Blkcor][Tel] = new G4Trd(csicrystal43[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector433[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal43[Blkcor][Tel] ,scintillator_mat ,csicrystal43[Blkcor][Tel]);

		PhysicalDetector433[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta43))*TMath::Cos(Blkphi43+(pi*Blkcor/2))*mm, 1051.058*(TMath::Sin(Blktheta43))*TMath::Sin(Blkphi43+(pi*Blkcor/2))*mm, 1051.058*(TMath::Cos(Blktheta43))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector433[Blkcor][Tel],
			csicrystal43[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal43Num[Blkcor][Tel],
			true);
	  }
	  {//4
		Silicon_441[Blkcor][Tel] = new G4Box(silicon_441[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector441[Blkcor][Tel] = new G4LogicalVolume(Silicon_441[Blkcor][Tel],Silicon_mat,silicon_441[Blkcor][Tel]);
		PhysicalDetector441[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta44))*TMath::Cos((Blkphi44+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Sin(Blktheta44))*TMath::Sin((Blkphi44+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Cos(Blktheta44))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector441[Blkcor][Tel],//theta4 == 0.0147
			silicon_441[Blkcor][Tel],
			logicWorld,
			false,
			silicon_441Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_442[Blkcor][Tel] = new G4Box(silicon_442[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector442[Blkcor][Tel] = new G4LogicalVolume(Silicon_442[Blkcor][Tel] ,Silicon_mat ,silicon_442[Blkcor][Tel]);

		PhysicalDetector442[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta44))*TMath::Cos((Blkphi44+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Sin(Blktheta44))*TMath::Sin((Blkphi44+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Cos(Blktheta44))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector442[Blkcor][Tel],
			silicon_442[Blkcor][Tel],
			logicWorld,
			false,
			silicon_442Num[Blkcor][Tel],
			true);

		//CsI crystal4		  
		CsIcrystal44[Blkcor][Tel] = new G4Trd(csicrystal44[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector443[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal44[Blkcor][Tel] ,scintillator_mat ,csicrystal44[Blkcor][Tel]);

		PhysicalDetector443[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta44))*TMath::Cos((Blkphi44+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Sin(Blktheta44))*TMath::Sin((Blkphi44+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Cos(Blktheta44))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.11352728)*TMath::Cos(0.66104317+(pi*Blkcor/2)),TMath::Sin(0.11352728)*TMath::Sin(0.66104317+(pi*Blkcor/2)),TMath::Cos(0.11352728))),
			logicalDetector443[Blkcor][Tel],
			csicrystal44[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal44Num[Blkcor][Tel],
			true);
	  }
	  //--------------------------------------------------------------------------------------------------------------------------------
	  //one block's geometry-OUTER---------------------------------------------------------------------------------------------------------		
	  {//1
		//Si1
		Silicon_511[Blkcor][Tel] = new G4Box(silicon_511[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector511[Blkcor][Tel] = new G4LogicalVolume(Silicon_511[Blkcor][Tel],Silicon_mat,silicon_511[Blkcor][Tel]);
		PhysicalDetector511[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta51))*TMath::Cos((Blkphi51+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Sin(Blktheta51))*TMath::Sin((Blkphi51+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Cos(Blktheta51))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector511[Blkcor][Tel],//theta == 0.0144
			silicon_511[Blkcor][Tel],
			logicWorld,
			false,
			silicon_511Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_512[Blkcor][Tel] = new G4Box(silicon_512[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector512[Blkcor][Tel] = new G4LogicalVolume(Silicon_512[Blkcor][Tel] ,Silicon_mat ,silicon_512[Blkcor][Tel]);

		PhysicalDetector512[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta51))*TMath::Cos((Blkphi51+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Sin(Blktheta51))*TMath::Sin((Blkphi51+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Cos(Blktheta51))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector512[Blkcor][Tel],
			silicon_512[Blkcor][Tel],
			logicWorld,
			false,
			silicon_512Num[Blkcor][Tel],
			true);

		//CsI crystal5		  
		CsIcrystal51[Blkcor][Tel] = new G4Trd(csicrystal51[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector513[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal51[Blkcor][Tel] ,scintillator_mat ,csicrystal51[Blkcor][Tel]);

		PhysicalDetector513[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta51))*TMath::Cos((Blkphi51+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Sin(Blktheta51))*TMath::Sin((Blkphi51+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Cos(Blktheta51))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector513[Blkcor][Tel],
			csicrystal51[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal51Num[Blkcor][Tel],
			true);
	  }
	  {//2
		Silicon_521[Blkcor][Tel] = new G4Box(silicon_521[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector521[Blkcor][Tel] = new G4LogicalVolume(Silicon_521[Blkcor][Tel],Silicon_mat,silicon_521[Blkcor][Tel]);
		PhysicalDetector521[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1000.*(TMath::Sin(Blktheta52))*TMath::Cos(Blkphi52+(pi*Blkcor/2))*mm, 1000*(TMath::Sin(Blktheta52))*TMath::Sin(Blkphi52+(pi*Blkcor/2))*mm, 1000*(TMath::Cos(Blktheta52))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector521[Blkcor][Tel],//theta == 0.0316
			silicon_521[Blkcor][Tel],
			logicWorld,
			false,
			silicon_521Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_522[Blkcor][Tel] = new G4Box(silicon_522[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector522[Blkcor][Tel] = new G4LogicalVolume(Silicon_522[Blkcor][Tel] ,Silicon_mat ,silicon_522[Blkcor][Tel]);

		PhysicalDetector522[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta52))*TMath::Cos(Blkphi52+(pi*Blkcor/2))*mm, 1000.625*(TMath::Sin(Blktheta52))*TMath::Sin(Blkphi52+(pi*Blkcor/2))*mm, 1000.625*(TMath::Cos(Blktheta52))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector522[Blkcor][Tel],
			silicon_522[Blkcor][Tel],
			logicWorld,
			false,
			silicon_522Num[Blkcor][Tel],
			true);

		//CsI crystal5		  
		CsIcrystal52[Blkcor][Tel] = new G4Trd(csicrystal52[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector523[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal52[Blkcor][Tel] ,scintillator_mat ,csicrystal52[Blkcor][Tel]);

		PhysicalDetector523[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta52))*TMath::Cos(Blkphi52+(pi*Blkcor/2))*mm, 1051.058*(TMath::Sin(Blktheta52))*TMath::Sin(Blkphi52+(pi*Blkcor/2))*mm, 1051.058*(TMath::Cos(Blktheta52))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector523[Blkcor][Tel],
			csicrystal52[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal52Num[Blkcor][Tel],
			true);
	  }
	  {//3
		Silicon_531[Blkcor][Tel] = new G4Box(silicon_531[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector531[Blkcor][Tel] = new G4LogicalVolume(Silicon_531[Blkcor][Tel],Silicon_mat,silicon_531[Blkcor][Tel]);
		PhysicalDetector531[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta53))*TMath::Cos(Blkphi53+(pi*Blkcor/2))*mm, 1000*(TMath::Sin(Blktheta53))*TMath::Sin(Blkphi53+(pi*Blkcor/2))*mm, 1000*(TMath::Cos(Blktheta53))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector531[Blkcor][Tel],//theta == 0.0147
			silicon_531[Blkcor][Tel],
			logicWorld,
			false,
			silicon_531Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_532[Blkcor][Tel] = new G4Box(silicon_532[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector532[Blkcor][Tel] = new G4LogicalVolume(Silicon_532[Blkcor][Tel] ,Silicon_mat ,silicon_532[Blkcor][Tel]);

		PhysicalDetector532[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta53))*TMath::Cos(Blkphi53+(pi*Blkcor/2))*mm, 1000.625*(TMath::Sin(Blktheta53))*TMath::Sin(Blkphi53+(pi*Blkcor/2))*mm, 1000.625*(TMath::Cos(Blktheta53))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector532[Blkcor][Tel],
			silicon_532[Blkcor][Tel],
			logicWorld,
			false,
			silicon_532Num[Blkcor][Tel],
			true);

		//CsI crystal5		  
		CsIcrystal53[Blkcor][Tel] = new G4Trd(csicrystal53[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector533[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal53[Blkcor][Tel] ,scintillator_mat ,csicrystal53[Blkcor][Tel]);

		PhysicalDetector533[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta53))*TMath::Cos(Blkphi53+(pi*Blkcor/2))*mm, 1051.058*(TMath::Sin(Blktheta53))*TMath::Sin(Blkphi53+(pi*Blkcor/2))*mm, 1051.058*(TMath::Cos(Blktheta53))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector533[Blkcor][Tel],
			csicrystal53[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal53Num[Blkcor][Tel],
			true);
	  }
	  {//4
		Silicon_541[Blkcor][Tel] = new G4Box(silicon_541[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector541[Blkcor][Tel] = new G4LogicalVolume(Silicon_541[Blkcor][Tel],Silicon_mat,silicon_541[Blkcor][Tel]);
		PhysicalDetector541[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta54))*TMath::Cos((Blkphi54+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Sin(Blktheta54))*TMath::Sin((Blkphi54+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Cos(Blktheta54))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector541[Blkcor][Tel],//theta == 0.0147
			silicon_541[Blkcor][Tel],
			logicWorld,
			false,
			silicon_541Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_542[Blkcor][Tel] = new G4Box(silicon_542[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector542[Blkcor][Tel] = new G4LogicalVolume(Silicon_542[Blkcor][Tel] ,Silicon_mat ,silicon_542[Blkcor][Tel]);

		PhysicalDetector542[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta54))*TMath::Cos((Blkphi54+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Sin(Blktheta54))*TMath::Sin((Blkphi54+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Cos(Blktheta54))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector542[Blkcor][Tel],
			silicon_542[Blkcor][Tel],
			logicWorld,
			false,
			silicon_542Num[Blkcor][Tel],
			true);

		//CsI crystal5		  
		CsIcrystal54[Blkcor][Tel] = new G4Trd(csicrystal54[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector543[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal54[Blkcor][Tel] ,scintillator_mat ,csicrystal54[Blkcor][Tel]);

		PhysicalDetector543[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta54))*TMath::Cos((Blkphi54+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Sin(Blktheta54))*TMath::Sin((Blkphi54+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Cos(Blktheta54))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.14921556)*TMath::Cos(1.5042282+(pi*Blkcor/2)),TMath::Sin(0.14921556)*TMath::Sin(1.5042282+(pi*Blkcor/2)),TMath::Cos(0.14921556))),
			logicalDetector543[Blkcor][Tel],
			csicrystal54[Blkcor][Tel],
			logicWorld,
			false,
			csicrystal54Num[Blkcor][Tel],
			true);
	  }
	  //----------------------------------------------------------------------------------------------------------------------------------
	  opTubeSurface[Blkcor][Tel] = new G4OpticalSurface(TubeSurface[Blkcor][Tel]);
	  opTubeSurface[Blkcor][Tel]->SetType(dielectric_metal);
	  opTubeSurface[Blkcor][Tel]->SetFinish(polished);
	  opTubeSurface[Blkcor][Tel]->SetModel(unified);

	  new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector13[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector23[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector33[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector43[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  
	  new G4LogicalBorderSurface(TubeSurface4[Blkcor][Tel],PhysicalDetector413[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface4[Blkcor][Tel],PhysicalDetector423[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface4[Blkcor][Tel],PhysicalDetector433[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface4[Blkcor][Tel],PhysicalDetector443[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  
	  new G4LogicalBorderSurface(TubeSurface5[Blkcor][Tel],PhysicalDetector513[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface5[Blkcor][Tel],PhysicalDetector523[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface5[Blkcor][Tel],PhysicalDetector533[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface5[Blkcor][Tel],PhysicalDetector543[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);

	  const G4int num = 2;
	  G4double ephoton[num] = {2.034*eV, 4.136*eV};

	  //Optical alu rap Surface
	  G4double refractiveIndex[num] = {0.48, 1.55};
	  G4double specularLobe[num]    = {0.3, 0.3};
	  G4double specularSpike[num]   = {0.2, 0.2};
	  G4double backScatter[num]     = {0.2, 0.2};
	  G4double reflecitivity[num]   = {0.90, 0.92};

	  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

	  myST1->AddProperty("RINDEX",                ephoton, refractiveIndex, num);
	  myST1->AddProperty("REFLECITIVITY",         ephoton, reflecitivity,   num);
	  myST1->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
	  myST1->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
	  myST1->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);

	  G4cout << "Tube Surface G4MaterialPropertiesTable" << G4endl;
	  myST1->DumpTable();

	  opTubeSurface[Blkcor][Tel]->SetMaterialPropertiesTable(myST1);

	  Silicon_11VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_11VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector11[Blkcor][Tel] -> SetVisAttributes (Silicon_11VisAtt[Blkcor][Tel]);

	  Silicon_12VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_12VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector12[Blkcor][Tel] -> SetVisAttributes (Silicon_12VisAtt[Blkcor][Tel]);

	  CsIcrystal1VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal1VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector13[Blkcor][Tel] -> SetVisAttributes (CsIcrystal1VisAtt[Blkcor][Tel]);
	  //
	  Silicon_21VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_21VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector21[Blkcor][Tel] -> SetVisAttributes (Silicon_21VisAtt[Blkcor][Tel]);

	  Silicon_22VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_22VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector22[Blkcor][Tel] -> SetVisAttributes (Silicon_22VisAtt[Blkcor][Tel]);

	  CsIcrystal2VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal2VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector23[Blkcor][Tel] -> SetVisAttributes (CsIcrystal2VisAtt[Blkcor][Tel]);
	  //
	  Silicon_31VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_31VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector31[Blkcor][Tel] -> SetVisAttributes (Silicon_31VisAtt[Blkcor][Tel]);

	  Silicon_32VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_32VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector32[Blkcor][Tel] -> SetVisAttributes (Silicon_32VisAtt[Blkcor][Tel]);

	  CsIcrystal3VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal3VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector33[Blkcor][Tel] -> SetVisAttributes (CsIcrystal3VisAtt[Blkcor][Tel]);
	  //		
	  Silicon_41VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_41VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector41[Blkcor][Tel] -> SetVisAttributes (Silicon_41VisAtt[Blkcor][Tel]);

	  Silicon_42VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_42VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector42[Blkcor][Tel] -> SetVisAttributes (Silicon_42VisAtt[Blkcor][Tel]);

	  CsIcrystal4VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal4VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector43[Blkcor][Tel] -> SetVisAttributes (CsIcrystal4VisAtt[Blkcor][Tel]);
	  //////////////////////////////////////
	  Silicon_411VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_411VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector411[Blkcor][Tel] -> SetVisAttributes (Silicon_411VisAtt[Blkcor][Tel]);

	  Silicon_412VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_412VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector412[Blkcor][Tel] -> SetVisAttributes (Silicon_412VisAtt[Blkcor][Tel]);

	  CsIcrystal41VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal41VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector413[Blkcor][Tel] -> SetVisAttributes (CsIcrystal41VisAtt[Blkcor][Tel]);
	  //
	  Silicon_421VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_421VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector421[Blkcor][Tel] -> SetVisAttributes (Silicon_421VisAtt[Blkcor][Tel]);

	  Silicon_422VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_422VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector422[Blkcor][Tel] -> SetVisAttributes (Silicon_422VisAtt[Blkcor][Tel]);

	  CsIcrystal42VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal42VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector423[Blkcor][Tel] -> SetVisAttributes (CsIcrystal42VisAtt[Blkcor][Tel]);
	  //
	  Silicon_431VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_431VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector431[Blkcor][Tel] -> SetVisAttributes (Silicon_431VisAtt[Blkcor][Tel]);

	  Silicon_432VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_432VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector432[Blkcor][Tel] -> SetVisAttributes (Silicon_432VisAtt[Blkcor][Tel]);

	  CsIcrystal43VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal43VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector433[Blkcor][Tel] -> SetVisAttributes (CsIcrystal43VisAtt[Blkcor][Tel]);
	  //		
	  Silicon_441VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_441VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector441[Blkcor][Tel] -> SetVisAttributes (Silicon_441VisAtt[Blkcor][Tel]);

	  Silicon_442VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_442VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector442[Blkcor][Tel] -> SetVisAttributes (Silicon_442VisAtt[Blkcor][Tel]);

	  CsIcrystal44VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal44VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector443[Blkcor][Tel] -> SetVisAttributes (CsIcrystal44VisAtt[Blkcor][Tel]);
	  //////////////////////////////
	  Silicon_511VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_511VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector511[Blkcor][Tel] -> SetVisAttributes (Silicon_511VisAtt[Blkcor][Tel]);

	  Silicon_512VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_512VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector512[Blkcor][Tel] -> SetVisAttributes (Silicon_512VisAtt[Blkcor][Tel]);

	  CsIcrystal51VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal51VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector513[Blkcor][Tel] -> SetVisAttributes (CsIcrystal51VisAtt[Blkcor][Tel]);
	  //
	  Silicon_521VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_521VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector521[Blkcor][Tel] -> SetVisAttributes (Silicon_521VisAtt[Blkcor][Tel]);

	  Silicon_522VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_522VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector522[Blkcor][Tel] -> SetVisAttributes (Silicon_522VisAtt[Blkcor][Tel]);

	  CsIcrystal52VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal52VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector523[Blkcor][Tel] -> SetVisAttributes (CsIcrystal52VisAtt[Blkcor][Tel]);
	  //
	  Silicon_531VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_531VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector531[Blkcor][Tel] -> SetVisAttributes (Silicon_531VisAtt[Blkcor][Tel]);

	  Silicon_532VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_532VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector532[Blkcor][Tel] -> SetVisAttributes (Silicon_532VisAtt[Blkcor][Tel]);

	  CsIcrystal53VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal53VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector533[Blkcor][Tel] -> SetVisAttributes (CsIcrystal53VisAtt[Blkcor][Tel]);
	  //		
	  Silicon_541VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_541VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector541[Blkcor][Tel] -> SetVisAttributes (Silicon_541VisAtt[Blkcor][Tel]);

	  Silicon_542VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_542VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector542[Blkcor][Tel] -> SetVisAttributes (Silicon_542VisAtt[Blkcor][Tel]);

	  CsIcrystal54VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal54VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector543[Blkcor][Tel] -> SetVisAttributes (CsIcrystal54VisAtt[Blkcor][Tel]);
	}
  }
  //G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
  //if (opticalSurface) opticalSurface->DumpInfo();
  // Generate & Add Material Properties Table attached to the optical surfaces
  //runManager -> SetSensitiveDetector(pvp);
  return physWorld;
}
