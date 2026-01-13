
#include "StepAction.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <numeric>
#include "Randomize.hh"
#include "G4EventManager.hh"
/*
StepAction::StepAction(G4String data_file, double lunghezza_collimatore)			 
{
 datai_file = data_file;
  collimatore = lunghezza_collimatore*mm;
  for(int n = 0; n < 33; n++ )
  {
    s_x1[n] = 0;
    s_x2[n] = 0;  
    s_y1[n] = 0;  
    s_y2[n] = 0;    
  }
}*/


StepAction::StepAction(MyEventAction *eventAction)
{
    fEventAction = eventAction;
}

StepAction::~StepAction()
{
  /*
  unsigned int n;
  std::ofstream data_scrivi;
  data_scrivi.open(datai_file.c_str(),std::fstream::in | std::fstream::out | std::fstream::app);
  for(int n = 0; n < 33; n++)
  {
    data_scrivi << s_x1[n] << " ";
    //G4cout << " ciccio " << s_x1[n]  << G4endl;
  }
  for(int n = 0; n < 33; n++ )
  {
    data_scrivi << s_x2[n] << " ";
    //G4cout << " ciccio " << s_x2[n]  << G4endl;
  }
  for(int n = 0; n < 33; n++ )
  {
    data_scrivi << s_y1[n] << " ";
    //G4cout << " ciccio " << s_y1[n]  << G4endl;
  }
  for(int n = 0; n < 33; n++ )
  {
    data_scrivi << s_y2[n] << " ";
    //G4cout << " ciccio " << s_y2[n]  << G4endl;
  }
  data_scrivi << " 9999 " << std::endl;
  data_scrivi.close(); */
}

void StepAction::UserSteppingAction(const G4Step* aStep)
{

  //G4cout << aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << G4endl;
if ( aStep->GetTrack()->GetDefinition()->GetParticleName() == "proton" && aStep->GetPostStepPoint()->GetProcessDefinedStep() != nullptr
    /* aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()== "PPAC_Vol" &&
     aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "expHall"*/)   
  {   

      G4AnalysisManager *man = G4AnalysisManager::Instance();   

      const DC* detectorConstruction = static_cast<const DC*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
      G4LogicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    
      G4LogicalVolume* fScoringVolume_1 = detectorConstruction->GetScoringVolume();
      G4Track *track = aStep->GetTrack();
      G4StepPoint *preStepPoint; 
      G4StepPoint *postStepPoint ;   
      if (volume != fScoringVolume_1) return;          

       
        if (aStep->IsFirstStepInVolume()){
           preStepPoint = aStep->GetPreStepPoint();
           postStepPoint = aStep->GetPostStepPoint();
           
           G4ThreeVector posPhoton2 = postStepPoint->GetPosition()/cm;

           //man->FillH2(1, posPhoton2[0], posPhoton2[1]);
          

           


        }

/*      
      //G4cout << aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << G4endl;
      G4float tg = aStep->GetPostStepPoint()->GetGlobalTime()/(0.001*ns);
      //G4float tl = aStep->GetPostStepPoint()->GetLocalTime()/(0.001*ns);
      G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
      G4float energia = (aStep->GetPostStepPoint()->GetKineticEnergy()/eV);
      G4float wave = 1240/energia;

      G4double edep = aStep->GetTotalEnergyDeposit();
      fEventAction->AddEdep(edep);     


      G4double X = position.getX();
      G4double Y = position.getY();
      G4cout << "X= " << X << " Y= "<< Y/mm << G4endl;
      double r_rand = G4UniformRand();
      int indice;

      man->FillNtupleDColumn(0, 0, X/mm);
      man->AddNtupleRow(0);

      man->FillNtupleDColumn(1, 0, Y/mm);
  
      man->AddNtupleRow(1);

      man->FillNtupleDColumn(2, 0, wave);
      man->AddNtupleRow(2);
 
      man->FillNtupleDColumn(3, 0, tg);
      man->AddNtupleRow(3);


      man->FillH2(0, X/mm, Y/mm);



      if (tg<7000) 
      { 
        if ((X==50*mm+collimatore) && Y>-50*mm &&  Y<50*mm) 
        {
    	  if (r_rand<=0.3)
	  {
            indice = hysto(Y);
            s_y1[indice]++;
          }
        }
        if ((X==-50*mm-collimatore) && Y>-50*mm &&  Y<50*mm) 
        {
    	  if (r_rand<=0.3)
	  {
            indice = hysto(Y);
            s_y2[indice]++;
          }
        }
        if ((Y==50*mm+collimatore) && X>-50*mm && X<50*mm)
        {
    	  if (r_rand<=0.3)
	  {
            indice = hysto(X);
            s_x1[indice]++;
          }
        }
        if ((Y==-50*mm-collimatore) && X>-50*mm &&  X<50*mm)      
        {
    	  if (r_rand<=0.3)
	  {
            indice = hysto(X);
            s_x2[indice]++;
          }
        }
      }*/
      //G4float tempo = tg-tl;
      //tot_photons++;
      //G4cout << "X= " << X/mm << " Y= "<< Y/mm << " Wavelenght= " << wave << G4endl;
      //mytempo.open(nome_data_file_t.c_str(),fstream::in | fstream::out | fstream::app);
  }
}

double StepAction::hysto(double valore)
{

  double calcola = 0;
  double segno;
  if (valore > 0 ) segno=1;
  if (valore < 0 ) segno=-1;
  calcola = (int) (fabs(valore)+0.9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999)/3;
  //G4cout << fabs(valore) << " " << calcola*segno << G4endl;
  calcola = (calcola*segno)+16;
  return calcola;
}
