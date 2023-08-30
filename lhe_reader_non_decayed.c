//&>/dev/null;x="${0%.*}";[ ! "$x" -ot "$0" ]||(rm -f "$x";g++ -o "$x" "$0" -I`root-config --incdir` `root-config --libs`);exit

// Build: g++ lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
// OR, g++ -std=c++11 lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
//#include </home/muhammad/root6/build/include/TLorentzVector.h>
#include <TLorentzVector.h>
using namespace std;

// Pour une description du format leshouches
// hep-ph/0609017
//
// pour chaque evenement
// une ligne générale : NbPart idprocess poids scale alpha_em alpha_s
// pour chaque particule : id status mere1 mere2 couleur1 couleur2 px py pz E m lifetime spin  
int main(int argc, char **argv) {

  if (argc != 2) {
    cout << "Missing argument. Expected LHE filename without '.lhe'"<< endl;
    exit(-5);
  }
  float DeltaR(float eta1, float eta2, float phi1, float phi2);
  float DeltaPhi(float phi1, float phi2);

  string basename=argv[1];

  string lhefname = basename+".lhe";
  string rootfname = basename+".root";

  string tt;
  int event=0;
  int npart,idprocess;
  double weight,scale,alpha_em,alpha_s;

  TH1::SetDefaultSumw2(true);
  
//--------------------- insert here histogram definition -----------------------------------//
//
  TH1F* hPt_tbq = new TH1F("tbq_pt","top_bquarks pt",50,0.0,400) ;
  TH1F* hEta_tbq = new TH1F("tbq_eta","top_bquarks eta",50,-6,6) ;
  TH1F* hPhi_tbq = new TH1F("tbq_phi","top_bquarks phi",50,-3.16,3.16) ;
  
  TH1F* hPt_bq = new TH1F("bq_pt","bqurk pt",50,0.0,400) ;
  TH1F* hEta_bq = new TH1F("bq_eta","bquark eta",50,-6,6) ;
  TH1F* hPhi_bq = new TH1F("bq_phi","bquark phi",50,-3.16,3.16) ;
  
  TH1F* hPt_lep = new TH1F("lep_pt","lepton pt",50,0.0,400) ;
  TH1F* hEta_lep = new TH1F("lep_eta","lepton eta",50,-6,6) ;
  TH1F* hPhi_lep = new TH1F("lep_phi","lepton phi",50,-3.16,3.16) ;

  TH1F* hPt_wboson = new TH1F("wboson_pt","wboson pt",50,0.0,400) ;
  TH1F* hEta_wboson = new TH1F("wboson_eta","wboson eta",50,-6,6) ;
  TH1F* hPhi_wboson = new TH1F("wboson_phi","wboson phi",50,-3.16,3.16) ;
  TH1F* hM_wboson = new TH1F("hM_wboson", "wboson mass", 50, 60, 100);
  
  TH1F* hPt_top = new TH1F("top_pt","top pt",50,0.0,400) ;
  TH1F* hEta_top = new TH1F("top_eta","top eta",50,-6,6) ;
  TH1F* hPhi_top = new TH1F("top_phi","top phi",50,-3.16,3.16) ;
  TH1F* hM_top = new TH1F("hM_top", "topquark mass", 50, 150, 200);
  
  TH1F* hPt_higgs = new TH1F("higgs_pt","Higgs pt",50,0.0,400) ;
  TH1F* hEta_higgs = new TH1F("higgs_eta","Higgs eta",50,-6,6) ;
  TH1F* hPhi_higgs = new TH1F("higgs_phi","Higgs phi",50,-3.16,3.16) ;
  TH1F* hM_higgs = new TH1F("hM_higgs", "Higgs mass", 50, 100, 150);
  
//---------------------------- end histogram definition -------------------------------------//
//
cout<<"I am here: "<<endl;

  int nlept=0, nsemi=0, nhadr=0;
  ifstream ff(lhefname.c_str(),ios::in); //ouverture du fichier .lhe
  //ifstream ff("test.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/s1/madevent/Events/zp4000_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/QCD/madevent/Events/qcd_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  int negativeWeight = 0;
  long long line = 0;
  while(!ff.eof()) {
    std::stringstream buffer;
    ff>>tt;
    buffer << tt << std::endl;
    line++;

    if(tt=="<event>") {
      ff>>npart>>idprocess>>weight>>scale>>alpha_em>>alpha_s; //event definition
      buffer << npart << " " << idprocess << " " << weight << " " << scale << " " << alpha_em << " " << alpha_s << std::endl;
      line++;
      event++;
      if (weight < 0) {
        negativeWeight++;
        weight = -1;
      } else {
        weight = 1;
      }
      /*weight = 1.;*/

      if (event%1==0) cout << "reading event "<< event << endl;
      int lmin=-1, lmax=-1, bmin=-1, met1=-1, met2=-1, bmax=-1, met=-1, topLep,topbJet, tbarLep, tbarbJet, top_jet, tbar_jet;
      int Part_Id, Moth_Id, Part_absId, Moth_absId;
int n_lep=0, n_alep=0, n_lep_nu=0, n_topjets=0, n_tbarjets=0, n_topbjets=0, n_bJets=0, n_tbar_nu=0, n_top_nu=0, n_wlep=0;
int n_aele=0, n_aele_nu=0, n_amuon=0, n_amuon_nu=0;

int  bj[2]={-1,-1};
int lp[2]={-1,-1};

      int q[4]={-1,-1,-1,-1,};
      int muon[5]={-1,-1,-1,-1,-1};
      int elec[5]={-1,-1,-1,-1,-1};
      int elec_nu[5]={-1,-1,-1,-1,-1};
      int lep[3]={-1,-1,-1};
      int muon_nu[5]={-1,-1,-1,-1,-1};
      int bq[4]={-1,-1,-1,-1};
      
      int top=-1,topbar=-1,zprime=-1;
      int *Id      = new int[npart+1];
      int *Status  = new int[npart+1];
      int *Mother1 = new int[npart+1];
      int *Mother2 = new int[npart+1];
      int *Color1  = new int[npart+1];
      int *Color2  = new int[npart+1];
      double *px = new double[npart+1];
      double *py = new double[npart+1];
      double *pz = new double[npart+1];
      double *E = new double[npart+1];
      double *m = new double[npart+1];
      double *lifetime = new double[npart+1];
      double *spin = new double[npart+1];
      
      TLorentzVector **v = new TLorentzVector*[npart+1];
      TLorentzVector v_muon1, v_muon2, v_muon3, v_muon4, v_muon5, v_mu_nu1, v_mu_nu2, v_mu_nu3, v_mu_nu4, v_mu_nu5;
      TLorentzVector v_elec1, v_elec2, v_elec3, v_elec4, v_elec5, v_elec_nu1, v_elec_nu2, v_elec_nu3, v_elec_nu4, v_elec_nu5;
      TLorentzVector v_lep1, v_lep2, v_lep3, v_lep4, v_lep5;
      TLorentzVector v_wboson, v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4, v_w_boson5, v_hboson;
      TLorentzVector v_bq, v_hbq, v_hbbarq, v_bq1, v_bq2, v_bq3, v_bq4;
      TLorentzVector v_top, v_top1, v_top2, v_top3, v_top4, v_top5;
      TLorentzVector v_top_lep, v_top_alep, v_lep_nu, v_top_bJet, v_tbar_lep, v_top_jet, v_tbar_jet, v_tbar_nu, v_top_nu, v_wminus_lep;
      TLorentzVector v_aele, v_aele_nu, v_amuon, v_amuon_nu, v_alep, v_alep_nu, v_lep, v_lep_anu; 
      
      // in lhe first line is number 1, so fill unused array [0] with a crazy value;
      Id[0]= -99999;
      Status[0]= -99999;
      Mother1[0]= -99999;
      Mother2[0]= -99999;
      Color1[0]= -99999;
      Color2[0]= -99999;
      px[0]= -99999;
      py[0]= -99999;
      pz[0]= -99999;
      E[0]= -99999;
      m[0]= -99999;
      lifetime[0]= -99999;
      spin[0]= -99999;
     for (int i=1 ; i<npart+1 ; i++) { //start at one
        ff >> Id[i] >> Status[i] >> Mother1[i] >> Mother2[i] >> Color1[i] >> Color2[i]
           >> px[i] >> py[i] >> pz[i] >> E[i] >> m[i] >> lifetime[i] >> spin[i] ;
        buffer << Id[i] << " " << Status[i] << " " << std::endl;
        line++;
        v[i] = new TLorentzVector(px[i], py[i], pz[i], E[i]);
        if (Status[i]==-1) continue; // status -1 = initial quark ==> skip
        if (Id[i]==6)  top=i;
        if (Id[i]==-6) topbar=i;
        if (Id[i]>6000000) zprime=i;


//------------------top and tbar lep, bJets-----------------------------------

	//b qaurks:
	if (abs(Id[Mother1[i]])==6)
	{ //Mother top
	  if (abs(Id[i])==5)
	  { // bquarks
		v_bq.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
           }
	}
	
	if (Id[Mother1[i]] == 25)
	{
		if (Id[i]==5)
		{
		v_hbq.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		}
		if (Id[i]==-5)
		{
		v_hbbarq.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		}
		v_hboson = v_hbq + v_hbbarq;
		
	}
	
	//electrons, muons, and neutrinos:
	if ( Id[Mother1[i]] == 6 || Id[Mother1[i]] == 24)
	{
          if ( abs(Id[i])==-11 ) 
          { // charged electrons
          	v_aele.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		n_aele++;
	  }
	  if ( abs(Id[i])==12 ) 
          { // charged electrons neutrinos
          	v_aele_nu.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		n_aele_nu++;
	  }
	  if ( abs(Id[i])==-13 ) 
          { // charged muons
          	v_amuon.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		n_amuon++;
	  }
	  if ( abs(Id[i])==14 ) 
          { // charged muons neutrinos
          	v_amuon_nu.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		n_amuon_nu++;
	  }
	  }
	  
	  
	  //leptons anti leptons:
	  if (Id[Mother1[i]]==-24)
	  {
	  if  (Id[i] == 11 || Id[i] == 13) 
          { // charged leptons
		v_lep.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		hPt_lep->Fill(v_lep.Pt());
		hEta_lep->Fill(v_lep.Eta());
		hPhi_lep->Fill(v_lep.Phi());
		n_lep++;
	  }
	  if (Id[i]== -12 || Id[i] == -14)
	  	{
	  	v_lep_anu.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
	  	}
	  }
	  
	  if (Id[Mother1[i]] == 24)
	  {
	  if  (Id[i] == -11 || Id[i] == -13) 
          { // charged leptons
		v_alep.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		n_alep++;
	  }
	  if (Id[i]== 12 || Id[i] == 14)
	  	{
	  	v_alep_nu.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
	  	}
	  	
	  	v_wboson = v_alep + v_alep_nu;
	  	v_top = v_wboson + v_bq;
	  }
	  
} //loop of i

//--------------------------------Filling Histograms----------------------------------//

hPt_bq->Fill(v_hbq.Pt());
hEta_bq->Fill(v_hbq.Eta());
hPhi_bq->Fill(v_hbq.Phi());

hPt_tbq->Fill(v_bq.Pt());
hEta_tbq->Fill(v_bq.Eta());
hPhi_tbq->Fill(v_bq.Phi());

hPt_wboson->Fill(v_wboson.Pt());
hEta_wboson->Fill(v_wboson.Eta());
hPhi_wboson->Fill(v_wboson.Phi());
hM_wboson->Fill(v_wboson.M());

hPt_top->Fill(v_top.Pt());
hEta_top->Fill(v_top.Eta());
hPhi_top->Fill(v_top.Phi());
hM_top->Fill(v_top.M());

hPt_higgs->Fill(v_wboson.Pt());
hEta_higgs->Fill(v_wboson.Eta());
hPhi_higgs->Fill(v_wboson.Phi());
hM_higgs->Fill(v_hboson.M());

// --- end filling the histograms

      ff>>tt;
      line++;
      //if (event==100)  break;
      delete Id;
      delete Status;
      delete Mother1;
      delete Mother2;
      delete Color1;
      delete Color2;
      delete px;
      delete py;
      delete pz;
      delete E;
      delete m;
      delete lifetime;
      delete spin;
      for (int k=1 ; k<npart+1 ; delete v[k++]);

    }
}

  cout << " Total number of events --> " << event << endl;
  TFile *rootfile = new TFile(rootfname.c_str(),"recreate");
  
  hPt_lep->Write();
  hEta_lep->Write();
  hPhi_lep->Write();
  
  hPt_bq->Write();
  hEta_bq->Write();
  hPhi_bq->Write();
  
  hPt_tbq->Write();
  hEta_tbq->Write();
  hPhi_tbq->Write();
  
  hPt_wboson->Write();
  hEta_wboson->Write();
  hPhi_wboson->Write();
  hM_wboson->Write();
  
  hPt_top->Write();
  hEta_top->Write();
  hPhi_top->Write();
  hM_top->Write();
  
  hPt_higgs->Write();
  hEta_higgs->Write();
  hPhi_higgs->Write();
  hM_higgs->Write();

  rootfile->Close();

  cout << "Events with negative weight: " << negativeWeight << endl;
  //cout << "lept decay = " << nlept*1.0/event << endl;
  //cout << "hadr decay = " << nhadr*1.0/event << endl;
  //cout << "semi decay = " << nsemi*1.0/event << endl;
  exit(0);

}
float DeltaPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  return dphi;
}
float DeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
  return DELTAR;
}

