


/* 

Copyright (c) 2007-2012, The Regents of the University of California. 
Produced at the Lawrence Livermore National Laboratory 
UCRL-CODE-227323. 
All rights reserved. 
 
For details, see http://nuclear.llnl.gov/simulations
Please also read this http://nuclear.llnl.gov/simulations/additional_bsd.html
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
 
1.  Redistributions of source code must retain the above copyright
notice, this list of conditions and the disclaimer below.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the disclaimer (as noted below) in
the documentation and/or other materials provided with the
distribution.

3. Neither the name of the UC/LLNL nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OF
THE UNIVERSITY OF CALIFORNIA, THE U.S. DEPARTMENT OF ENERGY OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "CRYGenerator.h"
#include "CRYSetup.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <math.h>
#include <stdlib.h> 
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraph.h" 
#include "TProfile.h"// For Ubuntu Linux

int main( int argc, const char *argv[]) {

  if (!TROOT::Initialized()) {
    static TROOT root("RooTuple", "RooTuple ROOT God in CosmicRaY simulation");
  }


  if ( argc < 2 ) {
    std::cout << "usage " << argv[0] << " <setup file name> <N events>\n";
    std::cout << "N events = 10k by default\n";
    return 0;
  }

  TFile *outputFile=new TFile("test_05.root","RECREATE");
  

  /*TH1F* xall         = new TH1F("xall",   "xall",   xybins, -box_size, box_size);
  TH1F* yall         = new TH1F("yall",   "yall",   xybins, -box_size, box_size);*/
 // TH1F* xmuon        = new TH1F("xmuon",  "xmuon",  xybins, -box_size, box_size);
  //TH1F* ymuon        = new TH1F("ymuon",  "ymuon",  xybins, -box_size, box_size);
  //TH2F* xyall        = new TH2F("xyall",  "xyall",  xybins, -box_size, box_size, xybins, -box_size, box_size);
  /*TH2F* xymuons      = new TH2F("xymuon", "xymuon", xybins, -box_size, box_size, xybins, -box_size, box_size);


  TH1F* costhe       = new TH1F("costhe",  "costhe plot;cos(theta)",   200, -1., 0. );

  std::map<CRYParticle::CRYId,TH1F*> keHistos;

  
    // Fill the histograms
    multiplicity->Fill(ev->size());
    kePrimary->Fill(log10(gen.primaryParticle()->ke());


 
}
      latLoc->Fill(sqrt(part->x()*part->x()+part->y()*part->y())); 
     // chargeHist->Fill(part->id(),part->charge());
      //multHist->Fill(part->id());
      //xall->Fill(part->x());
      //yall->Fill(part->y());
      //xyall->Fill(part->x(), part->y());
      costhe->Fill(part->w());*/


// Define log scale
  Int_t nbins = 9*5;
  Double_t xmin=1e0;
  Double_t xmax=1e9;
  Double_t *xbins    = new Double_t[nbins+1];
  Double_t xlogmin = log10(xmin);
  Double_t xlogmax = log10(xmax);
  Double_t dlogx   = (xlogmax-xlogmin)/((Double_t)nbins);
  for (int i=0;i<=nbins;i++) { 
            Double_t xlog = xlogmin+ i*dlogx;
            xbins[i] = exp( log(10) * xlog ); 
  }

 
  TH1I* multiplicity = new TH1I("mult",       "mult",      100,     0,  100);
  TH1F* latLoc       = new TH1F("latLoc",     "latLoc",    100, -150.0, 150.0);
  TH1F* kePrimary    = new TH1F("kePrimary",  "kePrimary",   nbins,xbins);
  TH1F* chargeHist   = new TH1F("chargeHist", "chargeHist",   10, -0.5, 9.5);
  TH1F* multHist     = new TH1F("multHist",   "multHist",     10, -02.5, 9.5);
  TH1F* time         = new TH1F("time",   "time",           5*3100, -0.5, 3100.5);
  
  //TH1F* keMuon       = new TH1F("keMuon",     "keMuon",    5*16, 0.0, 1e4);
  //TH1F* keNeutron    = new TH1F("keNeutron",  "keNeutron",  5*16, 0.0, 1e4);
  //TH1F* kePion       = new TH1F("kePion",     "kePion",     5*16, 0.0, 1e4);
  //TH1F* keKaon       = new TH1F("keKaon",     "keKaon",     5*16, 0.0, 1e4);
  //TH1F* keGamma      = new TH1F("keGamma",    "keGamma",    5*16, 0.0, 1e4);
  //TH1F* keElectron   = new TH1F("keElectron", "keElectron", 5*16, 0.0, 1e4);
  //TH1F* keProton     = new TH1F("keProton",   "keProton",   5*16, 0.0, 1e4);
  TH1F* keMuon       = new TH1F("keMuon",     "keMuon",   nbins,xbins);
  TH1F* keMuonp      = new TH1F("keMuonp",     "keMuonp", nbins,xbins);
  TH1F* keMuonn      = new TH1F("keMuonn",     "keMuonn", nbins,xbins);

  TH1F* keNeutron    = new TH1F("keNeutron",  "keNeutron",nbins,xbins);

  TH1F* kePion       = new TH1F("kePion",     "kePion",     nbins,xbins);
  TH1F* kePionP       = new TH1F("kePionP",     "kePionP",     nbins,xbins);
  TH1F* kePion0       = new TH1F("kePion0",     "kePion0",     nbins,xbins);
  TH1F* kePionN       = new TH1F("kePionN",     "kePionN",     nbins,xbins);

  TH1F* keKaon       = new TH1F("keKaon",     "keKaon",     nbins,xbins);
  TH1F* keKaonP      = new TH1F("keKaonP",     "keKaonP",     nbins,xbins);  
  TH1F* keKaon0      = new TH1F("keKaon0",     "keKaon0",     nbins,xbins);
  TH1F* keKaonN      = new TH1F("keKaonN",     "keKaonN",     nbins,xbins);

  TH1F* keGamma      = new TH1F("keGamma",    "keGamma",    nbins,xbins);

  TH1F* keElectron   = new TH1F("keElectron", "keElectron", nbins,xbins);
  TH1F* kePositron   = new TH1F("kePositron", "kePositron", nbins,xbins);
  TH1F* keElectronN   = new TH1F("keElectronN", "keElectronN", nbins,xbins);

  TH1F* keProton     = new TH1F("keProton",   "keProton",  nbins,xbins);
  TH1F* keProtonp     = new TH1F("keProtonp",   "keProtonp",  nbins,xbins);
  TH1F* keAntiProton      = new TH1F(" keAntiProton",     " keAntiProton", nbins,xbins);
   
  TH1F* keMuon5      = new TH1F("keMuon5",       "keMuon5", nbins,xbins);
  TH1F* keMuon15     = new TH1F("keMuon15",     "keMuon15", nbins,xbins);
  TH1F* keMuon25     = new TH1F("keMuon25",     "keMuon25", nbins,xbins);
  TH1F* keMuon35     = new TH1F("keMuon35",     "keMuon35", nbins,xbins);
  TH1F* keMuon45     = new TH1F("keMuon45",     "keMuon45", nbins,xbins);
  TH1F* keMuon55     = new TH1F("keMuon55",     "keMuon55", nbins,xbins);
  TH1F* keMuon65     = new TH1F("keMuon65",     "keMuon65", nbins,xbins);
  TH1F* keMuon75     = new TH1F("keMuon75",     "keMuon75", nbins,xbins);
  TH1F* keMuon85     = new TH1F("keMuon85",     "keMuon85", nbins,xbins);
  
 // TH2F* ratio     = new TH2F("ratio",     "ratio", nbins,xbins,100,1,2);
  auto ratio  = new TProfile("ratio","ratio", nbins,xbins,1,2," ");
  TH1F* PMuon     = new TH1F("PMuon",     "PMuon", nbins,xbins);
  TH1F* PMuonp     = new TH1F("PMuonp",     "PMuonp", nbins,xbins);
  TH1F* PMuonn     = new TH1F("PMuonn",     "PMuonn", nbins,xbins);
//histo of muon for different energy of primary
  TH1F* keMuonPrm1     = new TH1F("keMuonPrm1",     "keMuonPrm1", nbins,xbins);
  TH1F* keMuonPrm2     = new TH1F("keMuonPrm2",     "keMuonPrm2", nbins,xbins);
  TH1F* keMuonPrm3     = new TH1F("keMuonPrm3",     "keMuonPrm3", nbins,xbins);
  TH1F* keMuonPrm4     = new TH1F("keMuonPrm4",     "keMuonPrm4", nbins,xbins);
  TH1F* keMuonPrm5     = new TH1F("keMuonPrm5",     "keMuonPrm5", nbins,xbins);
  
  TH1F* keMuonPrm1p     = new TH1F("keMuonPrm1p",     "keMuonPrm1p", nbins,xbins);
  TH1F* keMuonPrm2p     = new TH1F("keMuonPrm2p",     "keMuonPrm2p", nbins,xbins);
  TH1F* keMuonPrm3p     = new TH1F("keMuonPrm3p",     "keMuonPrm3p", nbins,xbins);
  TH1F* keMuonPrm4p     = new TH1F("keMuonPrm4p",     "keMuonPrm4p", nbins,xbins);
  TH1F* keMuonPrm5p     = new TH1F("keMuonPrm5p",     "keMuonPrm5p", nbins,xbins);
  
  TH1F* keMuonPrm1n     = new TH1F("keMuonPrm1n",     "keMuonPrm1n", nbins,xbins);
  TH1F* keMuonPrm2n     = new TH1F("keMuonPrm2n",     "keMuonPrm2n", nbins,xbins);
  TH1F* keMuonPrm3n     = new TH1F("keMuonPrm3n",     "keMuonPrm3n", nbins,xbins);
  TH1F* keMuonPrm4n     = new TH1F("keMuonPrm4n",     "keMuonPrm4n", nbins,xbins);
  TH1F* keMuonPrm5n     = new TH1F("keMuonPrm5n",     "keMuonPrm5n", nbins,xbins);
  

  float box_size =1.5;
  int xybins = 150;
  TH1F* xall         = new TH1F("xall",   "xall",   xybins, -box_size, box_size);
  TH1F* yall         = new TH1F("yall",   "yall",   xybins, -box_size, box_size);
  TH1F* xmuon        = new TH1F("xmuon",  "xmuon",  xybins, -box_size, box_size);
  TH1F* ymuon        = new TH1F("ymuon",  "ymuon",  xybins, -box_size, box_size);
  TH2F* xyall        = new TH2F("xyall",  "xyall",  xybins, -box_size, box_size, xybins, -box_size, box_size);
  TH2F* xymuons      = new TH2F("xymuon", "xymuon", xybins, -box_size, box_size, xybins, -box_size, box_size);

  TH1F* costhe       = new TH1F("costhe",  "costhe plot;cos(theta)",   200, -1., 0. );
  TH1F* costhemu      = new TH1F("costhemu",  "costhemu",   200, 0, 1 );
 
//distribution angulaire different primary

  TH1F* costhemuPrm1      = new TH1F("costhemuPrm1",  "costhemuPrm1",   nbins, 0, 1);
  TH1F* costhemuPrm2      = new TH1F("costhemuPrm2",  "costhemuPrm2",   nbins, 0, 1 );
  TH1F* costhemuPrm3      = new TH1F("costhemuPrm3",  "costhemuPrm3",   nbins, 0, 1 );
  TH1F* costhemuPrm4      = new TH1F("costhemuPrm4",  "costhemuPrm4",   nbins, 0, 1 );
  TH1F* costhemuPrm5      = new TH1F("costhemuPrm5",  "costhemuPrm5",   nbins, 0, 1 );
  
  TH1F* costhePrimary       = new TH1F("costhePrimary",  "costhePrimary",   nbins, -1,1 );


   TH1F* histo       = new TH1F("histo",  "histo",   nbins, 0,1 );
   TH1F* histoPrm1       = new TH1F("histoPrm1",  "histoPrm1",   nbins, 0,1 );
   TH1F* histoPrm2       = new TH1F("histoPrm2",  "histoPrm2",   nbins, 0,1 );
   TH1F* histoPrm3       = new TH1F("histoPrm3",  "histoPrm3",   nbins, 0,1 );
   TH1F* histoPrm4      = new TH1F("histoPrm4",  "histoPrm4",   nbins, 0,1 );
   TH1F* histoPrm5       = new TH1F("histoPrm5",  "histoPrm5",   nbins, 0,1 );
  std::map<CRYParticle::CRYId,TH1F*> keHistos;

  keHistos[CRYParticle::Neutron]=keNeutron;
  keHistos[CRYParticle::Muon]=keMuon;
  keHistos[CRYParticle::Proton]=keProton;
  keHistos[CRYParticle::Electron]=keElectron;
  keHistos[CRYParticle::Gamma]=keGamma;
  keHistos[CRYParticle::Kaon]=keKaon;
  keHistos[CRYParticle::Pion]=kePion;
  //keHistos[CRYParticle::Positron]=kePositron;
  int nEv = 100000;
  if (argc > 2 ) nEv = atoi(argv[2]);

  // Read the setup file into setupString
  std::ifstream inputFile;
  inputFile.open(argv[1],std::ios::in);
  char buffer[1000];

  std::string setupString("");
  while (!inputFile.getline(buffer,1000).eof()) {
    setupString.append(buffer);
    setupString.append(" ");
  }

  // Parse the contents of the setup file
  CRYSetup *setup=new CRYSetup(setupString,"./data");
  double detectorsize = setup->param(CRYSetup::subboxLength);

  // Setup the CRY event generator
  CRYGenerator gen(setup);
 
  // Generate the events
  int ibin;
  Double_t wbin;




  int nMuon = 0;
  int nMuonp=0;
  double timeFlux =0.0;
  double timeFlux0 =0.0;
  

  std::vector<CRYParticle*> *ev=new std::vector<CRYParticle*>;
  for (int i = 0; i < nEv; i++) {
    ev->clear();
    gen.genEvent(ev);

    std::cout << "Primary " <<    " " <<  CRYUtils::partName(gen.primaryParticle()->id()) 
          << " ke=" << gen.primaryParticle()->ke() 
   	      << "  w=" << gen.primaryParticle()->w() << "\n"<< std::endl;
    //
    //if (i % 100000 == 0)
      // std::cout << "Event: " << i << "/" << nEv << "  " << 1.0*i/nEv*100 << " %" << std::endl;
      
    
    // Fill the histograms
    multiplicity->Fill(ev->size());
   
    //kePrimary->Fill(gen.primaryParticle()->ke());
       ibin = kePrimary->FindBin(gen.primaryParticle()->ke());
       wbin = 1.0/kePrimary->GetBinWidth(ibin);

      kePrimary->Fill(gen.primaryParticle()->ke());
      costhePrimary->Fill(gen.primaryParticle()->w()); 


    for (unsigned j = 0; j < ev->size(); j++) {
      CRYParticle *part=(*ev)[j];
      
      //....printout all secondaries every 1000 events just for fun
      //if (i % 1000 == 0) {
      //  std::cout << "Secondary " << j << " " <<
      //	  CRYUtils::partName(part->id()) << " ke=" << part->ke() << " " << log10( part->ke()) << "\n";
      //}
      
      //keHistos[part->id()]->Fill(log10(part->ke()));
     // ibin = keHistos[part->id()]->FindBin(part->ke());
      //wbin = 1.0/keHistos[part->id()]->GetBinWidth(ibin);
      //keHistos[part->id()]->Fill(part->ke(),wbin);

      latLoc->Fill(sqrt(part->x()*part->x()+part->y()*part->y())); 
      chargeHist->Fill(part->id(),part->charge());
      multHist->Fill(part->id());
      xall->Fill(part->x());
      yall->Fill(part->y());
      xyall->Fill(part->x(), part->y());
      costhe->Fill(part->w());
      time->Fill(part->t());

  if (part->id() == CRYParticle::Gamma) {
	ibin = keGamma->FindBin(part->ke());
	wbin = 1.0/keGamma->GetBinWidth(ibin);     
	keGamma->Fill(part->ke(),wbin);
      }







  if (part->id() == CRYParticle::Proton){ 

        ibin = keProton->FindBin(part->ke());
	wbin = 1.0/keProton->GetBinWidth(ibin);     
	keProton->Fill(part->ke(),wbin);
	if (part->charge() > 0)  keProtonp->Fill(part->ke(),wbin);
	if (part->charge() < 0)   keAntiProton->Fill(part->ke(),wbin);

        
    }

 if (part->id() == CRYParticle::Muon) {
	ibin = keMuon->FindBin(part->ke());
	wbin = 1.0/keMuon->GetBinWidth(ibin);     
	keMuon->Fill(part->ke(),wbin);
	if (part->charge() > 0)  keMuonp->Fill(part->ke(),wbin);
	if (part->charge() < 0)  keMuonn->Fill(part->ke(),wbin);
        xmuon->Fill(part->x());
        ymuon->Fill(part->y());
        costhemu->Fill(-part->w());
	if (part->ke()>0){
	  xymuons->Fill(part->x(), part->y());
	  nMuon++;
	  nMuonp++;}
}

   
  if (part->id() == CRYParticle::Electron) 
    {
        ibin = keElectron->FindBin(part->ke());
	wbin = 1.0/keElectron->GetBinWidth(ibin);     
	keElectron->Fill(part->ke(),wbin);
	if (part->charge() > 0)  kePositron->Fill(part->ke(),wbin);
	if (part->charge() < 0)  keElectronN->Fill(part->ke(),wbin);
       
}


 if (part->id() == CRYParticle::Pion) {
	ibin = kePion->FindBin(part->ke());
	wbin = 1.0/kePion->GetBinWidth(ibin);     
	kePion->Fill(part->ke(),wbin);
	if (part->charge() >  0)  kePionP->Fill(part->ke(),wbin);
	if (part->charge() == 0)  kePion0->Fill(part->ke(),wbin);
	if (part->charge() <  0)  kePionN->Fill(part->ke(),wbin);
 
}


 if (part->id() == CRYParticle::Kaon) {
	ibin = keKaon->FindBin(part->ke());
	wbin = 1.0/keKaon->GetBinWidth(ibin);     
	keKaon->Fill(part->ke(),wbin);
	if (part->charge() >  0)  keKaonP->Fill(part->ke(),wbin);
	if (part->charge() == 0)  keKaon0->Fill(part->ke(),wbin);
	if (part->charge() <  0)  keKaonN->Fill(part->ke(),wbin);
      }

  if (part->id() == CRYParticle::Neutron) {
	ibin = keNeutron->FindBin(part->ke());
	wbin = 1.0/keNeutron->GetBinWidth(ibin);     
	keNeutron->Fill(part->ke(),wbin);
      }


if (part->id() == CRYParticle::Muon) {

        ibin = keMuon->FindBin(part->ke());
	wbin = 1.0/keMuon->GetBinWidth(ibin);
        PMuon->Fill(sqrt(pow(part->ke(),2)-pow(105.65,2)),wbin);
        if (part->charge() > 0)  PMuonp->Fill(sqrt(pow(part->ke(),2)-pow(105.65,2)),wbin);
	if (part->charge() < 0)  PMuonn->Fill(sqrt(pow(part->ke(),2)-pow(105.65,2)),wbin);
//eflux vs costheta
       histo->Fill(-part->w(),wbin);
// ratio->Fill(sqrt(pow(part->ke(),2)-pow(105.65,2)),keMuonp->GetBinContent(nbins)/keMuonn->GetBinContent(nbins));
 }
  //condition à chaq angle 180-the(cos <0) voir ke
 if (part->id() == CRYParticle::Muon) {


        ibin = keMuon->FindBin(part->ke());
	wbin = 1.0/keMuon->GetBinWidth(ibin);     
	//keMuon->Fill(part->ke(),wbin);
         
     if (part->w()>-1 && part->w()<-0.98) {
        keMuon5->Fill(part->ke(),wbin);
     
      }

     if (part->w()>-0.98 && part->w()<-0.93) {
        keMuon15->Fill(part->ke(),wbin);
      }

    if (part->w()>-0.93 && part->w()<-0.86) {
        keMuon25->Fill(part->ke(),wbin);
      }

    if (part->w()>-0.86 && part->w()<-0.76) {
        keMuon35->Fill(part->ke(),wbin);
      }

   if (part->w()>-0.76 && part->w()<-0.64) {
        keMuon45->Fill(part->ke(),wbin);
      }

  if (part->w()>-0.64 && part->w()<-0.5) {
        keMuon55->Fill(part->ke(),wbin);
      }

  if (part->w()>-0.5 && part->w()<-0.3) {
        keMuon65->Fill(part->ke(),wbin);
      }

  if (part->w()>-0.3 && part->w()<-0.17) {
        keMuon75->Fill(part->ke(),wbin);
      }

  if (part->w()>-0.17 && part->w()<0) {
        keMuon85->Fill(part->ke(),wbin);
      }
}



 //condition sur Primary energy


  if (part->id() == CRYParticle::Muon) {
	ibin = keMuon->FindBin(part->ke());
	wbin = 1.0/keMuon->GetBinWidth(ibin);     
	
   if (gen.primaryParticle()->ke()>1e2 && gen.primaryParticle()->ke()<1e3){
       keMuonPrm1->Fill(part->ke(),wbin);
        if (part->charge() > 0)  keMuonPrm1p->Fill(part->ke(),wbin);
	if (part->charge() < 0)  keMuonPrm1n->Fill(part->ke(),wbin);
        costhemuPrm1->Fill(-part->w()); 
       histoPrm1->Fill(-part->w(),wbin); 
       }

  if (gen.primaryParticle()->ke()>1e3 && gen.primaryParticle()->ke()<1e4){
       keMuonPrm2->Fill(part->ke(),wbin);
        if (part->charge() > 0)  keMuonPrm2p->Fill(part->ke(),wbin);
	if (part->charge() < 0)  keMuonPrm2n->Fill(part->ke(),wbin);
       costhemuPrm2->Fill(-part->w()); 
       histoPrm2->Fill(-part->w(),wbin); 
     }

  if (gen.primaryParticle()->ke()>1e4 && gen.primaryParticle()->ke()<1e5){
       keMuonPrm3->Fill(part->ke(),wbin);
      if (part->charge() > 0)  keMuonPrm3p->Fill(part->ke(),wbin);
      if (part->charge() < 0)  keMuonPrm3n->Fill(part->ke(),wbin);
      costhemuPrm3->Fill(-part->w()); 
      histoPrm3->Fill(-part->w(),wbin); 
     }

   if (gen.primaryParticle()->ke()>1e5 && gen.primaryParticle()->ke()<1e6){
       keMuonPrm4->Fill(part->ke(),wbin);
      if (part->charge() > 0)  keMuonPrm4p->Fill(part->ke(),wbin);
      if (part->charge() < 0)  keMuonPrm4n->Fill(part->ke(),wbin);
      costhemuPrm4->Fill(-part->w()); 
     histoPrm4->Fill(-part->w(),wbin); 
    }

   if (gen.primaryParticle()->ke()>1e6 && gen.primaryParticle()->ke()<1e7){
       keMuonPrm5->Fill(part->ke(),wbin);
        if (part->charge() > 0)  keMuonPrm5p->Fill(part->ke(),wbin);
	if (part->charge() < 0)  keMuonPrm5n->Fill(part->ke(),wbin);
       costhemuPrm5->Fill(-part->w()); 
      histoPrm5->Fill(-part->w(),wbin); 
      }

}
     /* if (part->id() == CRYParticle::Muon) {
        xmuon->Fill(part->x());
        ymuon->Fill(part->y());      
        xymuons->Fill(part->x(), part->y());
        nMuon++;
	nMuonp++;
      }*/
      
    if (i % 10000 == 0){
      timeFlux = gen.timeSimulated() - timeFlux0;
      std::cout << "Event: " << i << "/" << nEv << "  " << 1.0*i/nEv*100 << " %" <<  "  " << nMuonp/(detectorsize*detectorsize*timeFlux) <<  std::endl;
      //muFlux->Fill(nMuonp/(0.16*0.16*timeFlux));
      nMuonp=0;
      timeFlux0= gen.timeSimulated();

     costhe->Fill(part->w()); 

   





     /* if (part->id() == CRYParticle::Muon) {
        xmuon->Fill(part->x());
        ymuon->Fill(part->y());      
        xymuons->Fill(part->x(), part->y());
        nMuon++;
     }*/
  
   
 }   
      delete (*ev)[j];
    }
  }

 

  chargeHist->Divide(multHist);
  std::cout << "Run completed.\n";
  std::cout << "Total time simulated: " << gen.timeSimulated() << " seconds\n";
  double muonsPerSecondPerm2=nMuon/(0.16*0.16*gen.timeSimulated());
  std::cout << "Muons per second per m2 " << muonsPerSecondPerm2 << std::endl;
//std::cout<<" ratio"<< keMuonp->GetBinContent(30)/keMuonn->GetBinContent(30)<< std::endl;
//normalisation
  double anglenorm = 1.0;
  double timenorm =  1.0/gen.timeSimulated();
  double allnorm = anglenorm*timenorm/(detectorsize*detectorsize);
  keMuon->Scale(allnorm);
  keNeutron->Scale(allnorm);
  kePion->Scale(allnorm);
  kePionP->Scale(allnorm);
  kePionN->Scale(allnorm);
  kePion0->Scale(allnorm);
  keKaon->Scale(allnorm);
  keGamma->Scale(allnorm);
  keElectron->Scale(allnorm);
  keProton->Scale(allnorm);
  kePositron->Scale(allnorm);
  keMuonn->Scale(allnorm);
  keMuonp->Scale(allnorm);
  keAntiProton->Scale(allnorm);
  keProtonp->Scale(allnorm);
  keElectronN->Scale(allnorm);
  keKaonP->Scale(allnorm);
  keKaon0->Scale(allnorm);
  keKaonN->Scale(allnorm);
  PMuon->Scale(allnorm);
 PMuonp->Scale(allnorm);
 PMuonn->Scale(allnorm);
// kePrimary->Scale(allnorm);
//µ
keMuon5->Scale(allnorm);
keMuon15->Scale(allnorm);
keMuon25->Scale(allnorm);
keMuon35->Scale(allnorm);
keMuon45->Scale(allnorm);
keMuon55->Scale(allnorm);
keMuon65->Scale(allnorm);
keMuon75->Scale(allnorm);
keMuon85->Scale(allnorm);


//KE
keMuonPrm1->Scale(allnorm);
keMuonPrm2->Scale(allnorm);
keMuonPrm3->Scale(allnorm);
keMuonPrm4->Scale(allnorm);
keMuonPrm5->Scale(allnorm);

keMuonPrm1p->Scale(allnorm);
keMuonPrm2p->Scale(allnorm);
keMuonPrm3p->Scale(allnorm);
keMuonPrm4p->Scale(allnorm);
keMuonPrm5p->Scale(allnorm);

keMuonPrm1n->Scale(allnorm);
keMuonPrm2n->Scale(allnorm);
keMuonPrm3n->Scale(allnorm);
keMuonPrm4n->Scale(allnorm);
keMuonPrm5n->Scale(allnorm);



histo->Scale(allnorm);
histoPrm1->Scale(allnorm);
histoPrm2->Scale(allnorm);
histoPrm3->Scale(allnorm);
histoPrm4->Scale(allnorm);
histoPrm5->Scale(allnorm);
//std::cout<<" ratio"<< keMuonp->GetBinContent(30)/keMuonn->GetBinContent(30)<< std::endl;

//printf(" %8.0f", keMuonp->GetBinContent(30)/keMuonn->GetBinContent(30));

//nom des axes
 keMuon->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
 keMuon->SetTitle("generated muon spectrum");


  keNeutron->GetXaxis()->SetTitle("KE_of_Neutron [Mev]");
  keNeutron->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keNeutron->SetTitle("generated neutron spectrum");

  keProton->GetXaxis()->SetTitle("KE_of_proton [Mev]");
  keProton->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keProton->SetTitle("generated proton spectrum");


  keProtonp->GetXaxis()->SetTitle("KE_of_proton [Mev]");
  keProtonp->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keProtonp->SetTitle("generated proton spectrum");
 
  keElectron->GetXaxis()->SetTitle("KE_of_electron [Mev]");
  keElectron->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keElectron->SetTitle("generated electron spectrum");


  keElectronN->GetXaxis()->SetTitle("KE_of_electron [Mev]");
  keElectronN->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keElectronN->SetTitle("generated electron spectrum");


  kePion->GetXaxis()->SetTitle("KE_of_pion[Mev]");
  kePion->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  kePion->SetTitle("generated pion spectrum");

  kePionP->GetXaxis()->SetTitle("KE_of_pion[Mev]");
  kePionP->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  kePionP->SetTitle("generated pion+ spectrum");

  kePionN->GetXaxis()->SetTitle("KE_of_pion[Mev]");
  kePionN->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  kePionN->SetTitle("generated pion- spectrum");

  kePion0->GetXaxis()->SetTitle("KE_of_pion[Mev]");
  kePion0->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  kePion0->SetTitle("generated pion0 spectrum");




 
  kePositron->GetXaxis()->SetTitle("KE_of_positron[Mev]");
  kePositron->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  kePositron->SetTitle("generated positron spectrum");


 keMuonn->GetXaxis()->SetTitle("KE_of_Muonn [Mev]");
  keMuonn->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuonn->SetTitle("generated muonn spectrum");



  keMuonp->GetXaxis()->SetTitle("KE_of_Muonp [Mev]");
  keMuonp->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuonp->SetTitle("generated muonp spectrum");

  keAntiProton->GetXaxis()->SetTitle("KE_of_AntiProton [Mev]");
  keAntiProton->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keAntiProton->SetTitle("generated AntiProton spectrum");


  keKaon->GetXaxis()->SetTitle("KE_of_kaon [Mev]");
  keKaon->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keKaon->SetTitle("generated kaon spectrum");


  keKaonP->GetXaxis()->SetTitle("KE_of_kaon [Mev]");
  keKaonP->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keKaonP->SetTitle("generated kaon+ spectrum");

  keKaonN->GetXaxis()->SetTitle("KE_of_kaon [Mev]");
  keKaonN->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keKaonN->SetTitle("generated kaon- spectrum");

  keKaon0->GetXaxis()->SetTitle("KE_of_kaon [Mev]");
  keKaon0->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keKaon0->SetTitle("generated kaon0 spectrum");






ratio->GetXaxis()->SetTitle("momentum_of_Muon [Mev/c]");
  ratio->GetYaxis()->SetTitle("ratio");
  ratio->SetTitle("ratio");

  kePrimary->GetXaxis()->SetTitle("KE_of_Primary [Mev]");
  //kePrimary->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  kePrimary->SetTitle(" ke of primary");


 keMuon5->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon5->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon5->SetTitle("generated muon at 5° spectrum");


  keMuon15->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon15->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon15->SetTitle("generated muon at 15° spectrum");


  keMuon25->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon25->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon25->SetTitle("generated muon at 25° spectrum");


  keMuon35->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon35->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon35->SetTitle("generated muon at 35° spectrum");


  keMuon45->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon45->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon45->SetTitle("generated muon at 45° spectrum");


  keMuon55->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon55->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon55->SetTitle("generated muon at 55° spectrum");

 
  keMuon65->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon65->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon65->SetTitle("generated muon at 65° spectrum");


  keMuon75->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon75->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon75->SetTitle("generated muon at 75° spectrum");


  keMuon85->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuon85->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuon85->SetTitle("generated muon at 85° spectrum");


  keMuonPrm1->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuonPrm1->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuonPrm1->SetTitle("ke of muon for 10^2<keprimary<10^3 Mev ");


  keMuonPrm2->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuonPrm2->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuonPrm2->SetTitle("ke of muon for 10^3<keprimary<10^4 Mev ");

  keMuonPrm3->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuonPrm3->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuonPrm3->SetTitle("ke of muon for 10^4<keprimary<10^5 Mev ");

  keMuonPrm4->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuonPrm4->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuonPrm4->SetTitle("ke of muon for 10^5<keprimary<10^6 Mev ");

  keMuonPrm5->GetXaxis()->SetTitle("KE_of_Muon [Mev]");
  keMuonPrm5->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  keMuonPrm5->SetTitle("ke of muon for 10^6<keprimary<10^7 Mev ");


  histo->GetXaxis()->SetTitle("costheta");
  histo->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");
  
  histoPrm1->GetXaxis()->SetTitle("costheta");
  histoPrm1->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");

  histoPrm2->GetXaxis()->SetTitle("costheta");
  histoPrm2->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");

  histoPrm3->GetXaxis()->SetTitle("costheta");
  histoPrm3->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");

  histoPrm4->GetXaxis()->SetTitle("costheta");
  histoPrm4->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");

  histoPrm5->GetXaxis()->SetTitle("costheta");
  histoPrm5->GetYaxis()->SetTitle("Vertical intensity[1/m^2/s/sr/Mev]");

  // Write the histogram file
  outputFile->Write();
  outputFile->Close();

  delete setup;

  return 0;
}
