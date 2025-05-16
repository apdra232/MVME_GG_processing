#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include "TStopwatch.h"
#include <cmath>

vector<Double_t> m={1.,1.,1.,1.};
vector<Double_t> b={0.,0.,0.,0.};

void Combine_Spectrum_Smoothly(TH1* spect, TH1 *additional_spectrum)
{
  int nbins = spect->GetNbinsX();
  for(int i=0; i<=nbins; i++)
  {
    double low_edge = additional_spectrum->GetBinLowEdge(i);
    double high_edge = low_edge + additional_spectrum->GetBinWidth(i);
    int start_bin = spect->FindFixBin(low_edge);
    int stop_bin = spect->FindFixBin(high_edge);
    int num_overlap = 1+ stop_bin-start_bin;
    if(num_overlap==1)
      spect->AddBinContent(start_bin,additional_spectrum->GetBinContent(i));
    else
    {
      double overlap = (spect->GetBinLowEdge(start_bin) + spect->GetBinWidth(start_bin))-low_edge;
      spect->AddBinContent(start_bin,additional_spectrum->GetBinContent(i)*overlap/(high_edge-low_edge));
      start_bin++;
      while(start_bin<stop_bin)
      {
        overlap = spect->GetBinWidth(start_bin);
        spect->AddBinContent(start_bin,additional_spectrum->GetBinContent(i)*overlap/(high_edge-low_edge));
        start_bin++;
      }
      overlap = high_edge - spect->GetBinLowEdge(stop_bin);
      spect->AddBinContent(stop_bin,additional_spectrum->GetBinContent(i)*overlap/(high_edge-low_edge));
    }
  }
}

void calib_dat(vector<Double_t>& m, vector<Double_t>& b, TString calibfile){
    string dummy, line;
    int k=0;
    ifstream data(calibfile.Data());
    while(getline(data,line)){
            istringstream s1(line);
            s1 >> dummy >> dummy >> b[k];
            getline(data,line);
            istringstream s2(line);
            s2 >> dummy >> dummy >> m[k];
            k++;
    }
    cout << calibfile << " found!" << endl;
    cout << "Using following calibration parameters!" << endl;
    for(int k=0; k<4; k++)
            cout << "m: " << m[k] << " b: " << b[k] << endl;
    data.close();
}

void check_calib(TString fileName){
	ifstream data(fileName.Data());
	if(!data){
		cerr << fileName << " not found!" << endl;
		cerr << "Histograms not calibrated!" << endl;
	}
	else
		calib_dat(m,b,fileName);
}

void process_mvmeTree_3ch(TString filename) {
    // Open the ROOT file
    TStopwatch totaltimer;
    totaltimer.Start();

    TString calibfile;
    int hbin = 16384;
    int fbin = 2;
    // check calibration and if it exists, implement
    calibfile = "calib.dat";
    check_calib(calibfile);
    
    // read input file to add to TChain
    string fname;
    ifstream evtfile(filename.Data());
    evtfile >> fname;
    cout << fname << endl;
    TFile *f = new TFile(fname.c_str(),"read");
    const string getstring = f->GetListOfKeys()->At(0)->GetName();
    // cout << getstring << endl;

    TChain *t1;
    t1 = new TChain(getstring.c_str());
    t1->Add(fname.c_str());

    // add ttree to tchain -> t1
    while(evtfile>>fname)
    {
            t1->Add(fname.c_str());
            cout << "Add file: " << fname << endl;
    }
    evtfile.close();

    // Define variables to hold branch data
    double mdpp16_scp_channel_time[7];
    double mdpp16_scp_module_timestamp; // Not used in subtraction, but can be read if needed
    double mdpp16_scp_amplitude[7];
    double mdpp16_scp_pileup[7];
    // Set branch addresses
    t1->SetBranchAddress("mdpp16_scp.channel_time", mdpp16_scp_channel_time);
    t1->SetBranchAddress("mdpp16_scp.module_timestamp", &mdpp16_scp_module_timestamp);
    t1->SetBranchAddress("mdpp16_scp.amplitude", mdpp16_scp_amplitude);
    t1->SetBranchAddress("mdpp16_scp.pileup", mdpp16_scp_pileup);

    // Create a histogram to store subtracted values
    TH1D *dt[4];
    // TH1D *hist = new TH1D("subtracted_channel_time", "Subtracted Channel Time;Subtracted Value;Counts", 16384, 0, 65530);
    TH1D *ghist[4];
    TH1D *CShist[4];
    TH1D *bhist[4];
    TH1D *chist[4];
    TH1D *hall;
    TH2D *h2d[4][4];
    TH2D *bh2d[4][4]; 
    TH1D *tstamp;
    dt[0] = new TH1D(Form("DT00"), Form("DT00"), 1000, -30000, 30000);
    dt[1] = new TH1D(Form("DT12"), Form("DT12"), 1000, -30000, 30000);
    dt[2] = new TH1D(Form("DT13"), Form("DT13"), 1000, -30000, 30000);
    dt[3] = new TH1D(Form("DT23"), Form("DT23"), 1000, -30000, 30000);
    hall = new TH1D(Form("coinc_all"), Form("coinc_all"), hbin,0.*m[1]+b[1],(hbin)*m[1]+b[1]);
    tstamp= new TH1D("Timestamp", "Timestamp", hbin, 0, 86400);
    for(int i=0; i<4; i++)
    {
        ghist[i]= new TH1D(Form("AllHist%d",i), Form("AllHist%d",i), hbin,0.*m[i]+b[i],(hbin)*m[i]+b[i]);
        CShist[i]= new TH1D(Form("ComptonSuppHist%d",i), Form("ComptonSuppHist%d",i), hbin,0.*m[i]+b[i],(hbin)*m[i]+b[i]);
        bhist[i]= new TH1D(Form("BackgroundHist%d",i), Form("BackgroundHist%d",i), hbin,0.*m[i]+b[i],(hbin)*m[i]+b[i]);
        chist[i]= new TH1D(Form("CoincHist%d",i), Form("CoincHist%d",i), hbin,0.*m[i]+b[i],(hbin)*m[i]+b[i]);
        for(int j=0; j<4; j++)
        {
            h2d[i][j] = new TH2D(Form("gg%i%i",i,j),Form("dt%i%i",i,j),hbin/4.,0.*m[i]+b[i],(hbin)*m[i]+b[i],hbin/fbin,0.*m[j]+b[j],(hbin+1.)*m[j]+b[j]);
            bh2d[i][j] = new TH2D(Form("bgg%i%i",i,j),Form("dt%i%i",i,j),hbin/4.,0.*m[i]+b[i],(hbin)*m[i]+b[i],hbin/fbin,0.*m[j]+b[j],(hbin+1.)*m[j]+b[j]);
        }
    }

    // Loop over all entries in the TTree
    Long64_t nEntries = t1->GetEntries();
    cout << "Total Entries: " << nEntries << endl;
    double mdpp16_scp_channel_time_old=0;
    double energy1,energy2;
    double tstamp_old=0.;
    double tstamp_prev=0.;
    for (Long64_t i = 0; i < nEntries; i++) {
        t1->GetEntry(i);
        if(i%(ULong_t)((double)(nEntries-1)*0.1)==0){
            printf("... processing %3.0f%% of data (%d s)... \n", (double)i/(double)nEntries*100, (int)totaltimer.RealTime());
            totaltimer.Continue();
        }

        if(mdpp16_scp_module_timestamp>tstamp_prev){
            tstamp->Fill((mdpp16_scp_module_timestamp+tstamp_old)/16.*1.e-6);
        }
        else{
            tstamp_old+=tstamp_prev;
            tstamp_prev=0;
            cout << (mdpp16_scp_module_timestamp+tstamp_old)/16.*1.e-6 << endl;
            tstamp->Fill((mdpp16_scp_module_timestamp+tstamp_old)/16.*1.e-6);
        }

        if(mdpp16_scp_amplitude[1]>0 && mdpp16_scp_pileup[1]==0){
            ghist[1]->Fill(mdpp16_scp_amplitude[1]/4.*m[1]+b[1]);
            if(mdpp16_scp_amplitude[4]<50)
                CShist[1]->Fill(mdpp16_scp_amplitude[1]/4.*m[1]+b[1]);
        }      
        if(mdpp16_scp_amplitude[2]>0 && mdpp16_scp_pileup[2]==0){ 
            ghist[2]->Fill(mdpp16_scp_amplitude[2]/4.*m[2]+b[2]);
            if(mdpp16_scp_amplitude[5]<50)
                CShist[2]->Fill(mdpp16_scp_amplitude[2]/4.*m[2]+b[2]);
        }      
        if(mdpp16_scp_amplitude[3]>0 && mdpp16_scp_pileup[3]==0){
            ghist[3]->Fill(mdpp16_scp_amplitude[3]/4.*m[3]+b[3]);
            if(mdpp16_scp_amplitude[6]<50)
                CShist[3]->Fill(mdpp16_scp_amplitude[3]/4.*m[3]+b[3]);
        }      
        if(mdpp16_scp_amplitude[1]>0 && mdpp16_scp_amplitude[2]>0 && mdpp16_scp_pileup[1]==0 && mdpp16_scp_pileup[2]==0){
            double subtraction = mdpp16_scp_channel_time[1] - mdpp16_scp_channel_time[2];
            energy1 = mdpp16_scp_amplitude[1]/4.*m[1]+b[1];
            energy2 = mdpp16_scp_amplitude[2]/4.*m[2]+b[2];
            dt[1]->Fill(subtraction);
            dt[1]->Fill(-1*subtraction);
            if(abs(subtraction)<5000){
                h2d[1][2]->Fill(energy1,energy2);
                h2d[2][1]->Fill(energy2,energy1);
                hall->Fill(energy1);
                hall->Fill(energy2);
                chist[1]->Fill(energy1);
                chist[1]->Fill(energy2);
            }
            else if((subtraction > -12000 && subtraction < -7000) || (subtraction > 7000 && subtraction < 12000 )){
                bh2d[1][2]->Fill(energy1,energy2);
                bh2d[2][1]->Fill(energy2,energy1);
                bhist[1]->Fill(energy1);
                bhist[1]->Fill(energy2);
            }
        }
        else if(mdpp16_scp_amplitude[1]>0 && mdpp16_scp_amplitude[3]>0 && mdpp16_scp_pileup[1]==0 && mdpp16_scp_pileup[3]==0){
            double subtraction = mdpp16_scp_channel_time[1] - mdpp16_scp_channel_time[3];
            energy1 = mdpp16_scp_amplitude[1]/4.*m[1]+b[1];
            energy2 = mdpp16_scp_amplitude[3]/4.*m[3]+b[3];
            dt[2]->Fill(subtraction);
            dt[2]->Fill(-1*subtraction);
            if(abs(subtraction)<5000){
                h2d[1][3]->Fill(energy1,energy2);
                h2d[3][1]->Fill(energy2,energy1);              
                hall->Fill(energy1);
                hall->Fill(energy2);  
                chist[2]->Fill(energy1);
                chist[2]->Fill(energy2);
            }
            else if((subtraction > -12000 && subtraction < -7000) || (subtraction > 7000 && subtraction < 12000)){
                bh2d[1][3]->Fill(energy1,energy2);
                bh2d[3][1]->Fill(energy2,energy1);
                bhist[2]->Fill(energy1);
                bhist[2]->Fill(energy2);
            }
        }
        else if(mdpp16_scp_amplitude[2]>0 && mdpp16_scp_amplitude[3]>0 && mdpp16_scp_pileup[2]==0 && mdpp16_scp_pileup[3]==0){
            double subtraction = mdpp16_scp_channel_time[2] - mdpp16_scp_channel_time[3];
            energy1 = mdpp16_scp_amplitude[2]/4.*m[2]+b[2];
            energy2 = mdpp16_scp_amplitude[3]/4.*m[3]+b[3];
            dt[3]->Fill(subtraction);
            dt[3]->Fill(-1*subtraction);
            if(abs(subtraction)<5000){
                h2d[2][3]->Fill(energy1,energy2);
                h2d[3][2]->Fill(energy2,energy1);             
                hall->Fill(energy1);
                hall->Fill(energy2);   
                chist[3]->Fill(energy1);
                chist[3]->Fill(energy2);
            }
            else if((subtraction > -12000 && subtraction < -7000) || (subtraction > 7000 && subtraction < 12000)){
                bh2d[2][3]->Fill(energy1,energy2);
                bh2d[3][2]->Fill(energy2,energy1);
                bhist[3]->Fill(energy1);
                bhist[3]->Fill(energy2);
            }
        }
        tstamp_prev=mdpp16_scp_module_timestamp;
    }

    TFile bam(Form("gg_%s.root",filename.Data()),"recreate");
    for(int i=0;i<4;i++){
        ghist[i]->Write();
        CShist[i]->Write();
        bhist[i]->Write();
        chist[i]->Write();
        dt[i]->Write();
        for(int j=0;j<4;j++){
            h2d[i][j]->Write();
            bh2d[i][j]->Write();
        }
    }
    hall->Write();
    tstamp->Write();
    bam.Close();

    totaltimer.Stop();
    cout << "Finished in " << totaltimer.RealTime() << " s" << endl;
}