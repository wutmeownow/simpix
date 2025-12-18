// simple example of using ROOT libraries in a C++ program with graphics
// and use of TASImage class

#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"
#include "TSystem.h"
#include <TStopwatch.h>
#include <TRandom3.h>
#include "assert.h"
#include <iostream>
#include <stdio.h>
#include <thread>

using namespace std;

// calculate apparent difference between two pixels
double PixDiff(const UInt_t p1, const UInt_t p2) {
  // R byte is byte 2, G is byte 1, B is byte 0
  int r1 = (p1 >> 16) & 0xFF;
  int g1 = (p1 >> 8) & 0xFF;
  int b1 = p1 & 0xFF;

  int r2 = (p2 >> 16) & 0xFF;
  int g2 = (p2 >> 8) & 0xFF;
  int b2 = p2 & 0xFF;

  // calculate differences
  int dr = r2 - r1;
  int dg = g2 - g1;
  int db = b2 - b1;

  // We scale the standard formula by 256 to eliminate floating point division.
  // Original: (2 + ravg/256) * dr^2 + 4 * dg^2 + (2 + (255-ravg)/256) * db^2
  // Scaled:   (512 + ravg) * dr^2   + 1024 * dg^2 + (767 - ravg) * db^2
  
  int ravg = (r1 + r2) >> 1; // Bitwise shift >> 1 is the same as / 2
  
  long long diff = (512 + ravg) * dr * dr
                  + 1024 * dg * dg
                  + (767 - ravg) * db * db;

  return (double)diff;
}

// calculate total difference in pixel colors between two pixel arrays
double CalcTotalDiff(int npix, const UInt_t *parr1, const UInt_t *parr2) {
  double total_diff = 0.;
  for (int i=0; i<npix; i++) {
    total_diff += PixDiff(parr1[i], parr2[i]);
  }
  return total_diff;
}

// carry out pixel swap according to metropolis algorithm
double Update(int npix, TRandom3 *r, UInt_t *psrc, const UInt_t *ptgt, double T=-1.) {
  int iRand = 0; // random pixel 1
  int jRand = 0; // random pixel 2
  // swapping the same pixel changes nothing, keep generating until they are different
  while (iRand == jRand) {
    iRand = r->Integer(npix-1); 
    jRand = r->Integer(npix-1); 
  }

  // change in energy is the difference between the pixel differences now and after swap
  double dold = PixDiff(psrc[iRand],ptgt[iRand]) + PixDiff(psrc[jRand],ptgt[jRand]);
  double dnew = PixDiff(psrc[jRand],ptgt[iRand]) + PixDiff(psrc[iRand],ptgt[jRand]);
  double dE = dnew-dold;


  if (T<0) {
    // hot, always take the swap
    UInt_t pi = psrc[iRand];
    psrc[iRand] = psrc[jRand];
    psrc[jRand] = pi;
    return dE;
  } 
  else {
    double p = r->Rndm();
    if (dE<0 || p < exp(-dE/T)) {
      UInt_t pi = psrc[iRand];
      psrc[iRand] = psrc[jRand];
      psrc[jRand] = pi;
      return 1.; // successful change
    }
    return 0.; // no change
  }
}

double AnnealingTMin(int npix, UInt_t *outPix, const UInt_t *tgtPix, int N, double Tf, double alpha, TRandom3 *r){
  // heat up and get Tmax
  double Tmax = 0;
  for (int i=0; i<npix*N;i++) {
    double dE = Update(npix, r, outPix, tgtPix);
    if (dE>Tmax) {Tmax=dE;}
  }
  
  // now do annealing, stopping only when the decrease in path is negligible OR we reach Tf
  double T = Tmax;

  // stop annealing when Tf is reached
  while (T>Tf) {
    // keep applying changes to the order until we hit the ncity*N changes at this T
    double nchanges = 0;
    while (nchanges<npix*N*1.) {
      nchanges += Update(npix, r, outPix, tgtPix, T);
    }

    // calculate the final energy for this T
    // double currE = CalcTotalDiff(numPix, outPix, tgtPix); 

    // store annealing history
    // history_T.push_back(T);
    // history_Dist.push_back(currL);

    // update T
    T = alpha*T;
  }

  return T;
}


double AnnealingLimit(int npix, UInt_t *outPix, const UInt_t *tgtPix, int N, double init_diff, double limit, double alpha, TRandom3 *r){
  // heat up and get Tmax
  double Tmax = 0;
  for (int i=0; i<npix*N;i++) {
    double dE = Update(npix, r, outPix, tgtPix);
    if (dE>Tmax) {Tmax=dE;}
  }

  // stop annealing when change is less than limit
  double T = Tmax;
  double prev_E = init_diff;
  double E_change = 1.0;
  while (E_change>limit) {
    // keep applying changes to the order until we hit the npix*N changes at this T
    double nchanges = 0;
    while (nchanges<npix*N*1.) {
      nchanges += Update(npix, r, outPix, tgtPix, T);
    }

    // calculate the final energy for this T
    double currE = CalcTotalDiff(npix, outPix, tgtPix); 
    E_change = abs(currE - prev_E)/prev_E;
    prev_E = currE;

    // store annealing history
    // history_T.push_back(T);
    // history_Dist.push_back(currL);

    // update T
    T = alpha*T;
  }
  return T;
}


int main(int argc, char **argv){
  TStopwatch timer; // timer to track how long code takes
  timer.Start();

  // We need two random generators so the threads don't fight over the same one
  TRandom3 *r1 = new TRandom3(0);
  TRandom3 *r2 = new TRandom3(0);

  thread t1, t2; // two thread handles

  // parameters
  TString fsrc="";
  TString ftgt="";
  TString fout="";
  double alpha = 0.8; // next T parameter
  int N = 10; // multiplicative factor for number of trials
  double limit = -1.; // limit factor for when to stop annealing
  double Tf = -1.; // final temp for when to stop annealing
  int opt;

  // Loop through all arguments
  // "n:f:" means we expect flags -n and -f, and both require values (:)
  while ((opt = getopt(argc, argv, "s:f:o:n:a:L:t:")) != -1) {
      switch (opt) {
          case 's':
            fsrc = optarg;
            break;
          case 'f':
            ftgt = optarg;
            break;
          case 'o':
            fout = optarg;
            break;
          case 'a':
            alpha = std::stof(optarg);
            break;
          case 'n':
            N = std::atoi(optarg);
            break;
          case 'L':
            limit = std::stof(optarg);
            break;
          case 't':
            Tf = std::stof(optarg);
            break;
          default:
            std::cerr << "Usage: " << argv[0] << " -s <filename1> -f <filename2> -o <outputfile> -a <alpha> -n <N> -t <Tmin> -L <limit>" << std::endl;
            return 1;
      }
  }

  if (fsrc==""||ftgt==""||fout=="") {
      std::cerr << "Error: Filenames (-s -f -o) are required!" << std::endl;
      return 1;
  }

  if (Tf == -1. && limit == -1.) {
    std::cerr << "Error: limit (-L) or Tmin (-t) required!" << std::endl;
    return 1;
  }

  if (Tf != -1. && limit != -1.) {
    std::cerr << "Error: please only provide either limit (-L) or Tmin (-t) to choose when annealing stops!" << std::endl;
    return 1;
  }


  if (argc<3) {
    cout << "Usage: simapix_start image1 image2 <output=out.png>" << endl;
    return 0; 
  }

  cout << "Reading images: source= " << fsrc << " target= " << ftgt << endl;
  cout << "Output= " << fout << endl;


  // create image objects
  TASImage *src = new TASImage(fsrc.Data());
  TASImage *tgt = new TASImage(ftgt.Data());
  TASImage *out = new TASImage(*src); // start with copy of source
  TASImage *out2 = new TASImage(*tgt); // start with copy of tgt

  // Test image geometry, exit if they are not the same dimensions
  cout << "SRC Pixel Geometry: " << src->GetWidth() << " x " << src->GetHeight() << endl;
  cout << "TGT Pixel Geometry: " << tgt->GetWidth() << " x " << tgt->GetHeight() << endl;
  assert ( src->GetWidth() == tgt->GetWidth() && src->GetHeight() == tgt->GetHeight() );
  Long_t numPix=src->GetWidth()*src->GetHeight();

  // *** The work happens here
  // access the pixels for the output image 
  // each pixel is a 32-bit word, 1 byte each for (alpha,red,green,blue)
  // don't touch alpha (bits 31:28)
  UInt_t *outPix = out->GetArgbArray();
  UInt_t *outPix2 = out2->GetArgbArray();
  const UInt_t *tgtPix = tgt->GetArgbArray();
  const UInt_t *srcPix = src->GetArgbArray();  

  // starting difference
  double diff1 = CalcTotalDiff(numPix, outPix, tgtPix);
  cout<<"Initial Diff: "<<diff1<<endl;

  // either do annealing for both until Tmin is reached or the difference is under the limit
  // use multithreading to speed things up 
  if (Tf!=-1.) {
    t1 = thread(AnnealingTMin, numPix, outPix, tgtPix, N, Tf, alpha, r1);
    t2 = thread(AnnealingTMin, numPix, outPix2, srcPix, N, Tf, alpha, r2);
  } else {
    t1 = thread(AnnealingLimit, numPix, outPix, tgtPix, N, diff1, limit, alpha, r1);
    t2 = thread(AnnealingLimit, numPix, outPix2, srcPix, N, diff1, limit, alpha, r2);
  }

  // Wait for both to finish before continuing
  t1.join();
  t2.join();

  // stop timer after annealing is finished
  timer.Stop();

  // print the time it took
  timer.Print(); // Prints Real time (wall clock) and Cpu time

  // ending difference
  double enddiff1 = CalcTotalDiff(numPix, outPix, tgtPix);
  double enddiff2 = CalcTotalDiff(numPix, outPix2, srcPix);
  cout<<"Ending Diff 1: "<<enddiff1<<endl;
  cout<<"Ending Diff 2: "<<enddiff2<<endl;

  // *************************


  // print the results
  TCanvas *c1 = new TCanvas("c1", "images", 640, 480);
  c1->Divide(2,2);

  c1->cd(1);
  c1->Draw();
  src->Draw("X");
  c1->cd(2);
  tgt->Draw("X");
  c1->cd(3);
  out->Draw("X");
  c1->cd(4);
  out2->Draw("X");
  c1->Print(fout); 
  // save as png first, then convert to pdf externally later

  // Clean up
  delete r1;
  delete r2;

  return 1;
}
