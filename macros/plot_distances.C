#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <istream>

void plot_distances(){

    TCanvas *c = new TCanvas();
    TH1D    *h = new TH1D( "h", "Distances between true and reconstructed 'primary' vertex", 200, 0, 200 );

    std::string line;
    std::vector< double > dists;
    std::stringstream stringtodouble; 


    ifstream myfile( "/sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/distances.txt" );

    if( myfile.is_open() ){

        while( std::getline( myfile, line ) ){

            stringtodouble.clear();
            stringtodouble.str( std::string() );
            double temp = 0.0;

            stringtodouble << line;
            stringtodouble >> temp;

            dists.push_back( temp );

        }
    
        myfile.close();
    }

    for( unsigned int i = 0; i < dists.size(); ++i ){
        h->Fill( dists[i] );
    }

    h->Draw();
    c->Update();

    TFile *f = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/hist_distances.root", "RECREATE" );
    h->Write();
    f->Close();

    c->SaveAs("/sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/hist_dist.pdf" );
    c->SaveAs("/sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/hist_dist.png" );

    delete h;
    delete c;
}
