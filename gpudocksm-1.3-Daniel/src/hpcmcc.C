//Daniel Case 11-16-12
//This program calculates the cross MCCs for a set of ligand coords


#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include "hpcmcc.h"

using namespace std;

//std::fstream& GotoLine(std::fstream&, unsigned int);

int main(int argc, char *argv[]){

string molid         = "MOLID";


 if ( argc < 4 )
 {
  cout << " gpudocksm -p  <target protein structure, PDB>" << endl
       << "           -l  <compound library in SD format>" << endl
      // << "           -s  <LHM potentials>" << endl
       << "           -o  <docked compounds in SD format>" << endl << endl;
/*
       << " additional options:" << endl
       << "           -i  <molecule id keyword (default MOLID)>" << endl << endl

       << " remc options:" << endl
       << "           -nr <number of replicas (default 6)>" << endl
       << "           -ns <number of iterations for each replica (default 50)>" << endl
       << "           -nc <number of replica swaps (default 50)>" << endl << endl;
*/
  exit(EXIT_SUCCESS);
 }

 string protein_name;
 bool protein_opt = false;

 string compounds_name;
 bool compounds_opt = false;

 string output_name;
 bool output_opt = false;

 string data_name;

 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-p")  && i < argc ) { protein_name   = string(argv[i+1]); protein_opt   = true; }
  if ( !strcmp(argv[i],"-l")  && i < argc ) { compounds_name = string(argv[i+1]); compounds_opt = true; }
//  if ( !strcmp(argv[i],"-s")  && i < argc ) { lhm_name       = string(argv[i+1]); lhm_opt       = true; }
  if ( !strcmp(argv[i],"-o")  && i < argc ) { output_name    = string(argv[i+1]); output_opt    = true; }
//  if ( !strcmp(argv[i],"-i")  && i < argc ) { molid          = string(argv[i+1]);                       }
//  if ( !strcmp(argv[i],"-nr") && i < argc ) { remc_replicas  = atoi(argv[i+1]);                         }
//  if ( !strcmp(argv[i],"-ns") && i < argc ) { remc_steps     = atoi(argv[i+1]);                         }
//  if ( !strcmp(argv[i],"-nc") && i < argc ) { remc_cycles    = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-d")  && i < argc ) { data_name      = string(argv[i+1]);}
 }
 
 char * path1;

 path1 = getenv("GPUDOCKSMDAT_FF"); if ( path1==NULL ) { cout << "GPUDOCKSMDAT_FF is not set" << endl; exit(EXIT_FAILURE); }

 string data_path;
 data_path = getenv("GPUDOCKSMDAT_FF");

 if ( !protein_opt )
 {
  cout << "Provide target protein structure" << endl;
  exit(EXIT_FAILURE);
 }

 if ( !compounds_opt )
 {
  cout << "Provide compound library in SD format" << endl;
  exit(EXIT_FAILURE);
 }
 else
 {
  ofstream outprot( output_name.c_str() );
  //cout<<"test3"<<endl;
  outprot.close();
 }

 list<string> l1_sdf;
 list<string>::iterator i1_sdf;

 string line1;

 ifstream compounds_file( compounds_name.c_str() );

 if ( !compounds_file.is_open() )  { cout << "Cannot open " << compounds_name << endl; exit(EXIT_FAILURE); }

 while (getline(compounds_file,line1)){
  l1_sdf.push_back( line1 );
//cout<<"test4"<<endl;
}
 compounds_file.close();

 int gpu_n = 0;

 string lib1[MAXSDF];
 int lib2 = 0;

 Complex * gpu_complex[MAXLIB];
//cout<<"test5"<<endl;
 for ( i1_sdf = l1_sdf.begin(); i1_sdf != l1_sdf.end(); i1_sdf++ )
 {
 //cout<<"test5.5"; 
  lib1[lib2++] = (*i1_sdf);
  //cout<<endl;
  //cout<<lib1[lib2]<<endl;
//cout<<"test6";
  if ( (*i1_sdf) == "$$$$" )
  {
   //cout<<"test7";
   if ( lib2 > 10 )
   {
    gpu_complex[gpu_n] = new Complex( 0, 0 );
    // gpu_complex[gpu_n]->clearMoves();
    //cout<<"test8";
    
    bool load1 = gpu_complex[gpu_n]->loadProtein(protein_name);
    //cout<<"test9";
    if ( load1 )
    {
     cout << "Cannot read target protein structure" << endl;
     exit(EXIT_FAILURE);
    }
    
    bool load2 = false;
    load2 = gpu_complex[gpu_n]->loadParams(data_path);
    //cout<<"test10";
    if ( load2 )
    {
     cout << "Cannot read parameter file" << endl;
     exit(EXIT_FAILURE);
    }
    
    bool load3 = gpu_complex[gpu_n]->loadLigand(lib1, lib2, molid);
    //cout<<"test11";
    if ( load3 )
    {
     cout << "Cannot load ligand data" << endl;
     exit(EXIT_FAILURE);
    }
   
    // cout<<"test12";
     gpu_n++;
   }
   
   lib2 =0;
   
  }
// cout<<"test5.75";
 }
//cout<<"test100"<<endl;
for(int in=0; in<gpu_n; in++){
// cout<<"THIS IS ANOTHER TEST"<<endl;

 gpu_complex[in]->createContactMatrix();
//cout<<"THIS IS ANOTHER TEST"<<endl;
 Decoy(gpu_complex[in], data_name.c_str(), output_name.c_str());

 }
 return 0;
}
/*
std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}
*/
