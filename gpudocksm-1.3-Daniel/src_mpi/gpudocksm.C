/*
==============================================================================================
     __________________ ____ ___              .___             __      _________   _____   
    /  _____/\______   \    |   \           __| _/____   ____ |  | __ /   _____/  /     \  
   /   \  ___ |     ___/    |   /  ______  / __ |/  _ \_/ ___\|  |/ / \_____  \  /  \ /  \ 
   \    \_\  \|    |   |    |  /  /_____/ / /_/ (  <_> )  \___|    <  /        \/    Y    \
    \______  /|____|   |______/           \____ |\____/ \___  >__|_ \/_______  /\____|__  /
           \/                                  \/           \/     \/        \/         \/ 

      GPU-accelerated hybrid-resolution ligand docking using Replica Exchange Monte Carlo

==============================================================================================
*/


// SEE complex.C FOR NOTES FROM DANIEL//



#include "gpudocksm.h"

using namespace std;

int main(int argc, char *argv[])
{
 int    remc_steps    = 10;
 int    remc_cycles   = 10;
 int    remc_replicas = 6;
////////////////////////////////////////////////////////////
 double step_t        = 0.5;
 double step_r        = 5.0;
 double dsphere       = 6.0;
 double boltzmann     = 1.0;

////////////////////////////////////////////////////////////
 string molid         = "MOLID";
 
 srand(time(0));
 
 time_t t_start1, t_end1, t_start2, t_end2;
 
 time(&t_start1);
 
 cout << "------------------------------------------------------------" << endl
      << "                         GPU-dockSM" << endl
      << "                         version 1.0" << endl << endl
      << "   GPU-accelerated mixed-resolution ligand docking using" << endl
      << "                Replica Exchange Monte Carlo" << endl
      << "------------------------------------------------------------" << endl << endl;
 
 if ( argc < 7 )
 {
  cout << " gpudocksm -p  <target protein structure, PDB>" << endl
       << "           -l  <compound library in SD format>" << endl
       << "           -s  <LHM potentials>" << endl
       << "           -o  <docked compounds in SD format>" << endl << endl
       
       << " additional options:" << endl
       << "           -i  <molecule id keyword (default MOLID)>" << endl << endl
       
       << " remc options:" << endl
       << "           -nr <number of replicas (default 6)>" << endl
       << "           -ns <number of iterations for each replica (default 50)>" << endl
       << "           -nc <number of replica swaps (default 50)>" << endl << endl;
  
  exit(EXIT_SUCCESS);
 }
 
 string protein_name;
 bool protein_opt = false;
 
 string compounds_name;
 bool compounds_opt = false;
 
 string lhm_name;
 bool lhm_opt = false;
 
 string output_name;
 bool output_opt = false;
 
 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-p")  && i < argc ) { protein_name   = string(argv[i+1]); protein_opt   = true; }
  if ( !strcmp(argv[i],"-l")  && i < argc ) { compounds_name = string(argv[i+1]); compounds_opt = true; }
  if ( !strcmp(argv[i],"-s")  && i < argc ) { lhm_name       = string(argv[i+1]); lhm_opt       = true; }
  if ( !strcmp(argv[i],"-o")  && i < argc ) { output_name    = string(argv[i+1]); output_opt    = true; }
  if ( !strcmp(argv[i],"-i")  && i < argc ) { molid          = string(argv[i+1]);                       }
  if ( !strcmp(argv[i],"-nr") && i < argc ) { remc_replicas  = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-ns") && i < argc ) { remc_steps     = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-nc") && i < argc ) { remc_cycles    = atoi(argv[i+1]);                         }
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
 
 if ( !lhm_opt )
 {
  cout << "Provide LHM potentials" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !output_opt )
 {
  cout << "Provide output filename" << endl;
  exit(EXIT_FAILURE);
 }
 else
 {
  ofstream outprot( output_name.c_str() );
  outprot.close();
 }
 
 list<string> l1_sdf;
 list<string>::iterator i1_sdf;
 
 string line1;
 
 ifstream compounds_file( compounds_name.c_str() );
 
 if ( !compounds_file.is_open() )  { cout << "Cannot open " << compounds_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(compounds_file,line1))
  l1_sdf.push_back( line1 );
 
 compounds_file.close();
 
 int gpu_n = 0;
 
 string lib1[MAXSDF];
 int lib2 = 0;
 
 Complex * gpu_complex[MAXLIB];
 
 for ( i1_sdf = l1_sdf.begin(); i1_sdf != l1_sdf.end(); i1_sdf++ )
 {
  lib1[lib2++] = (*i1_sdf);
  
  if ( (*i1_sdf) == "$$$$" )
  {
   if ( lib2 > 10 )
   {
    gpu_complex[gpu_n] = new Complex( 0, 0 );
     gpu_complex[gpu_n]->clearMoves();
   
    bool load1 = gpu_complex[gpu_n]->loadProtein(protein_name);
    
    if ( load1 )
    {
     cout << "Cannot read target protein structure" << endl;
     exit(EXIT_FAILURE);
    }
    
    bool load2 = gpu_complex[gpu_n]->loadParams(data_path);
    
    if ( load2 )
    {
     cout << "Cannot read parameter file" << endl;
     exit(EXIT_FAILURE);
    }
    
    bool load3 = gpu_complex[gpu_n]->loadLigand(lib1, lib2, molid);
    
    if ( load3 )
    {
     cout << "Cannot load ligand data" << endl;
     exit(EXIT_FAILURE);
    }
    
    bool load4 = gpu_complex[gpu_n]->loadLHM(lhm_name);
    
    if ( load4 )
    {
     cout << "Cannot read LHM potentials" << endl;
     exit(EXIT_FAILURE);
    }
    gpu_complex[gpu_n]->createTrialCoords();
    gpu_complex[gpu_n]->calculateEnergy();
    
    gpu_n++;
   }
   
   lib2 = 0;
  }
 }
 
 for ( int in = 0; in < gpu_n; in++ )
 {
  time(&t_start2);
  gpu_complex[in]->clearMoves();
  cout << "Ligand ID: " << gpu_complex[in]->getLigandID() << endl;
  gpu_complex[in]->createTrialCoords(); 
  cout << "Initial conformation:";
  cout<< fixed << setprecision(7) << setw(12) << gpu_complex[in]->getEnergy(1) << endl;
 
  gpu_complex[in]->restoreCoords(1);
///////////////////////////////////////////////////////////
  REMC( gpu_complex[in], remc_replicas, remc_steps, remc_cycles, step_t, step_r, dsphere, boltzmann, false ); // the variable simplex_steps not included
  gpu_complex[in]->clearMoves();
  gpu_complex[in]->setBest();
  gpu_complex[in]->createTrialCoords();
  gpu_complex[in]->calculateEnergy();
  cout<<"Best Energy: ";
  for(int k=1;k<8;k++)
   cout<<gpu_complex[in]->getEnergy(k)<<" "<<endl;
  

   for ( int rn = 1; rn <= remc_replicas; rn++ )
    cout << "Replica" << setw(4) << rn << " acceptance ratio:" << fixed << setprecision(6) << setw(10) << gpu_complex[in]->getAcceptanceENE(rn) << endl;

   cout << "Swap acceptance ratio:" << fixed << setprecision(6) << setw(17) << gpu_complex[in]->getAcceptanceSWP() << endl;
  
 
 
 
 
 
  time(&t_end2);
  
  double dif2 = difftime(t_end2, t_start2);
  
  cout << "Docking time: " << fixed << setprecision(0) << dif2 << " s" << endl;
 }
 
 
 
 
 
 
 time(&t_end1);
 
 cout << endl;
 
 printTime( difftime(t_end1, t_start1) );
 
 return 0;
}

