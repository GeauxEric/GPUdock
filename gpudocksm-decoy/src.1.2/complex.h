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


#ifndef __COMPLEX_H_
#define __COMPLEX_H_

#include<cstring>
#include<cstdio>
#include<vector>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<cmath>
#include<algorithm>
#include<list>
#include<bitset>
#include<map>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_randist.h>


#include "size.h"
#include "coords.h"
#include "rmsd.h"
#include "data.h"

using namespace std;


struct Param{
	int remc_steps;
	int remc_cycles;	    
	int remc_replicas;	    // total replica numbers
	double mc_temp;	    // temperature in each replica
	double step_t;
	double step_r;
	double dsphere;
	double boltzmann;

	int data_pt_total;
};

struct mcs_restraints
{
 double          tcc;
 int             total;
 map<int,double> xcoord;
 map<int,double> ycoord;
 map<int,double> zcoord;
};

class Complex {
        
  private:
    
    /* RECEPTOR */
    
    vector<CoordsProtein>      _protein_xyz;          // protein effective points coords
    vector<CoordsKDE>          _kde_xyz;              // kde coords
    
    double                     _pocket_center[3];     // pocket center
    
    int                        _pnp;                  // number of protein effective points
    int                        _pnr;                  // number of protein residues
    int                        _pnk;                  // total number of kde points
    
    int                        _pns[MAXTP2];          // number of specific kde points
    
    std::string                _protein_seq1;         // aa sequence
    char                       _protein_seq2[MAXPRO]; // aa sequence
    int                        _protein_seq3[MAXPRO]; // aa sequence numbering
    
    int                        _pens_total;           // total ensemble conformations
    int                        _pens_current;         // current conformation
    
    map<pair<int,int>, double> _complex_psp;          // pocket-specific potential
    
    
    
    
    
    
    
     /* LIGAND */
    
    vector<CoordsLigand> _ligand_xyz;        		// ligand heavy-atom coordinates
    //////////////////////////////////////////////////////////////////////////////////ADDED BY DANIEL
    vector<CoordsLigand> _ligand_bkp[2];            	// backup coordinates
    vector<CoordsLigandTrial> _ligand_trial;        	// A set of trial coordinates, to be used for energy calculation
    vector<CoordsLigand> _ligand_best;              	// lowest-energy coordinates
    double               _ligand_tmatrix[3][3];     	// Matrix which holds the change in coordinates due to a rotation
    double               _perturb_delta[MAXLIG][3]; 	// Holds the changes in coords due to a perturbation  
    double               _ref_coords[MAXLIG][3][MAXEN2];// Original conformation coordinates    
    double		 _ebest;		    	// best energy
    double   	 	 _cbest;		    	// best ensemble number
    
    double     		 _cdist;             		// distance between geom centers
    double		 _crmsd;	     		// RMSD
    double		 _cmcc;		     		// MCC
 
    /////////////////////////////////////////////////////////////////////////////////////////////////
    double               _ligand_center[3];  // ligand geometric center
    
    int                  _lna;               // number of ligand atoms
    int                  _lnb;               // number of ligand bonds
    double		 _temp;              // temperature               
    int                  _lens_total;        // total ensemble conformations
    int                  _lens_current;      // current conformation
    
    double               _lens_rmsd[MAXEN2]; // ensemble rmsd
    
    std::string          _ligand_id;         // ligand id
    std::string          _ligand_smiles;     // ligand smiles
    
    bitset<MAXFP1>       _ligand_fpt_smiles; // fingerprint smiles
    bitset<MAXFP2>       _ligand_fpt_maccs;  // fingerprint maccs
    
    double               _ligand_prop_mw;    // ligand molecular weight
    double               _ligand_prop_logp;  // ligand water/octanol partition coeff
    double               _ligand_prop_psa;   // ligand polar surface area
    double               _ligand_prop_mr;    // ligand molar refractivity
    
    int                  _ligand_prop_hbd;   // ligand hydrogen bond donors
    int                  _ligand_prop_hba;   // ligand hydrogen bond acceptors
    
    double               _ligand_tra[3];     // ligand translation
    double               _ligand_rot[3];     // ligand rotation
    
    list<mcs_restraints> _ligand_mcs;        // lhm position restraints
   /////////////////////////////////////////////////////////////////////////////////DANIEL
    double    		 _acc_ene[MAXREP][2];// MC acceptance Ratio
    double 		 _acc_swp[2]; 	     // swap acceptance ratio
   ///////////////////////////////////////////////////////////////////////////////// 
    
    
    
    
    
    
    
     /* COMPLEX */
    
    double _complex_vdw[MAXTP1][MAXTP2][2]; // L-J potential
    double _complex_ele[MAXTP3];            // electrostatic potential
    double _complex_pmf[MAXTP1][MAXTP2][2]; // contact potential
    double _complex_hpp[MAXTP4];            // protein hydrophobicity
    double _complex_hpl[MAXTP2][2];         // ligand hydrophobicity
    double _complex_hdb[MAXTP1][MAXTP2][2]; // ligand hydrophobicity
    
    double _weights[MAXWEI];                // weights for energy terms
    
    double _etot;                           // total energy
    double _evdw;                           // soft L-J potential
    double _eele;                           // soft electrostatic potential
    double _epmf;                           // contact potential
    double _ehpc;                           // hydrophobic potential
    double _ehdb;                           // hydrogen bond potential
    double _edst;                           // distance restraints
    double _epsp;                           // pocket-specific potential
    double _ekde;                           // kde potential
    double _elhm;                           // position restraints
 
    double _par_lj[3];                      // L-J params
    double _par_el[2];                      // electrostatic params
    double _par_kde;                        // kde bandwidth
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  public:
    
    Complex( int, int );
    
    Complex( void );
    
    ~Complex();
    
     /* RECEPTOR */
    
    int getProteinPointsTotal( void );
    
    int getProteinResiduesTotal( void );
    
    int getProteinEnsembleTotal( void );
    
    int getProteinEnsembleCurrent( void );
    
    std::string getProteinSequence( void );
    
    bool loadProtein( std::string );
    
    bool loadParams( std::string );
    
    bool loadLHM( std::string );
    //////////////////////////////////////////
    void protEnsemble(int);

    void setProteinEnsembleCurrent( int );

    void setLigandEnsembleCurrent( int );
    //////////////////////////////////////////
    
    
    
     /* LIGAND */
    
    int getLigandAtomsTotal( void );
    
    int getLigandBondsTotal( void );
    
    int getLigandEnsembleTotal( void );
    
    int getLigandEnsembleCurrent( void );
    
    bool loadLigand( std::string [], int, std::string );
    
    std::string getLigandID( void );
    
    /////////////////////////////////// Added by DANIEL//////
    void ligTranslate(int, double );

    void ligRotate(double, double, double);

    void ligPerturb(int);

    void ligEnsemble(int);
    ////////////////////////////////////////////////////////   
    
    
    ////////////////////////////////////////////////////////   by Eric
    ////////////////////////////////////////////////////////   

    
    
    
    /* COMPLEX */
    
    void calculateEnergy( void );
    
    double getEnergy( int );

    ////////////////////////////////// Added by Daniel//////
    double getMCC( void);

    void backupCoords ( int );

    void restoreCoords ( int );

    void createTrialCoords ( void );

    void getConfCoords ( double [][3]);

    void getConfCoordsPerm ( double [][3]);

    void getConfParams ( double [] );
 
    void setConfParams ( double [] );
  
    void setConfCoords ( double [][3] );

    void setConfCoordsPerm ( double [][3] );

    void periodicBoundary ( double );

    double getTemperature ( void );

    void setTemperature ( double );

    void clearMoves ( void );

    double getDistance (int);

    void addAcceptanceENE( int, bool );

    void addAcceptanceSWP( bool );

    double getAcceptanceENE( int );

    double getAcceptanceSWP( void );
    
    void setBest( void);

    void printLigandInfo(ofstream&, ofstream&, ofstream&);

    void printEnerInfo(ofstream & eout);

    void printProteinInfo(ofstream&);
 
    ////////////////////////////////////////////////////////
    
    
    
};

struct OutStream{
    int all;	    // number of points in gapless data set
    int low;	    // number of points in low data set
    int high;	    // number of points in high data set

    bool all_stop;
    bool low_stop;
    bool high_stop;

    ofstream all_decoys;
    ofstream all_MCcoef;
    ofstream all_decoyenergies;

    ofstream L_decoys;
    ofstream L_MCcoef;
    ofstream L_decoyenergies;

    ofstream H_decoys;
    ofstream H_MCcoef;
    ofstream H_decoyenergies;
};

void initOutStream( OutStream *outfiles, std::string lig_id );

void closeOutStream(OutStream *outfiles);

#endif
