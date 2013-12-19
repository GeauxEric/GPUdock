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

#include<string>
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
#include<cstring>

#include "size.h"
#include "coords.h"
#include "rmsd.h"
#include "data.h"

using namespace std;

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
    int                        _pens_best;            // best conformation
    
    map<pair<int,int>, double> _complex_psp;          // pocket-specific potential
    
    
    
    
    
    
    
    /* LIGAND */
    
    vector<CoordsLigand> _ligand_xyz;         // ligand heavy-atom coordinates
    
    double               _ligand_center[3];   // ligand geometric center
    
    int                  _lna;                // number of ligand atoms
    int                  _lnb;                // number of ligand bonds
    
    int                  _lens_total;         // total ensemble conformations
    int                  _lens_current;       // current conformation
    int                  _lens_best;          // best conformation
    
    double               _lens_rmsd[MAXEN2];  // ensemble rmsd
    
    std::string          _ligand_id;          // ligand id
    std::string          _ligand_smiles;      // ligand smiles
    
    bitset<MAXFP1>       _ligand_fpt_smiles;  // fingerprint smiles
    bitset<MAXFP2>       _ligand_fpt_maccs;   // fingerprint maccs
    
    double               _ligand_prop_mw;     // ligand molecular weight
    double               _ligand_prop_logp;   // ligand water/octanol partition coeff
    double               _ligand_prop_psa;    // ligand polar surface area
    double               _ligand_prop_mr;     // ligand molar refractivity
    
    int                  _ligand_prop_hbd;    // ligand hydrogen bond donors
    int                  _ligand_prop_hba;    // ligand hydrogen bond acceptors
    
    double               _ligand_tra[3];      // ligand translation
    double               _ligand_rot[3];      // ligand rotation
    
    double               _ligand_tra_best[3]; // best ligand translation
    double               _ligand_rot_best[3]; // best ligand rotation
    
    list<mcs_restraints> _ligand_mcs;         // lhm position restraints
    
    
    
    
    
    
    
    
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
    
    double _ebst;                           // best energy
    
    double _par_lj[3];                      // L-J params
    double _par_el[2];                      // electrostatic params
    double _par_kde;                        // kde bandwidth
    
    double _temp;                           // current temperature
    
    
    
    
    
    
    
    
    
    
    /* STATS */
    
    int _stats_evals; // function evaluations
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  public:
    
    Complex( int, int );
    
    Complex( void );
    
    ~Complex();
    
    /* RECEPTOR */
    
    int getProteinPointsTotal( void );
    
    int getProteinResiduesTotal( void );
    
    int getProteinEnsembleTotal( void );
    
    int getProteinEnsembleCurrent( void );
    
    void setProteinEnsembleCurrent( int );
    
    std::string getProteinSequence( void );
    
    bool loadProtein( std::string );
    
    bool loadParams( std::string );
    
    bool loadLHM( std::string );
    
    
    
    
    /* LIGAND */
    
    int getLigandAtomsTotal( void );
    
    int getLigandBondsTotal( void );
    
    int getLigandEnsembleTotal( void );
    
    int getLigandEnsembleCurrent( void );
    
    void setLigandEnsembleCurrent( int );
    
    bool loadLigand( std::string [], int, std::string );
    
    std::string getLigandID( void );
    
    
    
    
    
    
    
    /* COMPLEX */
    
    void calculateEnergy( void );
    
    double getEnergy( int );
    
    void setTemperature( double );
    
    double getTemperature( void );
    
    void setConfiguration( double [] );
    
    void setBest( void );
    
    
    
    
    
    
    /* STATS */
    
    int getFunctionEvals( void );
    
    
    
    
    
    
    
    
};

#endif
