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


#include "complex.h"

using namespace std;

Complex::Complex( int ap, int al )
{
 _protein_xyz.clear();
 _kde_xyz.clear();
 _complex_psp.clear();
 
 _pnp = ap;
 _pnr = ap;
 _lna = al;
 _lnb = al;
 
 _pnk = 0;
 
 for ( int ai = 0; ai < MAXTP2; ai++ )
  _pns[ai] = 0;
 
 for ( int ai = 0; ai < 3; ai++ )
 {
  _pocket_center[ai] = 0.0;
  _ligand_center[ai] = 0.0;
 }
 
 _pens_total = 0;
 _pens_current = 0;
 _pens_best = 0;
 
 _lens_total = 0;
 _lens_current = 0;
 _lens_best = 0;
 
 _ligand_id = "";
 _ligand_smiles = "";
 
 _ligand_prop_mw = 0.0;
 _ligand_prop_logp = 0.0;
 _ligand_prop_psa = 0.0;
 _ligand_prop_mr = 0.0;
 
 _ligand_prop_hbd = 0;
 _ligand_prop_hba = 0;
 
 _temp = 0.0;
 
 for ( int ai = 0; ai < MAXEN2; ai++ )
  _lens_rmsd[ai] = 0.0;
 
 _ligand_fpt_smiles.reset();
 _ligand_fpt_maccs.reset();
 
 _protein_seq1 = "";
 
 for ( int ai = 0; ai < MAXTP1; ai++ )
  for ( int aj = 0; aj < MAXTP2; aj++ )
   for ( int ak = 0; ak < 2; ak++ )
   {
    _complex_vdw[ai][aj][ak] = 0.0;
    _complex_pmf[ai][aj][ak] = 0.0;
    
    _complex_hdb[ai][aj][ak] = 0.0;
   }
 
 for ( int ai = 0; ai < MAXTP3; ai++ )
  _complex_ele[ai] = 0.0;
 
 for ( int ai = 0; ai < MAXTP4; ai++ )
  _complex_hpp[ai] = 0.0;
 
 for ( int ai = 0; ai < MAXTP2; ai++ )
  for ( int aj = 0; aj < 2; aj++ )
   _complex_hpl[ai][aj] = 0.0;
 
 for ( int ai = 0; ai < MAXWEI; ai++ )
  _weights[ai] = 0.0;
 
 for ( int ai = 0; ai < 3; ai++ )
  _par_lj[ai] = 0.0;
 
 for ( int ai = 0; ai < 2; ai++ )
  _par_el[ai] = 0.0;
 
 _par_kde = 0.0;
 
 for ( int ai = 0; ai < 3; ai++ )
 {
  _ligand_tra[ai] = 0.0;
  _ligand_rot[ai] = 0.0;
  
  _ligand_tra_best[ai] = 0.0;
  _ligand_rot_best[ai] = 0.0;
 }
 
 _ebst = 1e6;
 
 _stats_evals = 0;
}

Complex::Complex( void )
{
 _protein_xyz.clear();
 _kde_xyz.clear();
 _complex_psp.clear();
 
 _pnp = 0;
 _pnr = 0;
 _lna = 0;
 _lnb = 0;
 
 _pnk = 0;
 
 for ( int ai = 0; ai < MAXTP2; ai++ )
  _pns[ai] = 0;
 
 for ( int ai = 0; ai < 3; ai++ )
 {
  _pocket_center[ai] = 0.0;
  _ligand_center[ai] = 0.0;
 }
 
 _pens_total = 0;
 _pens_current = 0;
 _pens_best = 0;
 
 _lens_total = 0;
 _lens_current = 0;
 _lens_best = 0;
 
 _ligand_id = "";
 _ligand_smiles = "";
 
 _ligand_prop_mw = 0.0;
 _ligand_prop_logp = 0.0;
 _ligand_prop_psa = 0.0;
 _ligand_prop_mr = 0.0;
 
 _ligand_prop_hbd = 0;
 _ligand_prop_hba = 0;
 
 _temp = 0.0;
 
 for ( int ai = 0; ai < MAXEN2; ai++ )
  _lens_rmsd[ai] = 0.0;
 
 _ligand_fpt_smiles.reset();
 _ligand_fpt_maccs.reset();
 
 _protein_seq1 = "";
 
 for ( int ai = 0; ai < MAXTP1; ai++ )
  for ( int aj = 0; aj < MAXTP2; aj++ )
   for ( int ak = 0; ak < 2; ak++ )
   {
    _complex_vdw[ai][aj][ak] = 0.0;
    _complex_pmf[ai][aj][ak] = 0.0;
    
    _complex_hdb[ai][aj][ak] = 0.0;
   }
 
 for ( int ai = 0; ai < MAXTP3; ai++ )
  _complex_ele[ai] = 0.0;
 
 for ( int ai = 0; ai < MAXTP4; ai++ )
  _complex_hpp[ai] = 0.0;
 
 for ( int ai = 0; ai < MAXTP2; ai++ )
  for ( int aj = 0; aj < 2; aj++ )
   _complex_hpl[ai][aj] = 0.0;
 
 _par_kde = 0.0;
 
 for ( int ai = 0; ai < MAXWEI; ai++ )
  _weights[ai] = 0.0;
 
 for ( int ai = 0; ai < 3; ai++ )
  _par_lj[ai] = 0.0;
 
 for ( int ai = 0; ai < 2; ai++ )
  _par_el[ai] = 0.0;
 
 for ( int ai = 0; ai < 3; ai++ )
 {
  _ligand_tra[ai] = 0.0;
  _ligand_rot[ai] = 0.0;
  
  _ligand_tra_best[ai] = 0.0;
  _ligand_rot_best[ai] = 0.0;
 }
 
 _ebst = 1e6;
 
 _stats_evals = 0;
}

Complex::~Complex() {}

// ==================================================================================   getProteinPointsTotal

int Complex::getProteinPointsTotal( void )
{
 return _pnp;
}

// ==================================================================================   getProteinResiduesTotal

int Complex::getProteinResiduesTotal( void )
{
 return _pnr;
}

// ==================================================================================   getProteinEnsembleTotal

int Complex::getProteinEnsembleTotal( void )
{
 return _pens_total;
}

// ==================================================================================   getProteinEnsembleCurrent

int Complex::getProteinEnsembleCurrent( void )
{
 return _pens_current;
}

// ==================================================================================   setProteinEnsembleCurrent

void Complex::setProteinEnsembleCurrent( int ei )
{
 _pens_current = ei;
}

// ==================================================================================   getProteinSequence

std::string Complex::getProteinSequence( void )
{
 return _protein_seq1;
}

// ==================================================================================   getLigandAtomsTotal

int Complex::getLigandAtomsTotal( void )
{
 return _lna;
}

// ==================================================================================   getLigandBondsTotal

int Complex::getLigandBondsTotal( void )
{
 return _lnb;
}

// ==================================================================================   getLigandEnsembleTotal

int Complex::getLigandEnsembleTotal( void )
{
 return _lens_total;
}

// ==================================================================================   getLigandEnsembleCurrent

int Complex::getLigandEnsembleCurrent( void )
{
 return _lens_current;
}

// ==================================================================================   setLigandEnsembleCurrent

void Complex::setLigandEnsembleCurrent( int ei )
{
 _lens_current = ei;
}

// ==================================================================================   loadProtein

bool Complex::loadProtein( std::string p1_name )
{
 list<string> p1_data;
 list<string>::iterator p1_i;
 
 string line1;
 
 ifstream p1_file( p1_name.c_str() );
 
 if ( !p1_file.is_open() )  { cout << "Cannot open " << p1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(p1_file,line1))
  p1_data.push_back( line1 );
 
 p1_file.close();
 
 for ( p1_i = p1_data.begin(); p1_i != p1_data.end(); p1_i++ )
 {
  if ( (*p1_i).size() > 53 )
  {
   if ( (*p1_i).substr(0,6) == "ATOM  " )
   {
    std::string atom1 = (*p1_i).substr(12,4);
    
    int residue1 = getResCode( (*p1_i).substr(17,3) );
    
    int residue2 = atoi((*p1_i).substr(22,4).c_str());
    
    if ( atom1 == " CA " )
    {
     _protein_seq1.append( three2oneS( (*p1_i).substr(17,3) ) );
     
     _protein_seq3[_pnr] = residue2;
     
     _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 0, residue1, 0 ) );
     
     if ( _pnr > 1 )
      _protein_xyz.push_back( CoordsProtein( _pnr - 1, _pnp++, 1, residue1, 1 ) );
     
     switch ( residue1 )
     {
      case  0 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++,  2, residue1, 2 ) ); break;
      case  1 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 21, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 22, residue1, 3 ) ); break;
      case  2 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++,  3, residue1, 2 ) ); break;
      case  3 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 18, residue1, 2 ) ); break;
      case  4 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 14, residue1, 2 ) ); break;
      case  5 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 11, residue1, 2 ) ); break;
      case  6 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++,  5, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++,  6, residue1, 3 ) ); break;
      case  7 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++,  9, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 10, residue1, 3 ) ); break;
      case  8 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 26, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 27, residue1, 3 ) ); break;
      case  9 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 24, residue1, 2 ) ); break;
      case 10 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 15, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 16, residue1, 3 ) ); break;
      case 11 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++,  4, residue1, 2 ) ); break;
      case 13 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++,  7, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++,  8, residue1, 3 ) ); break;
      case 14 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 12, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 13, residue1, 3 ) ); break;
      case 15 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 19, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 20, residue1, 3 ) ); break;
      case 16 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 23, residue1, 2 ) ); break;
      case 17 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 17, residue1, 2 ) ); break;
      case 18 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 28, residue1, 2 ) ); _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 29, residue1, 3 ) ); break;
      case 19 : _protein_xyz.push_back( CoordsProtein( _pnr, _pnp++, 25, residue1, 2 ) ); break;
     }
     
     _pnr++;
    }
   }
  }
  else if ( (*p1_i).substr(0,6) == "ENDMDL" )
  {
   break;
  }
 }
 
 strcpy(_protein_seq2, _protein_seq1.c_str());
 
 vector<CoordsProtein>::iterator p2_i;
 
 for ( p2_i = _protein_xyz.begin(); p2_i < _protein_xyz.end(); p2_i++ )
 {
  _pens_total = 0;
  
  int residue3 = _protein_seq3[(*p2_i).getResidueNumber()];
  int residue5 = (*p2_i).getResidueNumber();
  int point1 = (*p2_i).getPointType();
  
  double tx1 = 0.0;
  double ty1 = 0.0;
  double tz1 = 0.0;
  double tn1 = 0.0;
  
  for ( p1_i = p1_data.begin(); p1_i != p1_data.end(); p1_i++ )
  {
   if ( (*p1_i).size() > 53 )
   {
    if ( (*p1_i).substr(0,6) == "ATOM  " )
    {
     int residue4 = atoi((*p1_i).substr(22,4).c_str());
     
     std::string atom2 = (*p1_i).substr(12,4);
     
     double tx2 = atof((*p1_i).substr(30,8).c_str());
     double ty2 = atof((*p1_i).substr(38,8).c_str());
     double tz2 = atof((*p1_i).substr(46,8).c_str());
     
     if ( residue4 == residue3 )
     {
      if ( point1 == 0 && atom2 == " CA " )
      {
       tx1 += tx2;
       ty1 += ty2;
       tz1 += tz2;
       
       tn1 += 1.0;
      }
      else if ( point1 == 2 || point1 == 3 || point1 == 4 || point1 == 11 || point1 == 14 || point1 == 17 || point1 == 18 || point1 == 23 || point1 == 24 || point1 == 25 )
      {
       if ( atom2 != " N  " && atom2 != " CA " && atom2 != " C  " && atom2 != " O  " )
       {
        tx1 += tx2;
        ty1 += ty2;
        tz1 += tz2;
        
        tn1 += 1.0;
       }
      }
      else if ( point1 == 5 || point1 == 7 || point1 == 9 || point1 == 12 || point1 == 15 || point1 == 19 || point1 == 21 || point1 == 26 || point1 == 28 )
      {
       if ( atom2 != " N  " && atom2 != " CA " && atom2 != " C  " && atom2 != " O  " && atom2 != " CB " && atom2 != " CG " )
       {
        tx1 += tx2;
        ty1 += ty2;
        tz1 += tz2;
        
        tn1 += 1.0;
       }
      }
      else if ( point1 == 6 || point1 == 8 || point1 == 10 || point1 == 13 || point1 == 16 || point1 == 20 || point1 == 22 || point1 == 27 || point1 == 29 )
      {
       if ( atom2 == " CB " || atom2 == " CG " )
       {
        tx1 += tx2;
        ty1 += ty2;
        tz1 += tz2;
        
        tn1 += 1.0;
       }
      }
     }
     
     if ( point1 == 1 && residue5 > 0 && residue5 < _pnr - 1 )
     {
      int residue6 = _protein_seq3[(*p2_i).getResidueNumber()-1];
      
      if ( ( residue4 == residue3 && atom2 == " N  " ) || ( residue4 == residue6 && atom2 == " C  " ) || ( residue4 == residue6 && atom2 == " O  " ) )
      {
       tx1 += tx2;
       ty1 += ty2;
       tz1 += tz2;
       
       tn1 += 1.0;
      }
     }
    }
   }
   else if ( (*p1_i).substr(0,6) == "ENDMDL" )
   {
    if ( tn1 > 0.0 )
    {
     tx1 /= tn1;
     ty1 /= tn1;
     tz1 /= tn1;
     
     (*p2_i).setCoords(tx1, ty1, tz1, _pens_total);
    }
    
    tx1 = 0.0;
    ty1 = 0.0;
    tz1 = 0.0;
    tn1 = 0.0;
    
    _pens_total++;
   }
  }
 }
 
 return EXIT_SUCCESS;
}

// ==================================================================================   loadLHM

bool Complex::loadLHM( std::string h1_name )
{
 string line1;
 
 ifstream h1_file( h1_name.c_str() );
 
 if ( !h1_file.is_open() )  { cout << "Cannot open " << h1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(h1_file,line1))
 {
  if ( line1.length() > 3 )
  {
   if ( line1.substr(0,3) == "KDE" )
   {
    std::string dat1[5];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    _kde_xyz.push_back( CoordsKDE( _pnk++, getLigCode(dat1[1]), atof(dat1[2].c_str()), atof(dat1[3].c_str()), atof(dat1[4].c_str()) ) );
    
    _pns[getLigCode(dat1[1])]++;
   }
   else if ( line1.substr(0,3) == "PSP" )
   {
    std::string dat1[4];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    _complex_psp[make_pair( atoi(dat1[1].c_str()), getLigCode(dat1[2]) )] = atof(dat1[3].c_str());
   }
   else if ( line1.substr(0,3) == "MCS" )
   {
    std::string dat1[MAXMCS];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    if ( dat1[1] == _ligand_id )
    {
     mcs_restraints tmp_mcs;
     
     tmp_mcs.tcc = atof(dat1[2].c_str());
     tmp_mcs.total = atoi(dat1[3].c_str());
     
     for ( int ia = 0; ia < tmp_mcs.total; ia++ )
     {
      tmp_mcs.xcoord.insert( pair<int,double>(atoi(dat1[ia * 4 + 4].c_str()), atof(dat1[ia * 4 + 5].c_str())) );
      tmp_mcs.ycoord.insert( pair<int,double>(atoi(dat1[ia * 4 + 4].c_str()), atof(dat1[ia * 4 + 6].c_str())) );
      tmp_mcs.zcoord.insert( pair<int,double>(atoi(dat1[ia * 4 + 4].c_str()), atof(dat1[ia * 4 + 7].c_str())) );
     }
     
     if ( _ligand_mcs.size() < (int) MAXPOS )
      _ligand_mcs.push_back( tmp_mcs );
    }
   }
  }
  
  if ( _pnk >= MAXKDE )
   break;
 }
 
 h1_file.close();
 
 return EXIT_SUCCESS;
}

// ==================================================================================   loadParams

bool Complex::loadParams( std::string d1_name )
{
 string line1;
 
 ifstream d1_file( d1_name.c_str() );
 
 if ( !d1_file.is_open() )  { cout << "Cannot open " << d1_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(d1_file,line1))
  if ( line1.length() > 3 )
  {
   if ( line1.substr(0,3) == "VDW" )
   {
    std::string dat1[5];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    _complex_vdw[getPntCode(dat1[1])][getLigCode(dat1[2])][0] = atof(dat1[3].c_str());
    _complex_vdw[getPntCode(dat1[1])][getLigCode(dat1[2])][1] = atof(dat1[4].c_str());
   }
   else if ( line1.substr(0,3) == "PLJ" )
   {
    std::string dat1[4];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    for ( int i1 = 0; i1 < 3; i1++ )
     _par_lj[i1] = atof(dat1[i1+1].c_str());
   }
   else if ( line1.substr(0,3) == "ELE" )
   {
    std::string dat1[3];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    if ( dat1[1].substr(1,1) != "C" )
     _complex_ele[getPntCode(dat1[1])] = atof(dat1[2].c_str());
    else
     _complex_ele[getResCodeOne(dat1[1].substr(0,1))+30] = atof(dat1[2].c_str());
   }
   else if ( line1.substr(0,3) == "PEL" )
   {
    std::string dat1[3];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    for ( int i1 = 0; i1 < 2; i1++ )
     _par_el[i1] = atof(dat1[i1+1].c_str());
   }
   else if ( line1.substr(0,3) == "PMF" )
   {
    std::string dat1[5];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    _complex_pmf[getPntCode(dat1[1])][getLigCode(dat1[2])][0] = atof(dat1[3].c_str());
    _complex_pmf[getPntCode(dat1[1])][getLigCode(dat1[2])][1] = atof(dat1[4].c_str());
   }
   else if ( line1.substr(0,3) == "HPP" )
   {
    std::string dat1[3];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    _complex_hpp[getResCodeOne(dat1[1])] = atof(dat1[2].c_str());
   }
   else if ( line1.substr(0,3) == "HPL" )
   {
    std::string dat1[4];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    _complex_hpl[getLigCode(dat1[1])][0] = atof(dat1[2].c_str());
    _complex_hpl[getLigCode(dat1[1])][1] = atof(dat1[3].c_str());
   }
   else if ( line1.substr(0,3) == "HDB" )
   {
    std::string dat1[5];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    _complex_hdb[getPntCode(dat1[1])][getLigCode(dat1[2])][0] = atof(dat1[3].c_str());
    _complex_hdb[getPntCode(dat1[1])][getLigCode(dat1[2])][1] = atof(dat1[4].c_str());
   }
   else if ( line1.substr(0,3) == "WEI" )
   {
    std::string dat1[10];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    for ( int wi = 0; wi < MAXWEI; wi++ )
     _weights[wi] = atof(dat1[wi+1].c_str());
   }
   else if ( line1.substr(0,3) == "KDE" )
   {
    std::string dat1[2];
    
    int dat2 = 0;
    
    istringstream dat3(line1);
    
    while (dat3)
     dat3 >> dat1[dat2++];
    
    _par_kde = atof(dat1[1].c_str());
   }
  }
 
 d1_file.close();
 
 /* normalize hydrophobic scale */
 
 double hpc1 = 1e6;
 double hpc2 = -1e6;
 
 for ( int ai = 0; ai < MAXTP4; ai++ )
 {
  if ( _complex_hpp[ai] < hpc1 )
   hpc1 = _complex_hpp[ai];
  
  if ( _complex_hpp[ai] > hpc2 )
   hpc2 = _complex_hpp[ai];
 }
 
 for ( int ai = 0; ai < MAXTP4; ai++ )
 {
  _complex_hpp[ai] = ( ( _complex_hpp[ai] - hpc1 ) / ( hpc2 - hpc1 ) ) * 2.0 - 1.0;
  
  if ( _complex_hpp[ai] < -1.0 )
   _complex_hpp[ai] = -1.0;
  
  if ( _complex_hpp[ai] >  1.0 )
   _complex_hpp[ai] =  1.0;
 }
 
 return EXIT_SUCCESS;
}

// ==================================================================================   loadLigand

bool Complex::loadLigand( std::string llib1[], int llib2, std::string llib3 )
{
 _lna = atoi(llib1[3].substr(0,3).c_str());
 _lnb = atoi(llib1[3].substr(3,3).c_str());
 
 double tmp1[MAXEN2][MAXLIG][3];
 double tmp2[MAXLIG];
 int tmp3[MAXLIG];
 
 int tmp7 = 0;
 
 for ( int i1 = 4 + _lna + _lnb; i1 < llib2 - 1; i1++ )
 {
  if ( llib1[i1].find(llib3) != string::npos )
   _ligand_id = llib1[i1+1];
  
  else if ( llib1[i1].find("SMILES_CANONICAL") != string::npos )
   _ligand_smiles = llib1[i1+1];
  
  else if ( llib1[i1].find("OB_MW") != string::npos )
   _ligand_prop_mw = atof(llib1[i1+1].c_str());
  
  else if ( llib1[i1].find("OB_logP") != string::npos )
   _ligand_prop_logp = atof(llib1[i1+1].c_str());
  
  else if ( llib1[i1].find("OB_PSA") != string::npos )
   _ligand_prop_psa = atof(llib1[i1+1].c_str());
  
  else if ( llib1[i1].find("OB_MR") != string::npos )
   _ligand_prop_mr = atof(llib1[i1+1].c_str());
  
  else if ( llib1[i1].find("MCT_HBD") != string::npos )
   _ligand_prop_hbd = atoi(llib1[i1+1].c_str());
  
  else if ( llib1[i1].find("MCT_HBA") != string::npos )
   _ligand_prop_hba = atoi(llib1[i1+1].c_str());
  
  else if ( llib1[i1].find("FINGERPRINT") != string::npos )
  {
   string::iterator i2;
   int i3;
   
   for ( i2 = llib1[i1+1].begin(), i3 = 0; i2 < llib1[i1+1].end(); i2++, i3++ )
    switch(*i2)
    {
     case '0': _ligand_fpt_smiles[MAXFP1-1-i3] = 0; break;
     case '1': _ligand_fpt_smiles[MAXFP1-1-i3] = 1; break;
    }
  }
  
  else if ( llib1[i1].find("MACCS166") != string::npos )
  {
   string::iterator i2;
   int i3;
   
   for ( i2 = llib1[i1+1].begin(), i3 = 0; i2 < llib1[i1+1].end(); i2++, i3++ )
    switch(*i2)
    {
     case '0': _ligand_fpt_maccs[MAXFP2-1-i3] = 0; break;
     case '1': _ligand_fpt_maccs[MAXFP2-1-i3] = 1; break;
    }
  }
  
  else if ( llib1[i1].find("OB_ATOM_TYPES") != string::npos )
  {
   int tmp4 = 0;
   
   istringstream tmp5(llib1[i1+1]);
   
   while (tmp5)
   {
    std::string tmp6;
    
    tmp5 >> tmp6;
    
    if ( tmp6.length() > 0 )
     tmp3[tmp4++] = getLigCode(tmp6);
   }
  }
  
  else if ( llib1[i1].find("OB_ATOMIC_CHARGES") != string::npos )
  {
   int tmp4 = 0;
   
   istringstream tmp5(llib1[i1+1]);
   
   while (tmp5)
   {
    std::string tmp6;
    
    tmp5 >> tmp6;
    
    if ( tmp6.length() > 0 )
     tmp2[tmp4++] = atof(tmp6.c_str());
   }
  }
  
  else if ( llib1[i1].find("ENSEMBLE_COORDS") != string::npos )
   while ( llib1[++i1].size() || tmp7 < MAXEN2 )
   {
    int tmp8 = 0;
    int tmp9 = 0;
    
    istringstream tmp5(llib1[i1]);
    
    while (tmp5)
    {
     std::string tmp6;
     
     tmp5 >> tmp6;
     
     if ( tmp6.length() > 0 )
     {
      tmp1[tmp7][tmp8][tmp9] = atof(tmp6.c_str());
      
      if ( ++tmp9 > 2 )
      {
       tmp8++;
       
       tmp9 = 0;
      }
     }
    }
    
    tmp7++;
   }
 }
 
 double tmp8[MAXLIG][3];
 
 for ( int i1 = 4; i1 < _lna + 4; i1++ )
 {
  _ligand_xyz.push_back( CoordsLigand( i1 - 4, llib1[i1].substr(31,24).c_str(), tmp3[i1 - 4], tmp2[i1 - 4] ) );
  
  tmp8[i1-4][0] = atof(llib1[i1].substr(0,10).c_str());
  tmp8[i1-4][1] = atof(llib1[i1].substr(10,10).c_str());
  tmp8[i1-4][2] = atof(llib1[i1].substr(20,10).c_str());
  
  for ( int i5 = 0; i5 < 3; i5++ )
  {
   _ligand_center[i5] += tmp8[i1-4][i5];
   _pocket_center[i5] += tmp8[i1-4][i5];
  }
 }
 
 for ( int i5 = 0; i5 < 3; i5++ )
 {
  _ligand_center[i5] /= (double) _lna;
  _pocket_center[i5] /= (double) _lna;
 }
 
 for ( int i1 = 0; i1 < _lna; i1++ )
  for ( int i5 = 0; i5 < 3; i5++ )
   tmp8[i1][i5] -= _ligand_center[i5];
 
 vector<CoordsLigand>::iterator i4;
 
 for ( i4 = _ligand_xyz.begin(); i4 < _ligand_xyz.end(); i4++ ){
  (*i4).setCoords( tmp8[(*i4).getAtomNumber()][0], tmp8[(*i4).getAtomNumber()][1], tmp8[(*i4).getAtomNumber()][2], _lens_total );
   
  // cout << tmp8[(*i4).getAtomNumber()][0] << endl;
  // cout << tmp8[(*i4).getAtomNumber()][1] << endl;
  // cout << tmp8[(*i4).getAtomNumber()][2] << endl;
 }
 
 _lens_rmsd[_lens_total++] = 0.0;
 
 for ( int i6 = 0; i6 < tmp7; i6++ )
 {
  double weights[MAXLIG];
  double ref_xyz[MAXLIG][3];
  double mob_xyz[MAXLIG][3];
  
  for ( int i7 = 0; i7 < _lna; i7++ )
  {
   weights[i7] = 1.0;
   
   for ( int i5 = 0; i5 < 3; i5++ )
   {
    ref_xyz[i7][i5] = tmp8[i7][i5];
    
    mob_xyz[i7][i5] = tmp1[i6][i7][i5];
   }
  }
  
  double mob_cen[3];
  
  for ( int i5 = 0; i5 < 3; i5++ )
   mob_cen[i5] = 0.0;
  
  for ( int i7 = 0; i7 < _lna; i7++ )
   for ( int i5 = 0; i5 < 3; i5++ )
    mob_cen[i5] += mob_xyz[i7][i5];
  
  for ( int i5 = 0; i5 < 3; i5++ )
   mob_cen[i5] /= (double) _lna;
  
  for ( int i7 = 0; i7 < _lna; i7++ )
   for ( int i5 = 0; i5 < 3; i5++ )
    mob_xyz[i7][i5] -= mob_cen[i5];
  
  int mode = 1;
  double rms1 = 0.0;
  double u[3][3];
  double t[3];
  int ier = 0;
  
  u3b_(&weights, &mob_xyz, &ref_xyz, &_lna, &mode, &rms1, &u, &t, &ier);
  
  double rms2 = sqrt( rms1 / (double) _lna );
  
  if ( rms2 > 0.1 )
  {
   for ( i4 = _ligand_xyz.begin(); i4 < _ligand_xyz.end(); i4++ )
   (*i4).setCoords( mob_xyz[(*i4).getAtomNumber()][0], mob_xyz[(*i4).getAtomNumber()][1], mob_xyz[(*i4).getAtomNumber()][2], _lens_total );
   
   _lens_rmsd[_lens_total++] = rms2;
  }
  
  if ( _lens_total >= MAXEN2 )
   break;
 }
 
 return EXIT_SUCCESS;
}

// ==================================================================================   getLigandID

std::string Complex::getLigandID( void )
{
 return _ligand_id;
}

// ==================================================================================   calculateEnergy

void Complex::calculateEnergy( void )
{
 _etot = 0.0;
 _evdw = 0.0;
 _eele = 0.0;
 _epmf = 0.0;
 _ehpc = 0.0;
 _ehdb = 0.0;
 _edst = 0.0;
 _epsp = 0.0;
 _ekde = 0.0;
 _elhm = 0.0;
 
 for ( int il5 = 0; il5 < 3; il5++ )
  _ligand_center[il5] = 0.0;
 
 list<double> rmsd1[MAXPOS];
 
 list<mcs_restraints>::iterator ir1;
 
 int ir2;
 
 for ( ir1 = _ligand_mcs.begin(), ir2 = 0; ir1 != _ligand_mcs.end(); ir1++, ir2++ )
  rmsd1[ir2].push_back( (*ir1).tcc );
 
 int ir4 = ir2;
 
 vector<CoordsProtein>::iterator ip1;
 vector<CoordsLigand>::iterator il1;
 
 for ( il1 = _ligand_xyz.begin(); il1 < _ligand_xyz.end(); il1++ )
 {
  double b_xyz[3], t_xyz[3], p_xyz[3], r_mat[3][3][3];
  
  r_mat[0][0][0] = 1.;                      r_mat[0][0][1] = 0.;                      r_mat[0][0][2] = 0.;
  r_mat[0][1][0] = 0.;                      r_mat[0][1][1] = cos(_ligand_rot[0]);     r_mat[0][1][2] = sin(_ligand_rot[0]);
  r_mat[0][2][0] = 0.;                      r_mat[0][2][1] = -1.*sin(_ligand_rot[0]); r_mat[0][2][2] = cos(_ligand_rot[0]);
  
  r_mat[1][0][0] = cos(_ligand_rot[1]);     r_mat[1][0][1] = 0.;                      r_mat[1][0][2] = -1.*sin(_ligand_rot[1]);
  r_mat[1][1][0] = 0.;                      r_mat[1][1][1] = 1.;                      r_mat[1][1][2] = 0.;
  r_mat[1][2][0] = sin(_ligand_rot[1]);     r_mat[1][2][1] = 0.;                      r_mat[1][2][2] = cos(_ligand_rot[1]);
  
  r_mat[2][0][0] = cos(_ligand_rot[2]);     r_mat[2][0][1] = sin(_ligand_rot[2]);     r_mat[2][0][2] = 0.;
  r_mat[2][1][0] = -1.*sin(_ligand_rot[2]); r_mat[2][1][1] = cos(_ligand_rot[2]);     r_mat[2][1][2] = 0.;
  r_mat[2][2][0] = 0.;                      r_mat[2][2][1] = 0.;                      r_mat[2][2][2] = 1.;
  
  for ( int il2 = 0; il2 < 3; il2++)
   b_xyz[il2] = (*il1).getCoords( il2 + 1, _lens_current );
  
  for ( int il4 = 0; il4 < 3; il4++)
  {
   for ( int il2 = 0; il2 < 3; il2++)
   {
    t_xyz[il2] = 0.0;
    
    for ( int il3 = 0; il3 < 3; il3++)
     t_xyz[il2] += b_xyz[il3] * r_mat[il4][il3][il2];
   }
   
   for ( int il2 = 0; il2 < 3; il2++)
    b_xyz[il2] = t_xyz[il2];
  }
  
  for ( int il2 = 0; il2 < 3; il2++)
  {
   b_xyz[il2] = b_xyz[il2] + _ligand_tra[il2] + _pocket_center[il2];
   
   _ligand_center[il2] += b_xyz[il2];
  }
  
  double hpc1 = 0.0;
  
  for ( ip1 = _protein_xyz.begin(); ip1 < _protein_xyz.end(); ip1++ )
  {
   for ( int ip2 = 0; ip2 < 3; ip2++)
    p_xyz[ip2] = (*ip1).getCoords( ip2 + 1, _pens_current );
   
   double dst = sqrt( pow( b_xyz[0] - p_xyz[0], 2 ) +
                      pow( b_xyz[1] - p_xyz[1], 2 ) +
                      pow( b_xyz[2] - p_xyz[2], 2 )   );
   
   /* L-J potential */
   
   double r1 = _complex_vdw[(*ip1).getPointType()][(*il1).getAtomType()][0];
   double e1 = _complex_vdw[(*ip1).getPointType()][(*il1).getAtomType()][1];
   
   double p1 = ( 2.0 * e1 * pow( _par_lj[2] * r1, 9 ) ) / ( pow( dst, 9 ) );
   
   double p2 = ( 3.0 * e1 * pow( _par_lj[2] * r1, 6 ) ) / ( pow( dst, 6 ) );
   
   double p3 = p1 - p2;
   
   double p4 = p1 * _par_lj[0] * ( 1.0 + _par_lj[1] * pow( dst, 2 ) ) + 1.0;
   
   _evdw += p3 / p4;
   
   /* electrostatic potential */
   
   double a1 = 4.0 - 3.0 * _par_el[0];
   double b1 = 2.0 * _par_el[0] - 3.0;
   
   double s1 = _par_el[1] * dst;
   
   double g1;
   
   if ( s1 < 1 )
    g1 = _par_el[0] + a1 * pow(s1, 2.0) + b1 * pow(s1, 3.0);
   else
    g1 = 1.0 / s1;
   
   double e2 = (*il1).getAtomCharge();
   
   double e3;
   
   if ( (*ip1).getPointType() == 0 )
    e3 = _complex_ele[(*ip1).getResidueCode()+30];
   else
    e3 = _complex_ele[(*ip1).getPointType()];
   
   _eele += e2 * e3 * g1;
   
   /* contact potential */
   
   _epmf += _complex_pmf[(*ip1).getPointType()][(*il1).getAtomType()][1] * 1.0 / ( 1.0 + exp( ( -0.5 * dst + 6.0 ) * ( dst - _complex_pmf[(*ip1).getPointType()][(*il1).getAtomType()][0] ) ) );
   
   /* hydrophobic potential */
   
   if ( (*ip1).getPointClass() == 2 || ( (*ip1).getPointClass() == 0 && (*ip1).getResidueCode() == 12 ) )
    if ( dst <= 9.0 )
     hpc1 += _complex_hpp[(*ip1).getResidueCode()] * ( 1.0 - 0.5 * ( 7.0 * pow( dst / 9.0, 2.0 ) - 9.0 * pow( dst / 9.0, 4.0 ) + 5.0 * pow( dst / 9.0, 6.0 ) - pow( dst / 9.0, 8.0 ) ) );
    
   /* hydrogen bond potential */
   
   if ( _complex_hdb[(*ip1).getPointType()][(*il1).getAtomType()][0] > 0.1 )
    _ehdb += ( -1.0 / ( _complex_hdb[(*ip1).getPointType()][(*il1).getAtomType()][1] * sqrt( 2.0 * PI ) ) ) * 
             exp( -0.5 * pow( ( dst - _complex_hdb[(*ip1).getPointType()][(*il1).getAtomType()][0] ) / _complex_hdb[(*ip1).getPointType()][(*il1).getAtomType()][1] , 2.0) );
   
   /* pocket-specific potential */
   
   if ( (*ip1).getPointClass() == 2 && dst <= _complex_pmf[(*ip1).getPointType()][(*il1).getAtomType()][0] )
   {
    map<pair<int,int>, double>::iterator is1;
    
    is1 = _complex_psp.find( make_pair( _protein_seq3[(*ip1).getResidueNumber()], (*il1).getAtomType() ) );
    
    if ( is1 != _complex_psp.end() )
     _epsp += (*is1).second;
   }
  }
  
  /* hydrophobic restraints */
  
  _ehpc += 0.5 * pow( ( hpc1 - _complex_hpl[(*il1).getAtomType()][0] ) / _complex_hpl[(*il1).getAtomType()][1], 2.0 ) - log( 1.0 / ( _complex_hpl[(*il1).getAtomType()][1] * sqrt( 2.0 * PI ) ) );
  
  /* kde potential */
  
  double kde1 = 0.0;
  int kde2 = 0;
  
  vector<CoordsKDE>::iterator ik1;
  
  for ( ik1 = _kde_xyz.begin(); ik1 < _kde_xyz.end(); ik1++ )
   if ( (*il1).getAtomType() == (*ik1).getPointType() )
   {
    kde1 +=  ( ( exp( -0.5 * pow( ( b_xyz[0] - (*ik1).getCoords(1) )/ _par_kde , 2.0 ) ) / ( _par_kde * sqrt( 2.0 * PI ) ) ) * 
               ( exp( -0.5 * pow( ( b_xyz[1] - (*ik1).getCoords(2) )/ _par_kde , 2.0 ) ) / ( _par_kde * sqrt( 2.0 * PI ) ) ) * 
               ( exp( -0.5 * pow( ( b_xyz[2] - (*ik1).getCoords(3) )/ _par_kde , 2.0 ) ) / ( _par_kde * sqrt( 2.0 * PI ) ) ) );
    
    kde2++;
   }
  
  if ( kde2 )
   _ekde += ( kde1 / ( (double) kde2 ) );
  
  /* position restraints */
  
  for ( ir1 = _ligand_mcs.begin(), ir2 = 0; ir1 != _ligand_mcs.end(); ir1++, ir2++ )
  {
   map<int,double>::iterator ir3x;
   map<int,double>::iterator ir3y;
   map<int,double>::iterator ir3z;
   
   ir3x = ((*ir1).xcoord).find( (*il1).getAtomNumber() + 1 );
   ir3y = ((*ir1).ycoord).find( (*il1).getAtomNumber() + 1 );
   ir3z = ((*ir1).zcoord).find( (*il1).getAtomNumber() + 1 );
   
   if ( ir3x != ((*ir1).xcoord).end() && ir3y != ((*ir1).ycoord).end() && ir3z != ((*ir1).zcoord).end() )
    rmsd1[ir2].push_back( pow( b_xyz[0] - (*ir3x).second, 2 ) + pow( b_xyz[1] - (*ir3y).second, 2 ) + pow( b_xyz[2] - (*ir3z).second, 2 ) );
  }
 }
 
 /* position restraints */
 
 list<double>::iterator ir5;
 
 for ( int ir6 = 0; ir6 < ir4; ir6++ )
 {
  double rmsd2 = 0.0;
  
  double rmsd3 = rmsd1[ir6].front();
  
  rmsd1[ir6].pop_front();
  
  for ( ir5 = rmsd1[ir6].begin(); ir5 != rmsd1[ir6].end(); ir5++ )
   rmsd2 += (*ir5);
  
  _elhm += ( rmsd3 * sqrt( rmsd2 / ( double ( rmsd1[ir6].size() ) ) ) );
 }
 
 if ( ir4 )
  _elhm = log( _elhm / (double) ir4 );
 
 for ( int il5 = 0; il5 < 3; il5++ )
  _ligand_center[il5] /= (double) _lna;
 
 _edst = sqrt( pow( _ligand_center[0] - _pocket_center[0], 2 ) +
               pow( _ligand_center[1] - _pocket_center[1], 2 ) +
               pow( _ligand_center[2] - _pocket_center[2], 2 )   );
 
 _evdw = _weights[0] * ( _evdw / (double) _lna );
 _eele = _weights[1] * ( _eele / (double) _lna );
 _epmf = _weights[2] * ( _epmf / (double) _lna );
 _ehpc = _weights[3] * ( _ehpc / (double) _lna );
 _ehdb = _weights[4] * ( _ehdb / (double) _lna );
 
 _edst = _weights[5] * _edst;
 
 _epsp = _weights[6] * ( _epsp / (double) _lna );
 _ekde = _weights[7] * ( _ekde / (double) _lna );
 
 _elhm = _weights[8] * _elhm;
 
 
 
 
 
 
 
 //cout << "E " << _evdw << " " << _eele << " " << _epmf << " " << _ehpc << " " << _ehdb << " " << _edst << " " << _epsp << " " << _ekde << " " << _elhm << endl;
 
 
 _etot = _evdw;
 
 
 
 
 if ( _etot < _ebst )
 {
  _ebst = _etot;
  
  for ( int bi = 0; bi < 3; bi++ )
  {
   _ligand_tra_best[bi] = _ligand_tra[bi];
   _ligand_rot_best[bi] = _ligand_rot[bi];
  }
 }
 
 _stats_evals++;
}

// ==================================================================================   getEnergy

double Complex::getEnergy( int ei )
{
 switch (ei)
 {
  case 1 : return _etot;
  case 2 : return _evdw;
  case 3 : return _eele;
  case 4 : return _epmf;
  case 5 : return _ehpc;
  case 6 : return _ehdb;
  case 7 : return _edst;
  
  default : return 0;
 }
}

// ==================================================================================   setTemperature

void Complex::setTemperature( double temp1 )
{
 _temp = temp1;
}

// ==================================================================================   getTemperature

double Complex::getTemperature( void )
{
 return _temp;
}

// ==================================================================================   setConfiguration

void Complex::setConfiguration( double conf1[] )
{
 for ( int ai = 0; ai < 3; ai++ )
 {
  _ligand_tra[ai] = conf1[ai];
  _ligand_rot[ai] = conf1[ai+3];
 }
}

// ==================================================================================   getFunctionEvals

int Complex::getFunctionEvals( void )
{
 return _stats_evals;
}

// ==================================================================================   setBest

void Complex::setBest( void )
{
 _etot = _ebst;
 
 for ( int bi = 0; bi < 3; bi++ )
 {
  _ligand_tra[bi] = _ligand_tra_best[bi];
  _ligand_rot[bi] = _ligand_rot_best[bi];
 }
}









    
    

































