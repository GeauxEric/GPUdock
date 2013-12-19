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


#ifndef __COORDS_H_
#define __COORDS_H_

#include<string>
#include<deque>

#include "size.h"
#include "data.h"

using namespace std;

// ==================================================================================   CoordsProtein

class CoordsProtein {
        
  private:
    
    int    _r;         // residue number
    int    _n;         // effective point number
    double _x[MAXEN1]; // x coords (ensemble)
    double _y[MAXEN1]; // y coords (ensemble)
    double _z[MAXEN1]; // z coords (ensemble)
    int    _t;         // effective point type
    int    _c;         // effective point class
    int    _d;         // residue code

  public:

    CoordsProtein( int, int, int, int, int );
    
    CoordsProtein( void );
    
    ~CoordsProtein( void );
    
    int getPointNumber( void );
    
    int getResidueNumber( void );
    
    int getResidueCode( void );
    
    double getCoords( int, int );
    
    void setCoords( double, double, double, int );
    
    int getPointType( void );
    
    int getPointClass( void );
};

// ==================================================================================   CoordsLigand

class CoordsLigand {
        
  private:
    
    int    _n;         // atom number
    double _x[MAXEN2]; // x coords (ensemble)
    double _y[MAXEN2]; // y coords (ensemble)
    double _z[MAXEN2]; // z coords (ensemble)
    string _a;         // atom name
    int    _t;         // atom type
    double _c;         // atom charge
    
  public:

    CoordsLigand( int, string, int, double );
    
    CoordsLigand( void );
    
    ~CoordsLigand( void );
    
    int getAtomNumber( void );
    
    double getCoords( int, int );
    
    void setCoords( double, double, double, int );
    
    string getAtomName( void );
    
    int getAtomType( void );
    
    void setAtomType( int );
    
    double getAtomCharge( void );
    
    void setAtomCharge( double );
};

// ==================================================================================   CoordsKDE

class CoordsKDE {
        
  private:
    
    int    _n; // point number
    double _x; // x coords (ensemble)
    double _y; // y coords (ensemble)
    double _z; // z coords (ensemble)
    int    _t; // atom type

  public:

    CoordsKDE( int, int, double, double, double );
    
    CoordsKDE( void );

    ~CoordsKDE( void );
    
    int getPointNumber( void );
    
    double getCoords( int );
    
    void setCoords( double, double, double );
    
    int getPointType( void );
};

#endif
