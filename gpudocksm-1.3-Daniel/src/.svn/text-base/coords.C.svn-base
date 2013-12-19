//Daniel Case 11-20-12
//Coordinates


#include "coords.h"
using namespace std;

// =========================================================	CoordsProtein
CoordsProtein::CoordsProtein( int ar, int an, int at, int ad, int ac )
{
 _r = ar;
 _n = an;
 _t = at;
 _d = ad;
 _c = ac;

 for ( int ai = 0; ai < MAXEN1; ai++ )
 {
  _x[ai] = 0.0;
  _y[ai] = 0.0;
  _z[ai] = 0.0;
 }
}

CoordsProtein::CoordsProtein( void )
{
 _r = 0;
 _n = 0;
 _t = 0;
 _d = 0;
 _c = 0;

 for ( int ai = 0; ai < MAXEN1; ai++ )
 {
  _x[ai] = 0.0;
  _y[ai] = 0.0;
  _z[ai] = 0.0;
 }
}

CoordsProtein::~CoordsProtein() {}

int CoordsProtein::getPointNumber( void )
{
 return _n;
}

int CoordsProtein::getResidueNumber( void )
{
 return _r;
}

int CoordsProtein::getResidueCode( void )
{
 return _d;
}
double CoordsProtein::getCoords( int an, int ai )
{
 switch ( an )
 {
  case  1  :  return _x[ai];
  case  2  :  return _y[ai];
  case  3  :  return _z[ai];

  default  :  return 0;
 }
}

void CoordsProtein::setCoords( double ax, double ay, double az, int ai)
{
 _x[ai] = ax;
 _y[ai] = ay;
 _z[ai] = az;
}

int CoordsProtein::getPointType( void )
{
 return _t;
}

int CoordsProtein::getPointClass( void )
{
 return _c;
}

// ==================================================================================   CoordsLigand

CoordsLigand::CoordsLigand( int an, string aa, int at, double ac )
{
 _n = an;
 _a = aa;
 _t = at;
 _c = ac;

 for ( int ai = 0; ai < MAXEN2; ai++ )
 {
  _x[ai] = 0.0;
  _y[ai] = 0.0;
  _z[ai] = 0.0;
 }
}

CoordsLigand::CoordsLigand( void )
{
 _n = 0;
 _a = "";
 _t = 0;
 _c = 0.0;

 for ( int ai = 0; ai < MAXEN2; ai++ )
 {
  _x[ai] = 0.0;
  _y[ai] = 0.0;
  _z[ai] = 0.0;
 }
}

CoordsLigand::~CoordsLigand() {}

int CoordsLigand::getAtomNumber( void )
{
 return _n;
}

double CoordsLigand::getCoords( int an, int ai )
{
 switch ( an )
 {
  case  1  :  return _x[ai];
  case  2  :  return _y[ai];
  case  3  :  return _z[ai];

  default  :  return 0;
 }
}

void CoordsLigand::setCoords( double ax, double ay, double az, int ai)
{
 _x[ai] = ax;
 _y[ai] = ay;
 _z[ai] = az;
}

string CoordsLigand::getAtomName( void )
{
 return _a;
}

int CoordsLigand::getAtomType( void )
{
 return _t;
}

void CoordsLigand::setAtomType( int at )
{
 _t = at;
}

double CoordsLigand::getAtomCharge( void )
{
 return _c;
}

void CoordsLigand::setAtomCharge( double ac )
{
 _c = ac;
}


