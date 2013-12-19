

#include "decoys.h"
using namespace std;

//static double mccMatrix[5000][5000];


void Decoy(Complex * d_complex, std::string iname, std::string oname)
{
 //d_complex->clearContactMatrix();
 ifstream infile;
 fstream outfile;
/* double mccMatrix[1000][1000];
 for(int i=0; i<5000;i++){
  for(int j=0; j<5000; j++)
    mccMatrix[i][j]=0.0;
  }
*/
 string linecount;
 int number_of_lines=0;
 infile.open(iname.c_str());
 if(infile.is_open()){
        while(!infile.eof()){
            getline(infile,linecount);
            number_of_lines++;
        }
 }

 cout<<number_of_lines<<endl;
 outfile.open(oname.c_str());
 int line=0;
 int atoms, newpens, newlens;
 int conf=0;
 int conf2=0;
 infile.clear();
 infile.seekg(0, ios::beg);
 atoms=d_complex->getLigandAtomsTotal();
 int total_confs=((number_of_lines-1)/(atoms+3));
 cout<<total_confs<<endl;
 double newcoords[MAXLIG][3];
 
 while(line<number_of_lines){ 			//Establishes new "native" conformation 
 infile>>newpens;
 line++;
 infile>>newlens;
 line++;
 cout<<newpens<<" "<<newlens<<endl;
 d_complex->setProteinEnsembleCurrent(newpens);
 d_complex->setLigandEnsembleCurrent(newlens);

 for(int i=0;i<atoms;i++ ){
  infile>>newcoords[i][0]>>newcoords[i][1]>>newcoords[i][2];
  line++;
 }
  
  d_complex->setConfCoords(newcoords);
  d_complex->clearContactMatrix();
  d_complex->createContactMatrix();
  d_complex->calculateEnergy();
  outfile<<d_complex->getEnergy(1)<<" ";
  //mccMatrix[conf][conf2]=d_complex->getEnergy(1);

  while(!infile.eof()){ //Calculate mcc for each conformation relative the the current "native"
  conf2++;
  infile.clear();
  infile.seekg(0, ios::beg);
  //cout<<"THE LINE IS: "<<line<<" "<<"The matrix entry is: ";
  line++;
  GotoLine(infile,line);
  infile>>newpens;
  //cout<<"newpens is: "<<newpens<<endl;
  line++;
  infile>>newlens;
  //cout<<"newlens is: "<<newlens<<endl;
  line++;

  d_complex->setProteinEnsembleCurrent(newpens);
  d_complex->setLigandEnsembleCurrent(newlens);

  for(int j=0;j<atoms;j++){
   infile>>newcoords[j][0]>>newcoords[j][1]>>newcoords[j][2];
   line++;
  }
  d_complex->setConfCoords(newcoords);
  d_complex->calculateEnergy();
  outfile<<d_complex->getEnergy(1)<<" ";
  if(conf2==total_confs)
   outfile<<endl;
  /* THIS IS GOOD CODE, TRYING SOMETHING ELSE
  mccMatrix[conf][conf2]=d_complex->getEnergy(1);
  mccMatrix[conf2][conf]=d_complex->getEnergy(1);
  cout<<mccMatrix[conf][conf2]<<endl;
  */
  }
 
 conf++;
 conf2=conf;
 line=(atoms+4)*conf-(conf-1); //Changed 5 to 4 because I changed the input file
 infile.clear();
 infile.seekg(0, ios::beg);
 GotoLine(infile,line); 
 cout<<conf<<endl; 
 }

/*
 for(int i=0;i<conf+1;i++){
  for(int j=0;j<conf+1;j++)
   outfile<<mccMatrix[i][j]<<" ";
  outfile<<endl;
 }
*/  
 infile.close();
 outfile.close();
 
}
std::ifstream& GotoLine(std::ifstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

