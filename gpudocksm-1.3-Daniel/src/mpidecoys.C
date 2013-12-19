

//"include <mpi.h>
#include "mpidecoys.h"

const int SIZE=70000;

//float static mcc_matrix[SIZE][SIZE];

using namespace std;


void Decoy(Complex * d_complex, std::string iname, std::string oname){
 ifstream infile;
 fstream outfile;


 vector<string> data;
 vector<float> mcc;

 string inputline;
 int number_of_lines=0;
 infile.open(iname.c_str());
 if(infile.is_open()){
        while(!infile.eof()){
            getline(infile,inputline);
   	    data.push_back(inputline);
            number_of_lines++;
        }
 }


 cout<<number_of_lines<<endl;
 outfile.open(oname.c_str());
 int line=0;
 int atoms, newpens, newlens;
 int conf=0;
 infile.clear();
 infile.seekg(0, ios::beg);
 atoms=d_complex->getLigandAtomsTotal();


 int total_confs=((number_of_lines-1)/(atoms+3));
 total_confs=11; //included for testing
 int mccMatrix_size=0;
 for (int i=total_confs;i>0;i--)
   mccMatrix_size+=i;

 float *mccMatrix = new float[mccMatrix_size];

 cout<<total_confs<<endl;
 double newcoords[MAXLIG][3];
 
 int n1; // index for fetching new native parameters

 while(conf<total_confs){ 			//Establishes new "native" conformation 

   n1=((atoms+3)*conf);

   newpens=atoi((data[n1]).c_str());
   newlens=atoi((data[n1+1]).c_str());

   d_complex->setProteinEnsembleCurrent(newpens);
   d_complex->setLigandEnsembleCurrent(newlens);


   for(int j=n1+2, k=0;k<atoms;j++, k++){
      stringstream vec1_entry((data[j]).c_str());
      vec1_entry>>newcoords[k][0]>>newcoords[k][1]>>newcoords[k][2];
   }

   d_complex->setConfCoords(newcoords);
   d_complex->clearContactMatrix();
   d_complex->createContactMatrix();

   int rank;
   int size;
  // MPI_Init();

  // MPI_Comm_size(MPI_COMM_WORLD,&size);
  // MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int pens_p;
    int lens_p;
    Complex complex_p = *d_complex;
    double coords_p[MAXLIG][3];
  
   

   int conf2_p;
   int n; //used for vector entries


   for(conf2_p=conf;conf2_p<total_confs;conf2_p++){
   //for (int cycle = 0; cycle <= total_confs / size + 1; ++cycle) {
	//int conf2_p = size * cycle + rank;
	if (conf2_p <= total_confs) {

	//Calculate mcc for each conformation relative the the current "native"

		n=((atoms+3)*conf2_p);

		pens_p=atoi((data[n]).c_str());
		lens_p=atoi((data[n+1]).c_str());
		
		complex_p.setProteinEnsembleCurrent(pens_p);
		complex_p.setLigandEnsembleCurrent(lens_p);

		for(int j=n+2, k=0;k<atoms;j++, k++){
			stringstream vec_entry((data[j]).c_str());
			vec_entry>>coords_p[k][0]>>coords_p[k][1]>>coords_p[k][2];
		}

		complex_p.setConfCoords(coords_p);
		complex_p.calculateEnergy();
		//mcc.push_back(complex_p.getEnergy(1));
		
		int mcc_index=0;

		if(conf==0){
		   mcc_index=conf2_p; //Does rank start at 0 or 1?
		}
		else{
		   for(int i=1;i<=conf;i++){
		     mcc_index+=total_confs-(i-1);
		   }
		   mcc_index+=conf2_p-conf;
		}
		mccMatrix[mcc_index]=complex_p.getEnergy(1);
	}
   }

  // MPI_Finalize();



   conf++;
   cout<<conf<<endl; 
 }
//================================================================ Print the entire matrix from 1D Array

//vector<float>::iterator l1;
 int l1,l2,l3,l5;

 for(l1 = 0, l2=1, l3=1; l1 < mccMatrix_size; l1++, l2++){
   
   if (l2==l3){
     l5=0;
     for(int l4=0; l4<(l3-1); l4++){
	if(l4>=2)
          l5+=l4-1;
      outfile<<mccMatrix[(((l4)*total_confs)-(l5))+(l2-l4)-1]<<" ";
     }
   }
  
   outfile<<mccMatrix[l1]<<" ";
   if(l2==total_confs){
   	outfile<<endl;
	l3++;
	l2=l3-1;
   }
 } 

 delete [] mccMatrix;

 infile.close();
 outfile.close();

}
std::ifstream& GotoLine(std::ifstream& file, int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

