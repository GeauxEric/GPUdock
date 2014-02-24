#ifndef  DOCK_H
#define  DOCK_H


#include <string>
#include "size.h"
#include "size_gpu.cuh"


struct LigandFile
{
  std::string path; // lig file path
  std::string id;	// ligand name
  std::string molid; // MOLID;

  int conf_total;	// number of conformation for one lig
  int raw_conf;		// raw total conf in the .sdf without small rmsd excluded

  int lna;	// total effective pts number
  int lnb;	// total bonds number

  //int ensemble_total;	// ENSEMBLE_TOTAL in the .sdf file
};


struct ProteinFile
{
  std::string path;  // prt file path
  std::string id;    // prt name

  int conf_total;	// number of conformations for one lig
  int pnp;	// total effective pts number
  int pnr;	// total bonds number
};


struct LhmFile
{
  std::string path;
  std::string ligand_id;

  int pos; // number of mcs positions
};


struct EneParaFile
{
  std::string path;
};

struct NorParaFile
{
  std::string path_a;
  std::string path_b;
};

struct WeightFile
{
  std::string path;
};

struct InputFiles
{
  LigandFile lig_file;
  ProteinFile prt_file;
  LhmFile lhm_file;
  EneParaFile enepara_file;
  WeightFile weight_file;
  NorParaFile norpara_file;
};






struct ComplexSize
{
  // replica numbers
  int n_lig; // number of ligand conf
  int n_prt; // number of protein conf
  int n_tmp; // number of temperature
  int n_rep; // n_rep = n_lig * n_prt * n_tmp;

  // residue numbers (of per replica)
  int lna; // number of ligand points
  int pnp; // number of protein points
  int pnk; // number of kde points
  int pos; // number of mcs positions
};







struct Energy
{
  float e[MAXWEI];
  /*
     0 - vdw
     1 - ele
     2 - pmf
     3 - psp
     4 - hdb

     5 - hpc
     6 - kde
     7 - lhm

     8 - dst

     9 - total
   */
  float cmcc;
};



struct TmpEnergy
{
  float vdw[TperB];
  float ele[TperB];
  float pmf[TperB];
  float psp[TperB];
  float hdb[TperB];

  float hpc_val[TperB];

  float kde_val[TperB];
  float kde_sz[TperB];

  float rmsd_val[TperB];
  float rmsd_sz[TperB];
};


struct ConfusionMatrix
{
  int matrix[MAXLIG][MAXPRO];
  int lig_conf, prt_conf;
  int lig_sz, prt_sz;
};


struct LigCoord
{
  float x[MAXLIG];		// Ligand x coords
  float y[MAXLIG];		// Ligand y coords
  float z[MAXLIG];		// Ligand z coords
  float center[3];		// ligand geometric center
};



struct Ligand0
{
  LigCoord coord_orig;           //                                      used

  int t[MAXLIG];		// atom type                            used
  float c[MAXLIG];		// atom charge                          used
  int n[MAXLIG];		// atom number                          used

  int lna;			// number of ligand atoms               used
  int lnb;			// number of ligand bonds               NOT USED

  float pocket_center[3];	// pocket center                        used
                                // should belong to "Protein" structure


  float lens_rmsd;		// ensemble rmsd                        NOT USED

  float mw;			// ligand molecular weight              NOT USED
  float logp;			// ligand water/octanol partition coeff NOT USED
  float psa;			// ligand polar surface area            NOT USED
  float mr;			// ligand molar refractivity            NOT USED

  int hbd;			// ligand hydrogen bond donors          NOT USED
  int hba;			// ligand hydrogen bond acceptors       NOT USED

  std::string a[MAXLIG];	// atom name                            NOT USED
  std::string id;		// ligand id                            NOT USED
  std::string smiles;		// ligand smiles                        NOT USED
};





struct Ligand
{
  // coord_center is always under lab system
  // coord_xyz for "orig" is under ligand_ceter system
  // coord_xyz for "new" is under lab system
  LigCoord coord_orig;
  LigCoord coord_new;

  // translation x y z, rotation x y z
  float movematrix_old[6];       // old matrix
  float movematrix_new[6];       // trail matrix

  Energy energy_old;		//                                      used
  Energy energy_new;		//                                      used

  int t[MAXLIG];		// atom type                            used
  float c[MAXLIG];		// atom charge                          used
  int n[MAXLIG];		// atom number                          used
                                // n == index???

  int lna;			// number of ligand atoms               used
  
  // confusion matrix
  int native_confusion_matx[MAXLIG][MAXPRO];
  int decoy_confusion_matx[MAXLIG][MAXPRO];
};











struct Protein0
{
  float x[MAXPRO];		// residue x coord                      used
  float y[MAXPRO];		// residue y coord                      used
  float z[MAXPRO];		// residue z coord                      used
  int n[MAXPRO];		// effective point number               NOT USED
  int t[MAXPRO];		// effective point type                 used
  int c[MAXPRO];		// effective point class                used
  int d[MAXPRO];		// redidue code                         used

  int pnp;			// number of protein effective points   used
  int pnr;			// number of protein residues           NOT USED

  int r[MAXPRO];		// residue number                       replaced
  int seq3[MAXPRO];		// aa sequence numbering                replaced


  //std::string protein_seq1;  // aa sequence
  //char protein_seq2[MAXPRO]; // aa sequence
};






struct Protein
{
  float x[MAXPRO];		// residue x coord
  float y[MAXPRO];		// residue y coord
  float z[MAXPRO];		// residue z coord

  int t[MAXPRO];		// effective point type
  int c[MAXPRO];		// effective point class


  float ele[MAXPRO];            // dt = prt->t[i] == 0 ? prt->d[i] + 30 : prt->t[i];
                                // enepara->ele[dt]

  int seq3r[MAXPRO];            // prt->seq3r[i] == prt->seq3[prt->r[i]];
  int c0_and_d12_or_c2[MAXPRO]; // (prt_c == 0 && prt_d == 12)) || (prt_c == 2)
  float hpp[MAXPRO];            // enepara->hpp[prt->d[i]]

  int pnp;			// number of protein effective points

  float pocket_center[3];
};






struct Psp0
{
  float psp[MAXPRO][MAXLIG];                                           // replaced
  int n;			// total number of PSP point           NOT USED
};



struct Psp
{
  float psp[MAXLIG][MAXPRO];                                           // replaced

  float sparse_psp[MAXLIG][MAXPRO];                                    // replaced

};






struct Kde0
{
  float x[MAXKDE];		// KDE x coord                          used
  float y[MAXKDE];		// KDE y coord                          used
  float z[MAXKDE];		// KDE z coord                          used
  int n[MAXKDE];		// KDE point number                     NOT USED
  int t[MAXKDE];		// KDE atom type                        used

  int pnk;			// number of kde points                 used
  int pns[MAXTP2];		// number of specific kde points        NOT USED
};




struct Kde
{
  float x[MAXKDE];		// KDE x coord                          used
  float y[MAXKDE];		// KDE y coord                          used
  float z[MAXKDE];		// KDE z coord                          used
  int t[MAXKDE];		// KDE atom type                        used

  int pnk;			// number of kde points                 used
};








struct Mcs0
{
  float x[MAXMCS];              //                                      used
  float y[MAXMCS];              //                                      used
  float z[MAXMCS];              //                                      used

  int total;                    //                                      NOT USED
  float tcc;                    //                                      used
};




struct Mcs
{
  float x[MAXMCS];              //                         used, mcs->y[lig_n]
  float y[MAXMCS];              //                         used
  float z[MAXMCS];              //                         used

  float tcc;                    //                         used
};








struct EnePara0
{
  float vdw[MAXTP1][MAXTP2][2];	// L-J potential                        vdw[prt_t][lig_t][]
  float ele[MAXTP3];		// electrostatic potential              ele[prt_d + 30], ele[prt_t]
  float pmf[MAXTP1][MAXTP2][2];	// contact potential                    pmf[prt_t][lig_t][]
  float hpp[MAXTP4];		// protein hydrophobicity               hpp[prt_d]
  float hpl[MAXTP2][2];		// ligand hydrophobicity                hpl[lig_t][]
  float hdb[MAXTP1][MAXTP2][2];	// ligand hydrophobicity                hdb[prt_t][lig_t][]

  float lj[3];			// L-J params
  float el[2];			// electrostatic params
  float kde;			// kde bandwidth

  float w[MAXWEI];		// weights for energy terms
  float a_para[MAXWEI];         // the a parameter in normalization
  float b_para[MAXWEI];         // the b parameter in normalization
};





// optimized EnePara
struct EnePara
{
  // L-J
  float p1a[MAXTP2][MAXTP1];
  float p2a[MAXTP2][MAXTP1];
  float lj0, lj1;

  // electrostatic
  float el0;
  float el1;
  float ele[MAXTP3];
  float a1; // 4.0f - 3.0f * el0;
  float b1; // 2.0f * el0 - 3.0f;

  // contact
  float pmf0[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float pmf1[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float hdb0[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float hdb1[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float hdb2[MAXTP2][MAXTP1];  // interchange to [lig][prt]

  // hydrophobic
  float hpp[MAXTP4];
  float hpl0[MAXTP2];
  float hpl1[MAXTP2];
  float hpl2[MAXTP2];

  // kde
  float kde2; // -0.5f / (kde * kde)
  float kde3; // powf (kde * sqrtf (2.0f * PI), 3.0f)

  // weights for energy terms
  float w[MAXWEI];
  float a_para[MAXWEI];         // the a parameter in normalization
  float b_para[MAXWEI];         // the b parameter in normalization
};







struct Temp
{
  float t;
  float minus_beta;
  int order;
};





// replica[n_rep]
// replica[n_prt][n_tmp][n_lig]

struct Replica
{
  int idx_rep; // n_rep, replica

  int idx_prt; // n_prt, protein
  int idx_tmp; // n_tmp, temperature
  int idx_lig; // n_lig, ligand
};


struct ExchgPara{
  float floor_temp;   // lowest temperature in all replicas
  float ceiling_temp;   // highest temperature in all replicas
  int num_temp;      // number of temperatures for the same ligand and protein conformations
};


struct McPara
{
  int steps_total;
  int steps_per_dump;
  int steps_per_exchange;

  float move_scale[6]; // translation x y z, rotation x y z

  char outputdir[MAXSTRINGLENG];
  char outputfile[MAXSTRINGLENG];
};

struct McLog
{
  double t0, t1, t2; // time escape
  int ac_mc;  // MC acceptance counter
  int acs_mc[MAXREP];   // MC acceptance counter for all replicas
  int ac_temp_exchg;
  int acs_temp_exchg[MAXREP]; 
  
  // int ac_lig_exchg;
  // int acs_lig_exchg[MAXREP];
};




struct LigRecordSingleStep
{
  Replica replica;
  Energy energy;
  float movematrix[6]; // // translation x y z, rotation x y z
  int step;
};




struct LigRecord
{
  LigRecordSingleStep step[STEPS_PER_DUMP];
};

struct LigMoveVector
{
  float ele[6]; // translation xyz, rotation xyz
};




#endif // DOCK_H

