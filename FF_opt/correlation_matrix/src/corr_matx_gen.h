#ifndef CORR_MATX_GEN_H
#define CORR_MATX_GEN_H

void
CorrMatxGen (float * corr_mat, 
	     const float * track_mat, 
	     const int total_rows, 
	     const Ligand * lig, 
	     const Protein * prt,
	     const EnePara * enepara,
	     const ComplexSize complexsize);

#endif // CORR_MATX_GEN_H
