#ifndef  GPU_CUH
#define  GPU_CUH


__global__ void InitCurand_d ();

__global__ void ExchangeReplicas_d (const int, const int);

__global__ void MonteCarlo_Init_d (const int, const int);

__global__ void MonteCarlo_d (const int, const int, const int, const int);







__device__ void Move_d (const int, Ligand * __restrict__); //, LigCoord * __restrict__);

__device__ void CalcEnergy_d (const int, Ligand * __restrict__, const Protein *);

__forceinline__ __device__ void Accept_d (const int, Ligand * __restrict__, const float, const int);




// util_d.cu

__device__ void InitAcs_d (const int);

__device__ void InitLigRecord_d (const int, const int, const int);

//__forceinline__ __device__ void BackupLigCoord_d (const int, Ligand *);

__device__ void ComputeMoveMatrix_d (const int, const int, Ligand *);

__device__ void RecordLigand_d (const int, const int, const int, const int, const int, const Ligand *);



__forceinline__ __device__ float MyRand_d ();

//__forceinline__ __device__ int Minimal_int_d (const int, const int);

__forceinline__ __device__ void SumReduction1D_d (const int, float *);

__forceinline__ __device__ void SumReduction1D_5_d (const int, float *, float *, float *, float *, float *);

__forceinline__ __device__ void SumReduction2D_d (float a[BDy][BDx]);

__forceinline__ __device__ void SumReduction2D_2_d (float a[BDy][BDx], float b[BDy][BDx]);




#endif
