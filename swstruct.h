#ifndef _SWSTRUCT_H_
#define _SWSTRUCT_H_

#include <assert.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

inline int isnan(double x) { return x != x; }
inline int isinf(double x) { return !isnan(x) && isnan(x - x); }
inline int iserr(double x) {return isnan(x)||isinf(x);}


typedef struct spikeTime_s{
	unsigned short time;//unit: resolution interval
	unsigned short nid;
}spikeTime_t;

typedef struct synInfo_s{
	unsigned short postId;//post neuron id
	unsigned short dl;//unit: resolution interval
	float wt;
}synInfo_t;

typedef struct neurInfo_s{

	float Izh_b,Izh_a,Izh_c,Izh_d;
	
	float voltage;
	float recovery;

	float gAMPA;
	float gNMDA_d;
	float gGABAa;
	float gGABAb_d;
	int nSpikeCnt;
	
}neurInfo_t;

typedef struct swInfo_s{
	int StartN;
	int SizeN;

	double dAMPA,dNMDA,rNMDA,sNMDA;
	double dGABAa,dGABAb,rGABAb,sGABAb;

	int sim_with_conductances;
	int WithSTDP;
	int WithSTP;
	int WithHomeostasis;
	int sim_with_homeostasis;
	int sim_with_NMDA_rise;
	int sim_with_GABAb_rise;
	int MaxDelay;
	
	//connectInfo
	int preN; //num of pre neurons (3rd dim)
	int Ndelay;
	int MaxN; //max num of connections in a Th for a preN (1st dim)
	int Ndma;
	int Nop,Ndt;
	float dt;
	int nSpikePoisAll;
	int fireCnt;
	synInfo_t *sInfoHost;
	neurInfo_t *nInfoHost;
	
	//mpi firingTable
	spikeTime_t *firingTableHost;
	spikeTime_t *firingTableAll;
	int NNgroup,NN,ND;
	int rank,nproc;
	int gStart,gSize;
	float currentfactor;
}swInfo_t;

#endif
