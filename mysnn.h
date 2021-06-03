#ifndef _MY_SNN_H_
#define _MY_SNN_H_
#include <stdint.h>
#include <stdlib.h> 
#include <stdio.h>
#include <athread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include "mpi.h"
#include "mypbuf.h"
typedef uint8_t bool;
#include "swstruct.h"
#define NTh 64

extern int timestep;
static inline unsigned long rpcc()
{
	unsigned long time;
	asm("rtc %0": "=r" (time) : );
	return time;
}

typedef struct PoissonRate_s
{
	int numPois;
	float *h_rates_;

}PoissonRate_t;
typedef struct snnInfo_s
{
	unsigned int	simTimeMs;
	uint64_t        simTimeSec;
	unsigned int	simTime;
	uint8_t         maxDelay;
	unsigned int    numGrp,numConn;	
	int	        	numN,numNReg,numNExcReg,numNInhReg,numNPois;
	unsigned int	maxSpikesD1;
	unsigned int	maxSpikesD2;
	unsigned int    simTimeRunStart;
	unsigned int    simTimeRunStop;

	double dAMPA,dNMDA,rNMDA,sNMDA;
	double dGABAa,dGABAb,rGABAb,sGABAb;

	unsigned int		*timeTableD1;
	unsigned int		*timeTableD2;
	unsigned int		*firingTableD1;
	unsigned int		*firingTableD2;
	unsigned int    MaxFiringRate;
	unsigned int	secD1fireCntHost;
	unsigned int	secD2fireCntHost;
	unsigned int	spikeCountAll1secHost;
	Pbuf_t 		*pbuf;
	bool sim_with_homeostasis;
	bool sim_with_conductances;
	bool sim_with_NMDA_rise;
	bool sim_with_GABAb_rise;
	int         	*nSpikeCnt;
	uint32_t    	*lastSpikeTime;

	/*****for SW******/
	swInfo_t swInfo[NTh];
	int preN,MaxN;
	synInfo_t* sInfoHost;
	neurInfo_t* nInfoHost;
	int Nsyn;
	int Ndelay;
	float dt;
	int Ndt, Nop;
	int randSeed;
	//mpi firingTable
	spikeTime_t *firingTableHost;
	spikeTime_t *firingTableAll;
	int NNgroup,NN,ND;
	int gStart,gSize;
	int *displs, *recvCount;
}snnInfo_t;

typedef struct grpInfo_s{
	int StartN;
	int EndN;
	int SizeN;
	int8_t		MaxDelay;
	int 		FiringCount1sec;

	bool	isSpikeGenerator;
	int			NewTimeSlice;
	int			CurrTimeSlice;
	uint32_t 	SliceUpdateTime;
	float   	RefractPeriod;
	PoissonRate_t*	RatePtr;

	bool	WithSTP;
	bool 		WithHomeostasis;
	unsigned int	Type;

	unsigned int numSynPre,numSynPost;

	float Izh_a,Izh_a_sd;
	float Izh_b,Izh_b_sd;
	float Izh_c,Izh_c_sd;
	float Izh_d,Izh_d_sd;

}grpInfo_t;

typedef struct connInfo_s{
	uint8_t connId;
	uint8_t maxDelay;
	int grpSrc,grpDest;
	float initWt,maxWt;
	float mulSynFast,mulSynSlow;
	int maxSynPost,maxSynPre;
	int numSyn;
	uint32_t connProp;//bit info & op
	float p;//conn probability
}connInfo_t;

#if 1
int runNetwork(snnInfo_t *snnInfo, grpInfo_t *grpInfo,
               connInfo_t *connInfo, neurInfo_t *neurInfo,
               synInfo_t *synInfo_t, swInfo_t *swInfo,
               int _nmsec, bool printRun);

void initNetwork(snnInfo_t *sInfo, grpInfo_t *gInfo,connInfo_t *cInfo,
                int numGrp,int numConn,int randSeed);
void setNeuronParameters(grpInfo_t *gInfo,int gid,float izh_a,
		float izh_a_sd,float izh_b,float izh_b_sd,float izh_c,
		float izh_c_sd,float izh_d,float izh_d_sd);
void setConductances(snnInfo_t *sInfo,bool isSet,int tdAMPA,int trNMDA,
		int tdNMDA,int tdGABAa,int trGABAb,int tdGABAb);
void createNetwork(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo);
void setupNetwork(snnInfo_t *sInfo,grpInfo_t *gInfo);
void setRates(PoissonRate_t *ratePtr,int numPois,float rate);
void setSpikeRate(grpInfo_t *gInfo,int gid,PoissonRate_t* rPtr,int refPeriod);
void freeNetwork();
#endif

#endif


