#include "mysnn.h"

#define true 1
#define false 0
#define PROPAGATED_BUFFER_SIZE  (1023)
#define ALL -1
#define MAX_SIMULATION_TIME     ((uint32_t)(0x7fffffff))

#define UNKNOWN_NEURON  (0)
#define POISSON_NEURON  (1 << 0)
#define TARGET_AMPA     (1 << 1)
#define TARGET_NMDA     (1 << 2)
#define TARGET_GABAa    (1 << 3)
#define TARGET_GABAb    (1 << 4)

volatile unsigned int dma[64],spike[64],NS_group,NSall,numSpike[64],nspikeall;

extern SLAVE_FUN(initSW)(swInfo_t*);
extern SLAVE_FUN(freeSW)(void*);
extern SLAVE_FUN(SnnSim)(void*);
extern SLAVE_FUN(StateUpdate)(void*);
extern SLAVE_FUN(SpikeDeliver)(void*);

static bool updateTime(snnInfo_t *snnInfo);
static void doSnnSim(snnInfo_t*,grpInfo_t*,connInfo_t*,neurInfo_t*,synInfo_t*, swInfo_t*);


double DO_T1,DO_T2,DO_T3,DO_T4,DO_T5,DO_T6,DO_T7,DO_T8,DO_TA;
long CC0,CC1;

int runNetwork(snnInfo_t *sInfo, grpInfo_t *gInfo,
	       connInfo_t *cInfo, neurInfo_t *nInfo,
	       synInfo_t *synInfo, swInfo_t *swInfo,
	       int _nmsec, bool printRun) 
{

	struct timeval t0,t4,t5,t6,t7,t8,t9;
	struct timeval *time1,*time2;
	gettimeofday( &t0, NULL );

	int rank, nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	assert(_nmsec >= 0);
	int runDurMs = _nmsec;

	sInfo->simTimeRunStart = sInfo->simTime;
	sInfo->simTimeRunStop  = sInfo->simTime+runDurMs;

	DO_T5=0.0; DO_T6=0.0; DO_T7=0.0; DO_T8=0.0; DO_TA=0.0;CC0=0;CC1=0;

	athread_init();//initialize athread!!!!
	athread_spawn(initSW,swInfo);
    athread_join();
	int i;
	nspikeall=0;
	for(i=0; i<runDurMs; i++) 
	{
		doSnnSim(sInfo,gInfo,cInfo,nInfo,synInfo,swInfo);
		updateTime(sInfo);
	}
	athread_spawn(freeSW,swInfo);
    athread_join();
	athread_halt();

	unsigned long cdma=0,cspike=0;
	for(i=0;i<63;i++)
	{
		cdma+=dma[i];
		cspike+=spike[i];
	}
	cdma/=63; cspike/=63;
	int nSpikeAll=0;
	for(i=0;i<64;i++)
	{
		nSpikeAll+=numSpike[i];
	}
	//printf("rank=%d nspikeall=%d nSpikeAll=%d cdma=%ld,cspike=%ld\n",rank,nspikeall,nSpikeAll,cdma,cspike);
	
	if (rank == 0)
	{
		printf("Number of pulses %d\n", nspikeall);
		printf("doSnnSim::T5: %f\n", DO_T5);
		printf("doSnnSim::T6: %f\n", DO_T6);
		printf("doSnnSim::T7: %f\n", DO_T7);
	}

	return 0;
}

static void doSnnSim(snnInfo_t *sInfo, grpInfo_t *gInfo,
               connInfo_t *cInfo, neurInfo_t *nInfo,
               synInfo_t *synInfo, swInfo_t *swInfo) 
{
	struct timeval t0,t1,t2,t3,t4,t5,t6,t7,t8,t9;
	struct timeval *tm1,*tm2;
	gettimeofday( &t0, NULL );

	long st,st0,ed;
	{{gettimeofday( &t4, NULL );st = rpcc();}
	{athread_spawn(StateUpdate,swInfo);//
	athread_join();}
	{gettimeofday( &t5, NULL );st0= rpcc();}}
	sInfo->simTime++;

	int rank, nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int i;
	int iND = sInfo->simTime%sInfo->Ndelay;
#if 1
	spikeTime_t *sendBuf, *recvBuf;
	int *displs, *recvCount;
	int root = 0;
	int senddatanum;

	sendBuf = &(sInfo->firingTableHost[0]);
	senddatanum= NS_group;

	displs = sInfo->displs; // displs
	recvCount = sInfo->recvCount; // send num

	MPI_Gather(&senddatanum,1,MPI_INT,recvCount,1, MPI_INT,root,MPI_COMM_WORLD);

	if(!rank)
	{
		displs[0] = 0;
		NSall = recvCount[0];
		for(i=1;i<nproc;i++)
		{
			displs[i] = displs[i-1]+recvCount[i-1];
			NSall += recvCount[i];
		}
	}

	recvBuf = &(sInfo->firingTableAll[iND*sInfo->NN]);

	MPI_Gatherv(sendBuf, 
		senddatanum, 
		MPI_INT, 
		recvBuf, 
		recvCount, 
		displs, 
		MPI_INT,
		root, 
		MPI_COMM_WORLD);
#if 1	
	MPI_Bcast(&NSall,1,MPI_INT,root,MPI_COMM_WORLD);
	MPI_Bcast(recvBuf,NSall,MPI_INT,root,MPI_COMM_WORLD);
	nspikeall += NSall;
#endif 

#endif

	gettimeofday( &t6, NULL );

	{athread_spawn(SpikeDeliver,swInfo);
	athread_join();
	{gettimeofday( &t7, NULL );ed = rpcc();}}

	tm1=&t4; tm2=&t5;
	DO_T5+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
	tm1=&t5; tm2=&t6;
	DO_T6+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
	tm1=&t6; tm2=&t7;
	DO_T7+=((double)(tm2->tv_sec-tm1->tv_sec)+(double)(tm2->tv_usec-tm1->tv_usec)*1e-6);
	return;
}

static bool updateTime(snnInfo_t *snnInfo) 
{
	bool finishedOneSec = false;
	if(++snnInfo->simTimeMs == 1000) 
	{
		snnInfo->simTimeMs = 0;
		snnInfo->simTimeSec++;
		finishedOneSec = true;
	}

	return finishedOneSec;
}

