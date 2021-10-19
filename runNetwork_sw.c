#include "mysnn.h"
#include "mpi.h"
#include <athread.h>

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
extern int rank,nproc;
volatile static long t1=0,t2=0,t3=0;

int runNetwork(snnInfo_t *sInfo, grpInfo_t *gInfo,connInfo_t *cInfo, neurInfo_t *nInfo,
    synInfo_t *synInfo, swInfo_t *swInfo,int _nmsec, bool printRun) 
{
	int runDurMs = _nmsec;

	sInfo->simTimeRunStart = sInfo->simTime;
	sInfo->simTimeRunStop  = sInfo->simTime+runDurMs;

	//initialize athread!!!!
	athread_spawn(initSW,swInfo);
    athread_join();
	int i;
	nspikeall=0;
	for(i = 0; i<runDurMs; i++) 
	{
		doSnnSim(sInfo,gInfo,cInfo,nInfo,synInfo,swInfo);
		updateTime(sInfo);
	}
	athread_spawn(freeSW,swInfo);
    athread_join();
	

	unsigned long cdma=0,cspike=0;
	for(i=0;i<63;i++){
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
		printf("neuron state update: %.4lf ms\n", (double)(t1)*1000/CLOCKRATE);
		printf("mpi time: %.4lf ms\n", (double)(t2)*1000/CLOCKRATE);
		printf("spike deliver: %.4lf ms\n", (double)(t3)*1000/CLOCKRATE);
		//printf("sInfo->NN%d",sInfo->NN);//20000神经元
	}

	return 0;
}

static void doSnnSim(snnInfo_t *sInfo, grpInfo_t *gInfo,connInfo_t *cInfo, 
    neurInfo_t *nInfo,synInfo_t *synInfo, swInfo_t *swInfo) 
{
	volatile long time0,time1,time2,time3;

	{time0 = rpcc();}
	athread_spawn(StateUpdate,swInfo);
	athread_join();
	{time1 = rpcc();}
	sInfo->simTime++;

	int iND = sInfo->simTime%sInfo->Ndelay;  //IND=0;


	spikeTime_t *sendBuf, *recvBuf;
	int *displs, *recvCount;
	int root = 0;
	int senddatanum;

	sendBuf = &(sInfo->firingTableHost[0]);
	senddatanum= NS_group;
	displs = sInfo->displs; // displs
	recvCount = sInfo->recvCount; // send num

	MPI_Gather(&senddatanum,1,MPI_INT,recvCount,1, MPI_INT,root,MPI_COMM_WORLD);

	if(!rank){
		displs[0] = 0;
		NSall = recvCount[0];
		int i;
		for(i=1;i<nproc;i++){	
			//printf("%d %d %d %d %d %d %d %d %d %d rank%d\n",recvCount[0],recvCount[1],recvCount[2],recvCount[3],recvCount[4],recvCount[5],recvCount[6],recvCount[7],recvCount[8],recvCount[9],rank);
			displs[i] = displs[i-1]+recvCount[i-1];
			NSall += recvCount[i];//总数量
		}
	}

	recvBuf = &(sInfo->firingTableAll[iND*sInfo->NN]);  //iND

	MPI_Gatherv(sendBuf,
		senddatanum,
		MPI_INT,
		recvBuf,
		recvCount, 
		displs, 
		MPI_INT,
		root,
		MPI_COMM_WORLD);

	// if(rank==0){
	// 	int m1=0;
	// 	for(m1=0;m1<10;m1++){
	// 		printf("%d ",sInfo->firingTableAll[m1].nid);
	// 		printf("%d ",sInfo->firingTableAll[m1].time);
	// 	}
	// 	printf("\n");
	// }
	if(rank == 0){
		//printf("NSall\n%d",NSall);
	}

	MPI_Bcast(&NSall,1,MPI_INT,root,MPI_COMM_WORLD);
	MPI_Bcast(recvBuf,NSall,MPI_INT,root,MPI_COMM_WORLD);
	nspikeall += NSall;

	{time2 = rpcc();}

	athread_spawn(SpikeDeliver,swInfo);
	athread_join();

	{time3 = rpcc();}

	t1 = time1-time0+t1;
	t2 = time2-time1+t2;
	t3 = time3-time2+t3;
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

