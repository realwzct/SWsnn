#include "mysnn.h"
#include <sys/time.h>    // for gettimeofday() chenged!!
#include <stdlib.h>
int timestep;
float currentfactor;

int main(int argc, char *argv[])
{

	int connectnumber = atoi(argv[1]);  			//connect number 1 0000
	int neuronnumber = atoi(argv[2]);				//neuron number 20 0000
	int delay = atoi(argv[3]);          			//max delay
	timestep = (int)(1/atof(argv[4]));              //time step
	int simulatetime = atoi(argv[5]);				//simulatetime
	currentfactor = atof(argv[6]);

	//printf("%f\n",currentfactor);
	int *send, *rec;
	unsigned long  size1, size2;
    int rank, nproc;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numGrp=3, numConn=2;
	int randSeed = 42;
	int NE = neuronnumber; //neuron number 20 0000
	int NI=0; int NP=0;

	int NReg=NE+NI; int NInput=NP; int NPre=NReg+NInput;

	snnInfo_t snnInfo;
	grpInfo_t *grpInfo=(grpInfo_t*)malloc(numGrp*sizeof(grpInfo_t));
	connInfo_t *connInfo=(connInfo_t*)malloc(numConn*sizeof(connInfo_t));
	initNetwork(&snnInfo,grpInfo,connInfo,numGrp,numConn,randSeed);

	struct timeval t0,t1,t2;
	struct timeval *time1,*time2;
	gettimeofday( &t0, NULL );


	float a,b,c,d,a_sd,b_sd,c_sd,d_sd;
	
	/*******grp para*****/
	int i=0;
	grpInfo[i].SizeN = NE;
	grpInfo[i].StartN= 0; 
	grpInfo[i].EndN = grpInfo[i].StartN+grpInfo[i].SizeN-1;
	a   =0.02; b   =0.2; c   =-65; d   =8.0;
	a_sd=0.0 ; b_sd=0.0; c_sd=0.0;  d_sd=0.0;
	setNeuronParameters(grpInfo,i,a,a_sd,b,b_sd,c,c_sd,d,d_sd);

	/*******grp para*****/
	i++;
	grpInfo[i].SizeN = NI;
	grpInfo[i].StartN= grpInfo[i-1].EndN+1;
	grpInfo[i].EndN = grpInfo[i].StartN+grpInfo[i].SizeN-1;
	a   =0.02; b   =0.2; c   =-65; d   =8.0;
	a_sd=0.0 ; b_sd=0.0; c_sd=0.0;  d_sd=0.0;
	setNeuronParameters(grpInfo,i,a,a_sd,b,b_sd,c,c_sd,d,d_sd);

	i++;
	grpInfo[i].SizeN = NP;
	grpInfo[i].StartN= grpInfo[i-1].EndN+1;
	grpInfo[i].EndN = grpInfo[i].StartN+grpInfo[i].SizeN-1;
	a   =0.02; b   =0.2; c   =-65; d   =8.0;
	a_sd=0.0 ; b_sd=0.05; c_sd=0.0;  d_sd=0.0;
	setNeuronParameters(grpInfo,i,a,a_sd,b,b_sd,c,c_sd,d,d_sd);
	grpInfo[i].isSpikeGenerator=1;//for poisson

	snnInfo.numGrp=i+1;
	snnInfo.numN=grpInfo[i].EndN+1;
	snnInfo.numNExcReg=grpInfo[0].SizeN;
	snnInfo.numNInhReg=grpInfo[1].SizeN;
	snnInfo.numNReg=grpInfo[0].SizeN+grpInfo[1].SizeN;
	snnInfo.numNPois=grpInfo[2].SizeN;
	assert(snnInfo.numNReg==NReg);
	assert(snnInfo.numNPois==NInput);
	assert(snnInfo.numNExcReg==NE);
	assert(snnInfo.numNInhReg==NI);
	assert(snnInfo.numN==NPre);

	snnInfo.Nsyn = connectnumber; //connection number

	int minD,maxD; double minW,maxW,initW;
	int gSrc,gDest;
	/*******conn para******/
	int DL = delay;
	assert(DL>=2);
	minD=1; maxD=DL; minW=0.0;initW=0.000001;maxW=0.000001;
	gSrc=0;gDest=0;
	int j=0;
	connInfo[j].connId=j;
	connInfo[j].grpSrc=gSrc;
	connInfo[j].grpDest=gDest;
	connInfo[j].maxDelay=maxD;
	connInfo[j].initWt=initW;
	connInfo[j].maxWt=maxW;
	connInfo[j].p=1.;
	connInfo[j].maxSynPost=grpInfo[gDest].SizeN-1;
	connInfo[j].maxSynPre=grpInfo[gSrc].SizeN-1;
	connInfo[j].numSyn=grpInfo[gSrc].SizeN*(grpInfo[gSrc].SizeN-1);
	
	/*******conn para******/
	//NP->NE
	minD=1; maxD=DL; minW=0.000;initW=0.0005;maxW=0.0005;
	gSrc=1;gDest=0;
	j++;
	connInfo[j].connId=j;
	connInfo[j].grpSrc=gSrc;
	connInfo[j].grpDest=gDest;
	connInfo[j].maxDelay=maxD;
	connInfo[j].initWt=initW;
	connInfo[j].maxWt=maxW;
	connInfo[j].p=1.;
	connInfo[j].maxSynPost=grpInfo[gDest].SizeN;
	connInfo[j].maxSynPre=grpInfo[gSrc].SizeN;
	connInfo[j].numSyn=grpInfo[gSrc].SizeN*grpInfo[gDest].SizeN;
	snnInfo.maxDelay=maxD;//

	//COBA or CUBA (true or false)
	int tdAMPA=5,tdNMDA=150,tdGABAa=6,tdGABAb=150;
	setConductances(&snnInfo,1,tdAMPA,0,tdNMDA,tdGABAa,0,tdGABAb);
	buildModel(&snnInfo,grpInfo,connInfo);
	createNetwork(&snnInfo,grpInfo,connInfo);
	snnInfo_t *s=&snnInfo;
	MPI_Barrier(MPI_COMM_WORLD);//mpi syn
	long tm0,tm1;
	{gettimeofday( &t1, NULL );
	tm0=rpcc();}
	int nmsec = simulatetime;              //simulatetime
	runNetwork(&snnInfo,grpInfo,connInfo,snnInfo.nInfoHost,
               snnInfo.sInfoHost,snnInfo.swInfo,nmsec,0);

	{gettimeofday( &t2, NULL );tm1=rpcc();}

	if(rank==0)
	{
		time1=&t1; time2=&t2;
		printf("running time: %f\n", (double)(time2->tv_sec-time1->tv_sec)+(double)(time2->tv_usec-time1->tv_usec)*1e-6);

		time1=&t0; time2=&t2;
		//printf("all time : %f\n", (double)(time2->tv_sec-time1->tv_sec)+(double)(time2->tv_usec-time1->tv_usec)*1e-6);

		//printf("main::cycle: %ld\n",tm1-tm0);
	}

    MPI_Finalize();
	return 0;
}
