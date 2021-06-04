#include "mysnn.h"
#include <sys/stat.h> // mkdir
#include <sys/time.h>    // for gettimeofday() chenged!!
#include <math.h>

#define COMPACTION_ALIGNMENT_PRE  16
#define COMPACTION_ALIGNMENT_POST 0

extern int timestep;
static float getWeights(grpInfo_t *gInfo,int connProp, float initWt, float maxWt, unsigned int nid, int gid);
static void connectFull(snnInfo_t *sInfo,connInfo_t* cInfo,grpInfo_t *gInfo); 
static void connectPart(snnInfo_t *sInfo,connInfo_t* cInfo,grpInfo_t *gInfo); 
static void resetNeuron(neurInfo_t *nInfo,int nid,grpInfo_t *gInfo,int gid);
static void resetSpikeCnt(snnInfo_t* snnInfo,grpInfo_t* grpInfo,int grpId);
static void resetSynapticConnections();
static void updateSpikeGeneratorsInit(grpInfo_t* grpInfo, int numGrp); 
static int updateSpikeTables(snnInfo_t *s);


void setConductances(snnInfo_t *sInfo, bool isSet, int tdAMPA, int trNMDA, int tdNMDA, int tdGABAa, int trGABAb, int tdGABAb)
{
	if (isSet) 
	{
		assert(tdAMPA>0); assert(tdNMDA>0); assert(tdGABAa>0); assert(tdGABAb>0);
		assert(trNMDA>=0); assert(trGABAb>=0); // 0 to disable rise times
		assert(trNMDA!=tdNMDA); assert(trGABAb!=tdGABAb); // singularity
	}

	// set conductances globally for all connections
	sInfo->sim_with_conductances  |= isSet;
	sInfo->dAMPA  = 1.0-1.0/tdAMPA;
	sInfo->dNMDA  = 1.0-1.0/tdNMDA;
	sInfo->dGABAa = 1.0-1.0/tdGABAa;
	sInfo->dGABAb = 1.0-1.0/tdGABAb;
}

// set Izhikevich parameters for group
void setNeuronParameters(grpInfo_t *gInfo,int gid,float izh_a,float izh_a_sd,float izh_b,float izh_b_sd,float izh_c,float izh_c_sd,float izh_d,float izh_d_sd)
{
	assert(gid>=0);
	assert(izh_a_sd>=0);assert(izh_b_sd>=0);
	assert(izh_c_sd>=0);assert(izh_d_sd>=0);

	gInfo[gid].Izh_a = izh_a;
	gInfo[gid].Izh_a_sd  =  izh_a_sd;
	gInfo[gid].Izh_b = izh_b;
	gInfo[gid].Izh_b_sd = izh_b_sd;
	gInfo[gid].Izh_c = izh_c;
	gInfo[gid].Izh_c_sd = izh_c_sd;
	gInfo[gid].Izh_d = izh_d;
	gInfo[gid].Izh_d_sd = izh_d_sd;
}

static void connectFull(snnInfo_t *sInfo,connInfo_t* cInfo,grpInfo_t *gInfo) {
	int gSrc=cInfo->grpSrc;
	int gDest=cInfo->grpDest;
	uint8_t minD=1;
	uint8_t maxD=cInfo->maxDelay;
	assert(maxD-1==sInfo->Ndelay);
	bool noDirect = 1;
	float synWt;
	int preN=sInfo->preN;
	int Nsyn=sInfo->Nsyn;
	int Ndelay=sInfo->Ndelay;
	int MaxN=sInfo->MaxN;
	int SizeN=sInfo->swInfo[0].SizeN;
	printf("connectFull::preN=%d Ndelay=%d NTh=%d MaxN=%d SizeN=%d\n",preN,Ndelay,NTh,MaxN,SizeN);
	/******init*****/
	long is;
	long len=sInfo->preN*Ndelay*NTh*MaxN;
	printf("connectFull::sInfoHost %d\n",sInfo->sInfoHost);
	printf("connectFull::%ld %d %d %d\n",len,preN,NTh,MaxN);
	for(is=0;is<len;is++){
		sInfo->sInfoHost[is].postId=0xffff;
		sInfo->sInfoHost[is].dl=0;
	}
	/***************/
	
	int MIND=minD<<sInfo->Nop;
	assert(MIND==sInfo->Ndt);
	int MAXD=maxD<<sInfo->Nop;
	assert((MAXD-MIND)%sInfo->Ndt==0);
	int i,j;
	int Ndma[NTh];for(i=0;i<NTh;i++)Ndma[i]=0;
	for(i=0;i<sInfo->preN;i++)  {
		int offset=i*Ndelay*NTh*MaxN;
		if(i==sInfo->numNReg) cInfo++;
		float fac=1.0;
		if(i>=gInfo[1].StartN && i<sInfo->numNReg) fac=-1.0;

		int iD=0,iTh=0,iN[32];
		for(iD=0;iD<Ndelay;iD++)iN[iD]=0;
		for(j=0;j<sInfo->numNReg;j++) {
			if(j>=sInfo->swInfo[iTh].StartN+sInfo->swInfo[iTh].SizeN){
				iTh++;
				for(iD=0;iD<Ndelay;iD++)iN[iD]=0;
				
			}
			if((noDirect)&&i==j) continue;

			uint8_t dVal;
			for(;;)
			{
				dVal=MIND+rand()%(MAXD-MIND);//for delay!!
				dVal=((i+j)%sInfo->Ndelay+1)<<sInfo->Nop;//for test
				assert((dVal>=MIND)&&(dVal<=MAXD));
				iD = (dVal-MIND)>>sInfo->Nop;
				assert(iD<sInfo->Ndelay && iD>=0);
				if(iN[iD]<MaxN) break;
				else {
				}
			}

			synWt=cInfo->initWt;
			assert(synWt>=0.);
			int addr=offset+iD*NTh*MaxN+iTh*MaxN+iN[iD];
			sInfo->sInfoHost[addr].postId=j;
			sInfo->sInfoHost[addr].wt=synWt;
			sInfo->sInfoHost[addr].dl=dVal;
			iN[iD]++;
			if(iN[iD]>Ndma[iTh]) Ndma[iTh]=iN[iD];
		}
	}
	for(i=0;i<NTh;i++)sInfo->swInfo[i].Ndma = Ndma[i];
	printf("connectFull::Ndma=%d Ndma2=%d MaxN=%d\n",Ndma[0],Ndma[NTh-1],MaxN);
	return;
}

void initNetwork(snnInfo_t *sInfo, grpInfo_t *gInfo,connInfo_t *cInfo,
		int numGrp,int numConn,int randSeed){
	int i,g;

	srand48(randSeed); sInfo->randSeed=randSeed;

	sInfo->simTimeMs = 0; sInfo->simTimeSec = 0; sInfo->simTimeMs = 0;

	sInfo->maxDelay = 1; //sInfo->maxDelay_= 1;

	sInfo->numGrp = numGrp;	sInfo->numConn= numConn;

	sInfo->numN = 0; sInfo->numNReg = 0; 
	sInfo->numNExcReg = 0; sInfo->numNInhReg = 0; sInfo->numNPois = 0;

	sInfo->maxSpikesD1 = 0; sInfo->maxSpikesD2 = 0;

	sInfo->simTimeRunStart = 0; sInfo->simTimeRunStop = 0;
	
	sInfo->secD1fireCntHost = 0; sInfo->secD2fireCntHost = 0;
	sInfo->spikeCountAll1secHost = 0;
	
	sInfo->MaxFiringRate = 60;//60Hz ????
	sInfo->sim_with_homeostasis = 0; sInfo->sim_with_conductances = 1;
	sInfo->sim_with_NMDA_rise = 0; sInfo->sim_with_GABAb_rise = 0;

	// some default decay and rise times
	sInfo->dAMPA = 1.0-1.0/5.0;
	sInfo->rNMDA = 1.0-1.0/10.0; sInfo->dNMDA = 1.0-1.0/150.0;
	sInfo->sNMDA = 1.0;
	sInfo->dGABAa = 1.0-1.0/6.0;
	sInfo->rGABAb = 1.0-1.0/100.0; sInfo->dGABAb = 1.0-1.0/150.0;
	sInfo->sGABAb = 1.0;

	sInfo->timeTableD1 = NULL; sInfo->timeTableD2 = NULL;
	sInfo->firingTableD1 = NULL; sInfo->timeTableD2 = NULL;

	sInfo->nSpikeCnt = NULL; sInfo->lastSpikeTime = NULL;
	
	for (i=0;i<numGrp;i++)
	{
		gInfo[i].Type = 0;
		gInfo[i].StartN = -1;gInfo[i].EndN = -1;gInfo[i].SizeN = -1;
		
		gInfo[i].MaxDelay = sInfo->maxDelay;
		gInfo[i].FiringCount1sec = 0;

		gInfo[i].isSpikeGenerator= 0;gInfo[i].NewTimeSlice = 0;
		gInfo[i].CurrTimeSlice = 0;gInfo[i].SliceUpdateTime = 0;
		gInfo[i].RefractPeriod = 0;gInfo[i].RatePtr = NULL;

		gInfo[i].WithSTP = 0;gInfo[i].WithHomeostasis = 0;
		
		gInfo[i].numSynPre = 0;	gInfo[i].numSynPost = 0;

		gInfo[i].Izh_a = 0.; gInfo[i].Izh_a_sd = 0.;
		gInfo[i].Izh_b = 0.; gInfo[i].Izh_b_sd = 0.;
		gInfo[i].Izh_c = 0.; gInfo[i].Izh_c_sd = 0.;
		gInfo[i].Izh_d = 0.; gInfo[i].Izh_d_sd = 0.;

	}
	
	for(i=0;i<numConn;i++){
		cInfo[i].connId = i;
		cInfo[i].maxDelay = sInfo->maxDelay;
		cInfo[i].grpSrc = -1;cInfo[i].grpDest = -1;
		cInfo[i].initWt = 0.;cInfo[i].maxWt = 0.;
		cInfo[i].mulSynFast = 1.;cInfo[i].mulSynSlow = 1.;
		cInfo[i].maxSynPost = -1;cInfo[i].maxSynPre = -1;
		cInfo[i].numSyn = -1;cInfo[i].connProp = 0;cInfo[i].p = 0.;
	}
	
	return;
}

static void connectPart(snnInfo_t *sInfo,connInfo_t* cInfo,grpInfo_t *gInfo) {
	int rank, nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	uint8_t minD=1;
	uint8_t maxD=cInfo->maxDelay;
	assert(maxD-1==sInfo->Ndelay);
	bool noDirect = 0;
	float synWt;
	int preN=sInfo->preN;
	int Nsyn=sInfo->Nsyn;
	int Ndelay=sInfo->Ndelay;
	int MaxN=sInfo->MaxN;
	int SizeN=sInfo->swInfo[0].SizeN;

	//printf("connectPart::preN=%d Ndelay=%d NTh=%d MaxN=%d SizeN=%d\n",preN,Ndelay,NTh,MaxN,SizeN);
	/******init*****/
	long is;
	long len=sInfo->preN*Ndelay*NTh*MaxN;

	for(is=0;is<len;is++){
		sInfo->sInfoHost[is].postId=0xffff;
		sInfo->sInfoHost[is].dl=0;
	}

	
	int MIND=minD<<sInfo->Nop;
	assert(MIND==sInfo->Ndt);
	int MAXD=maxD<<sInfo->Nop;
	assert((MAXD-MIND)%sInfo->Ndt==0);
	int i,j;
	int Ndma[NTh];for(i=0;i<NTh;i++)Ndma[i]=0;
	for(i=0;i<sInfo->preN;i++)  {
		int offset=i*Ndelay*NTh*MaxN;
		if(i==sInfo->numNReg) cInfo++;

		float fac=1.0;
		if(i>=gInfo[1].StartN && i<sInfo->numNReg) fac=-1.0;

		int iD=0,iTh=0,iN[20];
		assert(sInfo->Ndelay<=20);
		for(iD=0;iD<Ndelay;iD++)iN[iD]=0;

		int NNN; //select 1,2,3,4
		NNN = sInfo->numNReg/sInfo->Nsyn;
		assert(NNN>0);
		for(j=i%NNN;j<sInfo->numNReg;j+=NNN) {
			if(j<sInfo->gStart||j>=sInfo->gStart+sInfo->gSize) continue;
			if(j>=sInfo->swInfo[iTh].StartN+sInfo->swInfo[iTh].SizeN){
				iTh++;
				for(iD=0;iD<Ndelay;iD++)iN[iD]=0;
				
			}
			if((noDirect)&&i==j) continue;

			uint8_t dVal;
			for(;;){

				dVal=((i+j/NNN)%sInfo->Ndelay+1)<<sInfo->Nop;//for test
				assert((dVal>=MIND)&&(dVal<MAXD));
				iD = (dVal-MIND)>>sInfo->Nop;
				assert(iD<sInfo->Ndelay && iD>=0);
				if(iN[iD]<MaxN) break;
				else {
					break;
				}
			}
			synWt=cInfo->initWt;
			assert(synWt>=0.);
			int addr=offset+iD*NTh*MaxN+iTh*MaxN+iN[iD];
			sInfo->sInfoHost[addr].postId=j;
			sInfo->sInfoHost[addr].wt=synWt;
			sInfo->sInfoHost[addr].dl=dVal;
			iN[iD]++;
			if(iN[iD]>Ndma[iTh]) Ndma[iTh]=iN[iD];
		}

	}
	for(i=0;i<NTh;i++)sInfo->swInfo[i].Ndma = Ndma[i];
	return;
}//connectPart

void buildModel(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo){//????
	int i,j;
	int rank, nproc;
	
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	sInfo->gSize = (sInfo->numNReg-1)/nproc+1;
	sInfo->gStart = sInfo->gSize*rank;
	if(rank=nproc-1)sInfo->gSize = sInfo->numNReg-(nproc-1)*sInfo->gSize;
	int SizeN=(sInfo->gSize-1)/NTh+1;
	assert(SizeN*NTh>=sInfo->gSize);

	if (sInfo->Nsyn==0) sInfo->Nsyn=10000;
	if (sInfo->Nsyn>sInfo->numNReg) sInfo->Nsyn=sInfo->numNReg;

	sInfo->Ndelay=sInfo->maxDelay-1;
	if(sInfo->Ndelay==0) assert(0);
	
	sInfo->Ndt = timestep;
	//sInfo->Nop = (int)(log2(timestep));
	sInfo->Nop = (int)(log(timestep)/log(2));
	sInfo->dt=1.0/sInfo->Ndt;
	

	sInfo->simTime = 0;
	sInfo->preN = sInfo->numN;

	sInfo->MaxN = (sInfo->Nsyn/NTh/sInfo->Ndelay/nproc/4+1)*4;
	//printf("preN=%d Nsyn=%d Ndelay=%d MaxN=%d\n",sInfo->preN,sInfo->Nsyn,sInfo->Ndelay,sInfo->MaxN);

	for(i=0;i<NTh;i++){
		int SizeN=(sInfo->gSize)/NTh;
		int least=sInfo->gSize%NTh;
		if(i<least){
			sInfo->swInfo[i].StartN = sInfo->gStart+i*(SizeN+1);
			sInfo->swInfo[i].SizeN = SizeN+1;
		}
		else{
			sInfo->swInfo[i].StartN = sInfo->gStart+i*(SizeN)+least;
			sInfo->swInfo[i].SizeN = SizeN;
		}
		sInfo->swInfo[i].preN = sInfo->preN;
		sInfo->swInfo[i].Ndelay = sInfo->Ndelay;
		sInfo->swInfo[i].MaxN = sInfo->MaxN;
		sInfo->swInfo[i].Ndt = sInfo->Ndt;
		sInfo->swInfo[i].Nop = sInfo->Nop;
		sInfo->swInfo[i].dt = sInfo->dt;
		sInfo->swInfo[i].rank = rank;
		sInfo->swInfo[i].nproc = nproc;
		sInfo->swInfo[i].gStart = sInfo->gStart;
		sInfo->swInfo[i].gSize = sInfo->gSize;
	}

	sInfo->nInfoHost=(neurInfo_t*)malloc(sInfo->gSize*sizeof(neurInfo_t));
	sInfo->sInfoHost=NULL;	
	sInfo->sInfoHost=(synInfo_t*)malloc((long)sInfo->preN*sInfo->Ndelay*sInfo->MaxN*NTh*sizeof(synInfo_t));

	//mpi mem alloc
	sInfo->NNgroup = sInfo->gSize;
	sInfo->NN = sInfo->NNgroup*nproc;
	sInfo->ND = sInfo->Ndelay;
	sInfo->firingTableHost=(spikeTime_t*)malloc(sInfo->NNgroup*sizeof(spikeTime_t));
	sInfo->firingTableAll=(spikeTime_t*)malloc(sInfo->NN*sInfo->ND*sizeof(spikeTime_t));

	sInfo->displs = (int *) malloc(sizeof(int) * nproc); // displs array
    sInfo->recvCount = (int *) malloc(sizeof(int) * nproc); // send num array

#if 1
{
	spikeTime_t *sendBuf, *recvBuf;
	int *displs, *recvCount;
	int root = 0;
	int senddatanum,nsall=0;

	sendBuf = &(sInfo->firingTableHost[0]);
	senddatanum= sInfo->gSize;

	displs = sInfo->displs; // displs
	recvCount = sInfo->recvCount; // send num

	MPI_Gather(&senddatanum,1,MPI_INT,recvCount,1, MPI_INT,root,MPI_COMM_WORLD);

	if(!rank)
	{
		displs[0] = 0;
		nsall = recvCount[0];
			for(i=1;i<nproc;i++)
			{
				displs[i] = displs[i-1]+recvCount[i-1];
				nsall += recvCount[i];
			}
	}

	recvBuf = &(sInfo->firingTableAll[0]);
	
	assert(sizeof(spikeTime_t)==sizeof(MPI_INT));
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
 	MPI_Bcast(&nsall,1,MPI_INT,root,MPI_COMM_WORLD);
	MPI_Bcast(recvBuf,nsall,MPI_INT,root,MPI_COMM_WORLD);
#endif 
}
#endif
	/*******set syn and neur*******/
	int gid;
	assert(gInfo[0].StartN==0);
	assert(gInfo[0].EndN+1==gInfo[1].StartN);
	assert(gInfo[1].EndN+1==sInfo->numNReg);
	assert(sInfo->numGrp==3);
	for(gid=0;gid<2;gid++){
		for(i=gInfo[gid].StartN;i<=gInfo[gid].EndN;i++){
			if(i<sInfo->gStart||i>=sInfo->gStart+sInfo->gSize) continue;
			resetNeuron(sInfo->nInfoHost,i-sInfo->gStart,gInfo,gid);//nInfo
		}
	}

	if(sInfo->Nsyn==sInfo->numNReg)	{
		connectPart(sInfo,&cInfo[0],gInfo);
	}else if(sInfo->Nsyn<sInfo->numNReg){
	
		connectPart(sInfo,&cInfo[0],gInfo);
	}else {printf("Error::Nsyn can't be over Nreg\n");assert(0);}
	return;
}

void createNetwork(snnInfo_t *sInfo,grpInfo_t *gInfo,connInfo_t *cInfo){
	int i;
	int SizeN=(sInfo->numNReg-1)/NTh+1;
	sInfo->preN = sInfo->numN;

	for(i=0;i<NTh;i++){
		if(sInfo->swInfo[i].StartN>=sInfo->numNReg){assert(0);}

		sInfo->swInfo[i].dAMPA = sInfo->dAMPA;
		sInfo->swInfo[i].dNMDA = sInfo->dNMDA;
		sInfo->swInfo[i].rNMDA = sInfo->rNMDA;
		sInfo->swInfo[i].sNMDA = sInfo->sNMDA;
		sInfo->swInfo[i].dGABAa = sInfo->dGABAa;
		sInfo->swInfo[i].dGABAb = sInfo->dGABAb;
		sInfo->swInfo[i].rGABAb = sInfo->rGABAb;
		sInfo->swInfo[i].sGABAb = sInfo->sGABAb;
		sInfo->swInfo[i].sim_with_conductances = sInfo->sim_with_conductances;
		sInfo->swInfo[i].WithSTDP = 0;
		sInfo->swInfo[i].WithSTP = 0;
		sInfo->swInfo[i].WithHomeostasis = 0;
		sInfo->swInfo[i].sim_with_homeostasis = sInfo->sim_with_homeostasis;
		sInfo->swInfo[i].sim_with_NMDA_rise = sInfo->sim_with_NMDA_rise;
		sInfo->swInfo[i].sim_with_GABAb_rise = sInfo->sim_with_GABAb_rise;
		sInfo->swInfo[i].MaxDelay = sInfo->maxDelay;
		sInfo->swInfo[i].nSpikePoisAll = 0;
		sInfo->swInfo[i].fireCnt = 0;
		sInfo->swInfo[i].nInfoHost = sInfo->nInfoHost;
		sInfo->swInfo[i].sInfoHost = sInfo->sInfoHost;

		//mpi firingTable
		sInfo->swInfo[i].NNgroup = sInfo->NNgroup;
		sInfo->swInfo[i].NN = sInfo->NN;
		sInfo->swInfo[i].ND = sInfo->ND;
		sInfo->swInfo[i].firingTableHost = sInfo->firingTableHost;
		sInfo->swInfo[i].firingTableAll = sInfo->firingTableAll;
	}
	return;
}

void setupNetwork(snnInfo_t *sInfo,grpInfo_t *gInfo) 
{
}

static void resetNeuron(neurInfo_t *nInfo,int n,grpInfo_t *gInfo,int g) 
{
	nInfo[n].Izh_a=gInfo[g].Izh_a+gInfo[g].Izh_a_sd*(float)drand48();
	nInfo[n].Izh_b=gInfo[g].Izh_b+gInfo[g].Izh_b_sd*(float)drand48();
	nInfo[n].Izh_c=gInfo[g].Izh_c+gInfo[g].Izh_c_sd*(float)drand48();
	nInfo[n].Izh_d=gInfo[g].Izh_d+gInfo[g].Izh_d_sd*(float)drand48();
	nInfo[n].voltage=nInfo[n].Izh_c;// initial values for new_v
	nInfo[n].recovery=nInfo[n].Izh_b*nInfo[n].voltage;//for new_u
	nInfo[n].gAMPA = 0.;
	nInfo[n].gNMDA_d = 0.;
	nInfo[n].gGABAa = 0.;
	nInfo[n].gGABAb_d = 0.;
}



