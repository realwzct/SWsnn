
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"
#include "my_slave.h"
#include <assert.h>

#include "swstruct.h"

#define COND_INTEGRATION_SCALE        2

static __inline long rpcc()
{
   long a;
   asm volatile ("memb");
   asm volatile ("rcsr %0,4":"=r"(a));
   return a;
}
 #define REG_SYNR(mask) \
  asm volatile ("synr %0"::"r"(mask))
 #define REG_SYNC(mask) \
  asm volatile ("sync %0"::"r"(mask))

__thread_local volatile unsigned int reply,rpl[8];
__thread_local volatile unsigned int reply_mpi,rpl_mpi[8];
__thread_local volatile unsigned int st0,ed0,result;
__thread_local volatile unsigned int cdma,cspike;
__thread_local volatile unsigned int cdma_mpi,cspike_mpi;
extern volatile unsigned int dma[64],spike[64],NS_group,NSall,numSpike[64];
extern float currentfactor;
__thread_local volatile unsigned int NSgroup_slave,NSall_slave;

__thread_local swInfo_t swInfo;
__thread_local neurInfo_t *nInfo;
__thread_local synInfo_t *sInfo;
__thread_local int time_test=0;
__thread_local dma_desc dma_get_syn;
__thread_local dma_desc dma_get;

static void generatePostSpike(spikeTime_t st,synInfo_t *sInfoLc);
static void generatePostSpike_simd(spikeTime_t st,synInfo_t *sInfoLc);
static int addSpikeToTable(int i);
static int addSpikeToTable_simd(int,int);
static int addSpikeToTable_mpi(int i);
static int addSpikeToTable_simd_mpi(int,int);

__thread_local spikeTime_t *firingTable;
__thread_local spikeTime_t *firingTable_mpi;
__thread_local float *ringBuffer;
__thread_local int simTime,sliceTime;
__thread_local int simTime_mpi,sliceTime_mpi;
__thread_local int lenRB,offset;
__thread_local int lenRB_mpi,offset_mpi;

__thread_local int lenST,topST,endST,usedST,numST[20];
__thread_local volatile int lenST_mpi,topST_mpi,endST_mpi,usedST_mpi,numST_mpi[20];

static void decayConduct(void *ptr);
static void neuronUpdate(int it);
static void neuronUpdate_simd();
static void CurrentUpdate(void *ptr);
static void CurrentUpdate_mpi(void *ptr);
static void PoisCurrentUpdate(void *ptr);

void initSW(swInfo_t *ptr){

	int i;
	cdma=0;cspike=0;
	cdma_mpi=0;cspike_mpi=0;
	simTime=0,sliceTime=0;lenRB=0;
	simTime_mpi=0,sliceTime_mpi=0;lenRB_mpi=0;
	lenST=0; topST=0; endST=0; usedST=0;
	lenST_mpi=0; topST_mpi=0; endST_mpi=0; usedST_mpi=0;
	for(i=0;i<20;i++) numST[i]=0;
	for(i=0;i<20;i++) numST_mpi[i]=0;

	reply=0;
	athread_get(PE_MODE,&ptr[_MYID],&swInfo,sizeof(swInfo_t),&reply,0,0,0);
	while(reply!=1);

	int SizeN=swInfo.SizeN;
	int MaxN=swInfo.MaxN;

	lenRB=2*swInfo.Ndt; offset=0;
	lenRB_mpi=2*swInfo.Ndt; offset_mpi=0;
	lenST_mpi=SizeN+64;//spikeTime buffer

	if(_MYID==NTh-1) lenST+=100;
	if(_MYID==NTh-1) lenST_mpi+=100;
	/*****allocate mem********/
	nInfo=(neurInfo_t*)ldm_malloc(SizeN*sizeof(neurInfo_t)); //?????
	sInfo=(synInfo_t*)ldm_malloc(8*swInfo.Ndma*sizeof(synInfo_t));
	firingTable_mpi=(spikeTime_t*)ldm_malloc(lenST_mpi*sizeof(spikeTime_t));
	ringBuffer=(float*)ldm_malloc(lenRB*SizeN*sizeof(float));
	assert(sizeof(spikeTime_t)==4);
	assert(sizeof(synInfo_t)==8);
	for(i=0;i<lenST;i++){
		firingTable_mpi[i].time=0xFFFF;
		firingTable_mpi[i].nid =0xFFFF;
	}
	for(i=0;i<lenRB*SizeN;i++){ringBuffer[i] =0.;}
	memset(nInfo,0,SizeN*sizeof(neurInfo_t));//init zero

	dma_set_size(&dma_get_syn,swInfo.Ndma*sizeof(synInfo_t));
        dma_set_op(&dma_get_syn,DMA_GET);
        dma_set_reply(&dma_get_syn,&reply);
        dma_set_mode(&dma_get_syn,PE_MODE);
        dma_set_bsize(&dma_get_syn,0);
        dma_set_stepsize(&dma_get_syn,0);

	dma_set_size(&dma_get,SizeN*sizeof(neurInfo_t));//need reset
        dma_set_op(&dma_get,DMA_GET);
        dma_set_reply(&dma_get,&reply);
        dma_set_mode(&dma_get,PE_MODE);
        dma_set_bsize(&dma_get,0);
        dma_set_stepsize(&dma_get,0);

	reply=0;
	dma(dma_get,&(swInfo.nInfoHost[swInfo.StartN-swInfo.gStart]),&nInfo[0]);
        dma_wait(&reply,1);

	for(i=0;i<swInfo.SizeN;i++){
		nInfo[i].nSpikeCnt=0;
	}
        swInfo.nSpikePoisAll=0;
	return;	
}

void freeSW(void *ptr){//????
	int nSpikeAll=0;
	int i;

	for(i=0;i<swInfo.SizeN;i++){
		nSpikeAll += nInfo[i].nSpikeCnt;
	}

	reply=0;
	athread_put(PE_MODE,&cdma,&dma[_MYID],sizeof(int),&reply,0,0);
	athread_put(PE_MODE,&cspike,&spike[_MYID],sizeof(int),&reply,0,0);
	athread_put(PE_MODE,&nSpikeAll,&numSpike[_MYID],sizeof(int),&reply,0,0);
	while(reply!=3);
	
	return;
}

static int SpikeDmaWrite_mpi(ptr);//mpi
static int SpikeDmaRead_mpi(ptr);//mpi

void StateUpdate(void *ptr){
	int it;
	decayConduct(NULL);
	endST_mpi=0; topST_mpi=0; usedST_mpi=0;
	neuronUpdate_simd();
	sliceTime+=swInfo.Ndt;
	simTime++;
	assert(simTime==(sliceTime>>swInfo.Nop));
	offset += swInfo.Ndt;//ringBuffer offset????????
	if(offset>=lenRB) offset -= lenRB;
	SpikeDmaWrite_mpi(ptr);//mpi++++
	return;
}

static void decayConduct(void *ptr){

	int i=0;
        for(i=0; i<swInfo.SizeN; i++) {

                if (swInfo.sim_with_conductances) {
                        nInfo[i].gAMPA*=swInfo.dAMPA;
                        nInfo[i].gGABAa*=swInfo.dGABAa;
                        nInfo[i].gNMDA_d*=swInfo.dNMDA;//instantaneous rise
                        nInfo[i].gGABAb_d*=swInfo.dGABAb;//instantaneous rise
                }
                else {
                        nInfo[i].gAMPA=0.0f; //in CUBA current,sum up all wts
                }
        }

	return;
}

inline float dvdtIzh(float v, float u, float tmpI, float h);
inline float dudtIzh(float v, float u, float a, float b, float h);
#if 1
static void neuronUpdate(int it){

	double tmpiNMDA, tmpI;
	double tmpgNMDA, tmpgGABAb;
	int i=0,j;
	for(i=0;i<swInfo.SizeN;i++) {
#if 0

#else

		if (swInfo.sim_with_conductances) {
			float volt = -60.0;
			tmpiNMDA=(volt+80.0)*(volt+80.0)/60.0/60.0;

			tmpgNMDA=nInfo[i].gNMDA_d;
			tmpgGABAb=nInfo[i].gGABAb_d;

			tmpI=-(nInfo[i].gAMPA*(volt-0)
				 +tmpgNMDA*tmpiNMDA/(1+tmpiNMDA)*(volt-0)
				 +nInfo[i].gGABAa*(volt+70)
				 +tmpgGABAb*(volt+90));
		} else {
			tmpI=nInfo[i].gAMPA;
		}

#endif
		/* 4th Runge-Kutta */

		float v = nInfo[i].voltage;
		float u = nInfo[i].recovery;
		float h = swInfo.dt;
		float a = nInfo[i].Izh_a;
		float b = nInfo[i].Izh_b;

		float k1 = dvdtIzh(v, u, tmpI, h);
		float l1 = dudtIzh(v, u, a, b, h);

		float k2 = dvdtIzh(v+0.5*k1, u+0.5*l1, tmpI, h);
		float l2 = dudtIzh(v+0.5*k1, u+0.5*l1, a, b, h);

		float k3 = dvdtIzh(v+0.5*k2, u+0.5*l2, tmpI, h);
		float l3 = dudtIzh(v+0.5*k2, u+0.5*l2, a, b, h);

		float k4 = dvdtIzh(v+k3, u+l3, tmpI, h);
		float l4 = dudtIzh(v+k3, u+l3, a, b, h);
		v += (k1+2.0*k2+2.0*k3+k4)/6.0;
		if (v > 30.0) v = 30.0;
		if (v < -90.0) v = -90.0;

		u += (l1+2.0*l2+2.0*l3+l4)/6.0;

		nInfo[i].voltage = v;
		nInfo[i].recovery= u;
		/*** findFiring ***/
		if (nInfo[i].voltage>= 30.0) {
			nInfo[i].voltage = nInfo[i].Izh_c;
			nInfo[i].recovery+=nInfo[i].Izh_d;
			if(addSpikeToTable_mpi(i)) assert(0);
		}
	} 
#if 1
	int dIndex=offset+it;
	assert(dIndex<lenRB);
	for(i=0;i<swInfo.SizeN;i++){
		int addr = i*lenRB + dIndex;
		if (swInfo.sim_with_conductances) {
			nInfo[i].gAMPA += ringBuffer[addr];
			nInfo[i].gNMDA_d += ringBuffer[addr];
		} else {
			assert(0);
		}
		ringBuffer[addr]=0.;
	}
#endif

}
#endif
static int addSpikeToTable(int i) {
	int spikeBufferFull = 0;
	if(i<swInfo.SizeN) {nInfo[i].nSpikeCnt++;}
	firingTable[endST].nid = i+swInfo.StartN;
	firingTable[endST].time= sliceTime;
	endST++; usedST++; numST[simTime%swInfo.Ndelay]++;
	if(endST>=lenST) endST -= lenST;
	if(endST==topST) assert(0);
	if(usedST==lenST) assert(0);
	swInfo.fireCnt++;
	return spikeBufferFull;
}
static int addSpikeToTable_mpi(int i){
	int spikeBufferFull = 0;
	if(i<swInfo.SizeN) {nInfo[i].nSpikeCnt++;}
	firingTable_mpi[endST_mpi].nid = i+swInfo.StartN;
	firingTable_mpi[endST_mpi].time= sliceTime;
	endST_mpi++; usedST_mpi++; 
	if(endST_mpi>=lenST_mpi) assert(0);//endST -= lenST;
	if(endST_mpi==topST_mpi) assert(0);
	if(usedST_mpi==lenST_mpi) assert(0);
	if(topST_mpi!=0) assert(0);
	return spikeBufferFull;
}

inline float dvdtIzh(float v, float u, float tmpI, float h) {
	return (((0.04*v+5.0)*v+140.0-u+tmpI)*h);
}

// single integration step for recovery equation of 4-param Izhikevich
inline float dudtIzh(float v, float u, float a, float b, float h) {
 	return (a*(b*v-u)*h);
}

#define dvdtIzh_simd(v,u,tmpI,h) (((v0_04*(v)+v5_0)*(v)+v140_0-(u)+(tmpI))*(h))
#define dudtIzh_simd(v,u,a,b,h) ((a)*((b)*(v)-(u))*(h))
#if 1
static void neuronUpdate_simd(){
//st5=rpcc();

	int i=0,j,it;
	floatv4 viNMDA, vtmpI;
	floatv4 vgNMDA, vgGABAb;

	for(i=0;i<(swInfo.SizeN/4)*4;i+=4) {
		floatv4 vgAMPA,vgGABAa;

		vgNMDA =simd_set_floatv4(nInfo[i].gNMDA_d,nInfo[i+1].gNMDA_d,nInfo[i+2].gNMDA_d,nInfo[i+3].gNMDA_d);
		vgGABAb=simd_set_floatv4(nInfo[i].gGABAb_d,nInfo[i+1].gGABAb_d,nInfo[i+2].gGABAb_d,nInfo[i+3].gGABAb_d);
		vgAMPA =simd_set_floatv4(nInfo[i].gAMPA,nInfo[i+1].gAMPA,nInfo[i+2].gAMPA,nInfo[i+3].gAMPA);
		vgGABAa=simd_set_floatv4(nInfo[i].gGABAa,nInfo[i+1].gGABAa,nInfo[i+2].gGABAa,nInfo[i+3].gGABAa);

		floatv4 vv = simd_set_floatv4(nInfo[i].voltage,nInfo[i+1].voltage,nInfo[i+2].voltage,nInfo[i+3].voltage);
		floatv4 vu = simd_set_floatv4(nInfo[i].recovery,nInfo[i+1].recovery,nInfo[i+2].recovery,nInfo[i+3].recovery);
		floatv4 vh = swInfo.dt;
		floatv4 va = simd_set_floatv4(nInfo[i].Izh_a,nInfo[i+1].Izh_a,nInfo[i+2].Izh_a,nInfo[i+3].Izh_a);
		floatv4 vb = simd_set_floatv4(nInfo[i].Izh_b,nInfo[i+1].Izh_b,nInfo[i+2].Izh_b,nInfo[i+3].Izh_b);

		float tmpgNMDA_d[4],tmpgGABAb[4];
		float tmpgAMPA[4],tmpgGABAa[4];
		float tmpv[4],tmpu[4];

	for(it=0;it<swInfo.Ndt;it++){
		if (swInfo.sim_with_conductances) {
			floatv4 vvolt = -60.0;
			floatv4 aa = 80.0,bb=60.0,cc=70.0,dd=90.0,ee=0.0,ff=1.;
			viNMDA = (vvolt+aa)*(vvolt+aa)/(bb*bb);

			vtmpI=ee-(vgAMPA*(vvolt-ee)
				 +vgNMDA*viNMDA/(ff+viNMDA)*(vvolt-ee)
				 +vgGABAa*(vvolt+cc)
				 +vgGABAb*(vvolt+dd));
			

		} else {
			vtmpI=vgAMPA;
		}
		/* 4th Runge-Kutta */
		//float k1,k2,k3,k4;
		//float l1,l2,l3,l4;
		floatv4 v0_5=0.5,v0_04=0.04,v5_0=5.0,v140_0=140.0;
		floatv4 vk1 = dvdtIzh_simd(vv, vu, vtmpI, vh);
		floatv4 vl1 = dudtIzh_simd(vv, vu, va, vb, vh);

		floatv4 vk2 = dvdtIzh_simd(vv+v0_5*vk1, vu+v0_5*vl1, vtmpI, vh);
		floatv4 vl2 = dudtIzh_simd(vv+v0_5*vk1, vu+v0_5*vl1, va, vb,vh);

		floatv4 vk3 = dvdtIzh_simd(vv+v0_5*vk2, vu+v0_5*vl2, vtmpI, vh);
		floatv4 vl3 = dudtIzh_simd(vv+v0_5*vk2, vu+v0_5*vl2, va, vb,vh);

		floatv4 vk4 = dvdtIzh_simd(vv+vk3, vu+vl3, vtmpI, vh);
		floatv4 vl4 = dudtIzh_simd(vv+vk3, vu+vl3, va, vb, vh);

		floatv4 v2_0=2.0,v6_0=6.0;
		vv += (vk1+v2_0*vk2+v2_0*vk3+vk4)/v6_0;

		simd_store(vv,&(tmpv[0]));

		if (tmpv[0] > 30.0) tmpv[0] = 30.0;
		if (tmpv[0] < -90.0) tmpv[0] = -90.0;
		if (tmpv[1] > 30.0) tmpv[1] = 30.0;
		if (tmpv[1] < -90.0) tmpv[1] = -90.0;
		if (tmpv[2] > 30.0) tmpv[2] = 30.0;
		if (tmpv[2] < -90.0) tmpv[2] = -90.0;
		if (tmpv[3] > 30.0) tmpv[3] = 30.0;
		if (tmpv[3] < -90.0) tmpv[3] = -90.0;


		vu += (vl1+v2_0*vl2+v2_0*vl3+vl4)/v6_0;
		simd_store(vu,&(tmpu[0]));

		for(j=0;j<4;j++){
			if (tmpv[j]>= 30.0) {
				tmpv[j] = nInfo[i+j].Izh_c;
				tmpu[j]+= nInfo[i+j].Izh_d;
				if(addSpikeToTable_simd_mpi(i+j,it)) assert(0);//????
			}
		}

		simd_load(vv,&(tmpv[0]));
		simd_load(vu,&(tmpu[0]));

		int dIndex=offset+it;////????????
		assert(dIndex<lenRB);
		int addr = i*lenRB + dIndex;

		simd_store(vgAMPA,&(tmpgAMPA[0]));
		simd_store(vgNMDA,&(tmpgNMDA_d[0]));

		if (swInfo.sim_with_conductances) {
			for(j=0;j<4;j++){
				int addr2 = addr+j*lenRB;
				tmpgAMPA[j] += ringBuffer[addr2];
				tmpgNMDA_d[j] += ringBuffer[addr2];
				ringBuffer[addr2] = 0.;
			}
		} else {
			assert(0);
		}
		simd_load(vgAMPA,&(tmpgAMPA[0]));
		simd_load(vgNMDA,&(tmpgNMDA_d[0]));
	} //end Ndt

	/****simd store******/
		simd_store(vgGABAa,&(tmpgGABAa[0]));
		simd_store(vgGABAb,&(tmpgGABAb[0]));

		for(j=0;j<4;j++){
			nInfo[i+j].voltage=tmpv[j];
			nInfo[i+j].recovery=tmpu[j];
			nInfo[i+j].gAMPA=tmpgAMPA[j];
			nInfo[i+j].gNMDA_d=tmpgNMDA_d[j];
			nInfo[i+j].gGABAa=tmpgGABAa[j];
			nInfo[i+j].gGABAb_d=tmpgGABAb[j];
		}
	
	} // end SizeN
	int nr = swInfo.SizeN-(swInfo.SizeN/4)*4;
	if(nr>0){
		int i0=(swInfo.SizeN/4)*4;

		float tmpgNMDA_d[4],tmpgGABAb[4];
		float tmpgAMPA[4],tmpgGABAa[4];
		float tmpv[4],tmpu[4];
		float tmpa[4],tmpb[4];
		
		for(j=0;j<nr;j++){
			tmpgNMDA_d[j]=nInfo[i0+j].gNMDA_d;
			tmpgGABAb[j]=nInfo[i0+j].gGABAb_d;
			tmpgAMPA[j]=nInfo[i0+j].gAMPA;
			tmpgGABAa[j]=nInfo[i0+j].gGABAa;

			tmpv[j]=nInfo[i0+j].voltage;
			tmpu[j]=nInfo[i0+j].recovery;
			tmpa[j]=nInfo[i0+j].Izh_a;
			tmpb[j]=nInfo[i0+j].Izh_b;

		}

		floatv4 vgAMPA,vgGABAa;
		floatv4 vv,vu,va,vb;

		simd_load(vgNMDA,&(tmpgNMDA_d[0]));
		simd_load(vgGABAb,&(tmpgGABAb[0]));
		simd_load(vgAMPA,&(tmpgAMPA[0]));
		simd_load(vgGABAa,&(tmpgGABAa[0]));
		
		simd_load(vv,&(tmpv[0]));
		simd_load(vu,&(tmpu[0]));
		simd_load(va,&(tmpa[0]));
		simd_load(vb,&(tmpb[0]));

		floatv4 vh = swInfo.dt;

	for(it=0;it<swInfo.Ndt;it++){
		if (swInfo.sim_with_conductances) {
			floatv4 vvolt = -60.0;
			floatv4 aa = 80.0,bb=60.0,cc=70.0,dd=90.0,ee=0.0,ff=1.;
			viNMDA = (vvolt+aa)*(vvolt+aa)/(bb*bb);

			vtmpI=ee-(vgAMPA*(vvolt-ee)
				 +vgNMDA*viNMDA/(ff+viNMDA)*(vvolt-ee)
				 +vgGABAa*(vvolt+cc)
				 +vgGABAb*(vvolt+dd));
			

		} else {
			vtmpI=vgAMPA;
		}

		/* 4th Runge-Kutta */
		//float k1,k2,k3,k4;
		//float l1,l2,l3,l4;
		floatv4 v0_5=0.5,v0_04=0.04,v5_0=5.0,v140_0=140.0;
		floatv4 vk1 = dvdtIzh_simd(vv, vu, vtmpI, vh);
		floatv4 vl1 = dudtIzh_simd(vv, vu, va, vb, vh);

		floatv4 vk2 = dvdtIzh_simd(vv+v0_5*vk1, vu+v0_5*vl1, vtmpI, vh);
		floatv4 vl2 = dudtIzh_simd(vv+v0_5*vk1, vu+v0_5*vl1, va, vb,vh);

		floatv4 vk3 = dvdtIzh_simd(vv+v0_5*vk2, vu+v0_5*vl2, vtmpI, vh);
		floatv4 vl3 = dudtIzh_simd(vv+v0_5*vk2, vu+v0_5*vl2, va, vb,vh);

		floatv4 vk4 = dvdtIzh_simd(vv+vk3, vu+vl3, vtmpI, vh);
		floatv4 vl4 = dudtIzh_simd(vv+vk3, vu+vl3, va, vb, vh);

		floatv4 v2_0=2.0,v6_0=6.0;
		vv += (vk1+v2_0*vk2+v2_0*vk3+vk4)/v6_0;

		simd_store(vv,&(tmpv[0]));

		if (tmpv[0] > 30.0) tmpv[0] = 30.0;
		if (tmpv[0] < -90.0) tmpv[0] = -90.0;
		if (tmpv[1] > 30.0) tmpv[1] = 30.0;
		if (tmpv[1] < -90.0) tmpv[1] = -90.0;
		if (tmpv[2] > 30.0) tmpv[2] = 30.0;
		if (tmpv[2] < -90.0) tmpv[2] = -90.0;
		if (tmpv[3] > 30.0) tmpv[3] = 30.0;
		if (tmpv[3] < -90.0) tmpv[3] = -90.0;


		vu += (vl1+v2_0*vl2+v2_0*vl3+vl4)/v6_0;
		simd_store(vu,&(tmpu[0]));

		/*** findFiring ***/
		for(j=0;j<nr;j++){
			if (tmpv[j]>= 30.0) {
				tmpv[j] = nInfo[i0+j].Izh_c;
				tmpu[j]+= nInfo[i0+j].Izh_d;
				if(addSpikeToTable_simd_mpi(i0+j,it)) assert(0);//????
			}
		}

		simd_load(vv,&(tmpv[0]));
		simd_load(vu,&(tmpu[0]));

		/******* current set *******/
		int dIndex=offset+it;////????????
		assert(dIndex<lenRB);
		int addr = i0*lenRB + dIndex;

		simd_store(vgAMPA,&(tmpgAMPA[0]));
		simd_store(vgNMDA,&(tmpgNMDA_d[0]));

		if (swInfo.sim_with_conductances) {
			for(j=0;j<nr;j++){
				int addr2 = addr+j*lenRB;
				tmpgAMPA[j] += ringBuffer[addr2];
				tmpgNMDA_d[j] += ringBuffer[addr2];
				ringBuffer[addr2] = 0.;
			}
		} else {
			assert(0);
		}
		simd_load(vgAMPA,&(tmpgAMPA[0]));
		simd_load(vgNMDA,&(tmpgNMDA_d[0]));
	} 
	/****simd store******/
		simd_store(vgGABAa,&(tmpgGABAa[0]));
		simd_store(vgGABAb,&(tmpgGABAb[0]));

		for(j=0;j<nr;j++){
			nInfo[i0+j].voltage=tmpv[j];
			nInfo[i0+j].recovery=tmpu[j];

			nInfo[i0+j].gAMPA=tmpgAMPA[j];
			nInfo[i0+j].gNMDA_d=tmpgNMDA_d[j];
			nInfo[i0+j].gGABAa=tmpgGABAa[j];
			nInfo[i0+j].gGABAb_d=tmpgGABAb[j];
		}
	}



}
#endif
static int addSpikeToTable_simd(int i,int it) {
	int spikeBufferFull = 0;
	if(i<swInfo.SizeN) {nInfo[i].nSpikeCnt++;}
	firingTable[endST].nid = i+swInfo.StartN;
	firingTable[endST].time= it+sliceTime;//????????
	endST++; usedST++; numST[simTime%swInfo.Ndelay]++;
	if(endST>=lenST) endST -= lenST;
	if(endST==topST) assert(0);
	if(usedST==lenST) assert(0);
	swInfo.fireCnt++;
	return spikeBufferFull;
}
static int addSpikeToTable_simd_mpi(int i,int it) {
	int spikeBufferFull = 0;
	if(i<swInfo.SizeN) {nInfo[i].nSpikeCnt++;}//mpi++++
	firingTable_mpi[endST_mpi].nid = i+swInfo.StartN;
	if (firingTable_mpi[endST_mpi].nid==0xffff) firingTable_mpi[endST_mpi].nid = 0;
	firingTable_mpi[endST_mpi].time= it+sliceTime;//????????
	endST_mpi++; usedST_mpi++;
	if(endST_mpi>=lenST_mpi) assert(0);//endST -= lenST;
	if(endST_mpi==topST_mpi) assert(0);
	if(usedST_mpi==lenST_mpi) assert(0);
	if(topST_mpi!=0) assert(0);
	return spikeBufferFull;
}

//write firingTable_mpi[] from slave to host memory

static int SpikeDmaWrite_mpi(ptr){//mpi++++

	int i,core_number,all_number=0;
	//gather slaves spike numbers
	intv8 numCol = 0;
	intv8 numRow = 0;
	int col = (_MYID)%8;
	int row = (_MYID)/8;
	intv8 _v = 0;
	((int*)(&numCol))[col] = (int)endST_mpi;
	
	if(col==1||col==3||col==5||col==7){
		REG_PUTR(numCol,col-1);
	}
	if(col==0||col==2||col==4||col==6){
		REG_GETR(_v);numCol |= _v;
	}
	if(col==2||col==6){
		REG_PUTR(numCol,col-2);
	}
	if(col==0||col==4){
		REG_GETR(_v);numCol |= _v;
	}
	if(col==4){
		REG_PUTR(numCol,col-4);
	}
	if(col==0){
		REG_GETR(_v);numCol |= _v;
		for(i=0;i<8;i++)
			((int*)(&numRow))[row] += ((int*)(&numCol))[i];
	}

	if(row==1||row==3||row==5||row==7){
		REG_PUTC(numRow,row-1);
	}
	if(row==0||row==2||row==4||row==6){
		REG_GETC(_v);numRow |= _v;
	}
	if(row==2||row==6){
		REG_PUTC(numRow,row-2);
	}
	if(row==0||row==4){
		REG_GETC(_v);numRow |= _v;
	}
	if(row==4){
		REG_PUTC(numRow,row-4);
	}
	if(row==0){
		REG_GETC(_v);numRow |= _v;
	}	

	if (row==0){REG_PUTC(numRow,8);}
	else {REG_GETC(numRow);}		
	
	if (col==0){REG_PUTR(numCol,8);REG_PUTR(numRow,8);}
	else {REG_GETR(numCol);REG_GETR(numRow);}
    
	//dma write
	long addr = 0,reply = 0;
	for(i=0;i<col;i++) addr += ((int*)(&numCol))[i];
	for(i=0;i<row;i++) addr += ((int*)(&numRow))[i];

	if(endST_mpi){
		athread_put(PE_MODE,&(firingTable_mpi[0]),&(swInfo.firingTableHost[addr]),sizeof(spikeTime_t)*endST_mpi,&reply,0,0);
		dmaWait(&reply,1);
	}
	if(_MYID==0){
		NSgroup_slave = 0;
		for(i=0;i<8;i++) NSgroup_slave += ((int*)(&numRow))[i];
		reply = 0;
		athread_put(PE_MODE,&(NSgroup_slave),&(NS_group),sizeof(int),&reply,0,0);//发送脉冲数量
		dmaWait(&reply,1);
	}
	return 0;
}

//read firingTableAll_mpi[] from host memory to slave
static int SpikeDmaRead_mpi(ptr){//mpi++++
	//get data for firingTableAll[]
	if(_MYID==0){
		reply = 0;
		NSall_slave=0;
		
		numST_mpi[simTime%swInfo.Ndelay] = 0xfffffff;
		athread_get(PE_MODE,
			&(NSall),
			&(numST_mpi[simTime%swInfo.Ndelay]),
			sizeof(int),
			&reply,0,0,0);
		
		while(numST_mpi[simTime%swInfo.Ndelay]==0xfffffff);
		
	}
	
	intv8 _v;
	((int*)(&_v))[0] = numST_mpi[simTime%swInfo.Ndelay];
	_v = put_get_intv8(_v,0);
	numST_mpi[simTime%swInfo.Ndelay]=((int*)(&_v))[0];
	return 0;
}

static void InputCurrent(float wt, float nspike);
void SpikeDeliver(void *ptr){
	SpikeDmaRead_mpi(ptr);//mpi++++
	CurrentUpdate_mpi(ptr);
	float wt=0.00085;
	float nspike;
	nspike = currentfactor;
	InputCurrent(wt,nspike);
	return;
}
static void InputCurrent(float wt, float nspike)
{
	int dIndex = offset+swInfo.Ndt/2;
        if (dIndex>=lenRB)dIndex-=lenRB;
        float change = wt*nspike;
        int i;
	for(i=0;i<swInfo.SizeN;i++){
		ringBuffer[i*lenRB+dIndex] += change;
	}
}
static void syndma(spikeTime_t st,synInfo_t *sInfoLc);
static void syndma2(spikeTime_t st,synInfo_t *sInfoLc);
static void put_get_syn(synInfo_t *sInfoLc);
static void put_get_syn(synInfo_t *sInfoLc){

	intv8 *_v=(intv8*)sInfoLc;
	int len = swInfo.MaxN/4;
	int ilc = COL(_MYID);
	int i,j,offset;

	if(_MYID&0x01){
		intv8 *_vsrc = _v-len;
		int dst=ilc-1;
			REG_PUTR(_vsrc[0],dst);
			REG_PUTR(_vsrc[1],dst);
			REG_PUTR(_vsrc[2],dst);
			REG_PUTR(_vsrc[3],dst);
			REG_PUTR(_vsrc[4],dst);
			REG_GETR(_v[0]);
			REG_GETR(_v[1]);
			REG_GETR(_v[2]);
			REG_GETR(_v[3]);
			REG_GETR(_v[4]);

	}
	else{
		intv8 *_vsrc = _v+len;
		int dst=ilc+1;
			REG_GETR(_v[0]);
			REG_GETR(_v[1]);
			REG_GETR(_v[2]);
			REG_GETR(_v[3]);
			REG_GETR(_v[4]);
			REG_PUTR(_vsrc[0],dst);
			REG_PUTR(_vsrc[1],dst);
			REG_PUTR(_vsrc[2],dst);
			REG_PUTR(_vsrc[3],dst);
			REG_PUTR(_vsrc[4],dst);
	
	}
	return;
}

static void CurrentUpdate_mpi(void *ptr){//????
	int i,nid,srcId,iSpike,ibreak;
	intv8 _v[8];
	srcId = 0;
	int idelay, iblock,nblock;
	for(idelay=0; idelay<swInfo.Ndelay; idelay++) {//swInfo.Ndelay==1
		iSpike=topST;
		ibreak=0;
		int ivarray=0;
		while(1){
			int j,jindex=0;
			if(srcId==_MYID){
				int nvarray = lenST_mpi/64;//lenst_mpi=size+64
				int length = lenST_mpi/64*64;
				assert(nvarray>0);
				if (ivarray%nvarray==0){
					int least = numST_mpi[idelay]-ivarray*64;
					iSpike=0;
					endST_mpi = MIN(least,length);
					if(endST_mpi>0){
					reply = 0;
					athread_get(PE_MODE,
                        			&(swInfo.firingTableAll[idelay*swInfo.NN+ivarray*64]),
                        			&(firingTable_mpi[0]),
                        			sizeof(spikeTime_t)*endST_mpi,
                        			&reply,0,0,0);
					dmaWait(&reply,1);
					}
					if (endST_mpi<0){fprintf(stderr,"warn::more than number\n");}
				}
				for(j=0;j<8;j++){
				for(i=0;i<8;i++){
					if(iSpike==endST_mpi) {
					int ii;
					for(ii=i;ii<8;ii++)
        				((spikeTime_t*)(&_v[j]))[ii].nid=0xffff;
					break;
					}
        				((spikeTime_t*)(&_v[j]))[i]=firingTable_mpi[iSpike];
					iSpike++;
					assert(iSpike<=lenST_mpi);
				}
				if(jindex) break;
				}
			} else {
			}
			ivarray++;
			for(j=0;j<8;j++){			
				_v[j] = put_get_intv8(_v[j],srcId);
				if(((spikeTime_t*)(&_v[j]))[7].nid==0xffff) break;
			}
			spikeTime_t st,st0;
			//reply=0;
	//twice hiding access, unused

#if 1	//hiding access & reducing synchronization
long tm0,tm1,tm2,tm3;
				st=((spikeTime_t*)(&_v[0]))[0];
				if(st.nid!=0xffff){
				rpl[0]=0;
{tm0=rpcc();{
        			dma_set_reply(&dma_get_syn,&rpl[0]);
				synInfo_t *sInfoLc=&sInfo[0];
				syndma(st,sInfoLc);
}tm1=rpcc();}
cdma+=tm1-tm0;
				}
			
			for(i=1;i<64;i++){
				st0=st;
				st=((spikeTime_t*)(&_v[0]))[i];
				if(st.nid==0xffff) {ibreak=1;break;}
				j=i&(0x07);
				rpl[j]=0;
        			dma_set_reply(&dma_get_syn,&rpl[j]);
				synInfo_t *sInfoLc=&sInfo[j*swInfo.Ndma];

{tm0=rpcc();{
				syndma(st,sInfoLc);
				j=(i-1)&(0x07);
				dma_wait(&rpl[j],1);
}tm1=rpcc();}
				sInfoLc=&sInfo[(j)*swInfo.Ndma];
{tm2=rpcc();{
				generatePostSpike_simd(st0,sInfoLc);
}tm3=rpcc();}
cdma+=tm1-tm0;cspike+=tm3-tm2;
			}
				st=((spikeTime_t*)(&_v)[0])[i-1];
				if(st.nid!=0xffff) {
				j=(i-1)&(0x07);
{tm0=rpcc();{
				dma_wait(&rpl[j],1);
}tm1=rpcc();}
cdma+=tm1-tm0;
				synInfo_t *sInfoLc=&sInfo[j*swInfo.Ndma];
{tm2=rpcc();{
				generatePostSpike_simd(st,sInfoLc);
}tm3=rpcc();}
cspike+=tm3-tm2;
				}
#endif
			if(ibreak) break;
		}		
		
	}
	/****chenge topST****/
	short imv=simTime%swInfo.Ndelay;
	topST += numST[imv];
  	if(topST>=lenST)topST-=lenST;
	usedST -= numST[imv];
	numST[imv]=0;
}

#if 0
static void CurrentUpdate(void *ptr){//????
//st3=rpcc();
	int i,nid,srcId,iSpike,ibreak;
	intv8 _v[8];

	for(srcId=0; srcId<64; srcId++) {
		int nblock, iblock;
		
		
		iSpike=topST;
		ibreak=0;
		while(1){
			int j,jindex=0;
			if(srcId==_MYID){
				for(j=0;j<8;j++){
					for(i=0;i<8;i++){
						if(iSpike==endST) {
						int ii;
							for(ii=i;ii<8;ii++)
								((spikeTime_t*)(&_v[j]))[ii].nid=0xffff;
							break;
						}
						//if(iSpike==endST) break;
							((spikeTime_t*)(&_v[j]))[i]=firingTable[iSpike];
						iSpike++;
						if(iSpike>=lenST)iSpike -= lenST;
					}
					if(jindex) break;
				}
			}

			for(j=0;j<8;j++){			
				_v[j] = put_get_intv8(_v[j],srcId);
				if(((spikeTime_t*)(&_v[j]))[7].nid==0xffff) break;
			}


			spikeTime_t st,st0;
			//reply=0;
			
	//twice hiding access, unused

#if 1	//hiding access & reducing synchronization
long tm0,tm1,tm2,tm3;
				st=((spikeTime_t*)(&_v[0]))[0];
				if(st.nid!=0xffff){
				rpl[0]=0;
{tm0=rpcc();{
        		dma_set_reply(&dma_get_syn,&rpl[0]);
				synInfo_t *sInfoLc=&sInfo[0];
				syndma(st,sInfoLc);
}tm1=rpcc();}
cdma+=tm1-tm0;
				}
			
			
			for(i=1;i<64;i++){
				st0=st;
				st=((spikeTime_t*)(&_v[0]))[i];
				if(st.nid==0xffff) {ibreak=1;break;}
				j=i&(0x07);
				rpl[j]=0;
        		dma_set_reply(&dma_get_syn,&rpl[j]);
				synInfo_t *sInfoLc=&sInfo[j*swInfo.Ndma];
//long tm0,tm1,tm2,tm3;
{tm0=rpcc();{
				syndma(st,sInfoLc);
				j=(i-1)&(0x07);
				dma_wait(&rpl[j],1);
}tm1=rpcc();}
				sInfoLc=&sInfo[(j)*swInfo.Ndma];
{tm2=rpcc();{
				//generatePostSpike(st0,sInfoLc);
				generatePostSpike_simd(st0,sInfoLc);
}tm3=rpcc();}
cdma+=tm1-tm0;cspike+=tm3-tm2;
//if(_MYID==0)printf("timedma=%d timecmp=%d \n",tm1-tm0,tm3-tm2);

			}
				st=((spikeTime_t*)(&_v)[0])[i-1];
				if(st.nid!=0xffff) {
				j=(i-1)&(0x07);
				//while(reply<=i);
{tm0=rpcc();{
				dma_wait(&rpl[j],1);
}tm1=rpcc();}
cdma+=tm1-tm0;
				synInfo_t *sInfoLc=&sInfo[j*swInfo.Ndma];
{tm2=rpcc();{
				//generatePostSpike(st,sInfoLc);
				generatePostSpike_simd(st,sInfoLc);
}tm3=rpcc();}
cspike+=tm3-tm2;
				}
#endif

			if(ibreak) break;
		}		
		
	}	
	/****chenge topST****/
	short imv=simTime%swInfo.Ndelay;
	
	topST += numST[imv];
  	if(topST>=lenST)topST-=lenST;
	usedST -= numST[imv];
	numST[imv]=0;

//ed3=rpcc();
}
#endif
static void PoisCurrentUpdate(void *ptr){//????

}

static void generatePostSpike(spikeTime_t st,synInfo_t *sInfoLc) 
{
	int pre_i=st.nid;
	int is;
	int MaxN=swInfo.MaxN;
	int Ndma=swInfo.Ndma;
	//printf("%ddam \n",dma);
	int Ndelay=swInfo.Ndelay; 
	short occurTime=st.time>>swInfo.Nop;
	short iD = simTime-occurTime-1;
	if(!(iD>=0&&iD<Ndelay))iD = 0;
	for (is=0;is<Ndma;is++){
		unsigned int post_i = sInfoLc[is].postId;
		if(post_i==0xffff) break;
		unsigned int dIndex = st.time+sInfoLc[is].dl-sliceTime+offset;
		dIndex &= (lenRB-1);
		float change = sInfoLc[is].wt;
		int i = post_i&0xf;
		ringBuffer[i*lenRB+dIndex] += change;
	}
	return;
}
static void syndma(spikeTime_t st,synInfo_t *sInfoLc) {
	int pre_i=st.nid;
	if(!(pre_i<swInfo.preN&&pre_i>=0)) pre_i = 0;
	int is;
	int MaxN=swInfo.MaxN;
	int Ndma=swInfo.Ndma;
	int Ndelay=swInfo.Ndelay;

	short occurTime=st.time>>swInfo.Nop;
	short iD = simTime-occurTime-1;
	if(!(iD>=0&&iD<Ndelay))iD=0;
	unsigned int addr = pre_i*Ndelay*NTh*MaxN+iD*NTh*MaxN+_MYID*MaxN;
	dma(dma_get_syn,&swInfo.sInfoHost[addr],sInfoLc);
	return;
}
static void syndma2(spikeTime_t st,synInfo_t *sInfoLc) {
	int pre_i=st.nid;
	int is;
	int MaxN=swInfo.MaxN;
	int Ndma=swInfo.Ndma;
	int Ndelay=swInfo.Ndelay;

	short occurTime=st.time>>swInfo.Nop;
	short iD = simTime-occurTime-1;
	if(!(iD>=0&&iD<Ndelay))iD=0;
	unsigned int addr = pre_i*Ndelay*NTh*MaxN+iD*NTh*MaxN+(_MYID/2)*2*MaxN;
	dma(dma_get_syn,&swInfo.sInfoHost[addr],sInfoLc);
	return;
}

#if 1
static void generatePostSpike_simd(spikeTime_t st,synInfo_t *sInfoLc) {
	int pre_i=st.nid;
	int is,j;
	int MaxN=swInfo.MaxN;
	int Ndma=swInfo.Ndma;
	int Ndelay=swInfo.Ndelay;

	short occurTime=st.time>>swInfo.Nop;
	short iD = simTime-occurTime-1;
	if(!(iD>=0&&iD<Ndelay))iD=0;

	intv8 vtime=simd_set_intv8((int)st.time,(int)st.time,(int)st.time,(int)st.time,(int)st.time,(int)st.time,(int)st.time,(int)st.time);
	intv8 vsliceTime=(intv8)sliceTime;
	intv8 voffset=(intv8)offset;
	intv8 vStartN=(intv8)swInfo.StartN;
	intv8 vlenRB=(intv8)lenRB;
	for (is=0;is<(Ndma/8*8);is+=8){
		if(sInfoLc[is+7].postId==0xffff)break;
		intv8 vdl,vpid;
		vdl=simd_set_intv8((int)(sInfoLc[is].dl),(int)(sInfoLc[is+1].dl),(int)(sInfoLc[is+2].dl),(int)(sInfoLc[is+3].dl),(int)(sInfoLc[is+4].dl),(int)(sInfoLc[is+5].dl),(int)(sInfoLc[is+6].dl),(int)(sInfoLc[is+7].dl));
		vpid=simd_set_intv8((int)(sInfoLc[is].postId),(int)(sInfoLc[is+1].postId),(int)(sInfoLc[is+2].postId),(int)(sInfoLc[is+3].postId),(int)(sInfoLc[is+4].postId),(int)(sInfoLc[is+5].postId),(int)(sInfoLc[is+6].postId),(int)(sInfoLc[is+7].postId));
		intv8 vdIndex = vtime+vdl-vsliceTime+voffset;
		intv8 vlen = lenRB-1;
		vdIndex = vdIndex & vlen;
		intv8 vff = 0xf;
		intv8 vi = vpid&vff;
		intv8 vNop = (intv8)(swInfo.Nop+1);
		int op = swInfo.Nop+1;
		intv8 vaddr=vi<<4;
		vaddr += vdIndex;
		ringBuffer[((int*)(&vaddr))[0]] += sInfoLc[is].wt;
		ringBuffer[((int*)(&vaddr))[1]] += sInfoLc[is+1].wt;
		ringBuffer[((int*)(&vaddr))[2]] += sInfoLc[is+2].wt;
		ringBuffer[((int*)(&vaddr))[3]] += sInfoLc[is+3].wt;
		ringBuffer[((int*)(&vaddr))[4]] += sInfoLc[is+4].wt;
		ringBuffer[((int*)(&vaddr))[5]] += sInfoLc[is+5].wt;
		ringBuffer[((int*)(&vaddr))[6]] += sInfoLc[is+6].wt;
		ringBuffer[((int*)(&vaddr))[7]] += sInfoLc[is+7].wt;
	}
/*************residule syn***********************/
#if 1
	for (;is<Ndma;is++){
		unsigned int post_i = sInfoLc[is].postId;
		if(post_i==0xffff) break;
		int dIndex = st.time+sInfoLc[is].dl-sliceTime+offset;
		dIndex &= (lenRB-1);
		float change = sInfoLc[is].wt;
		int i = post_i&0xf;
		ringBuffer[i*lenRB+dIndex] += change;
	}
#endif

	return;
}
#endif
