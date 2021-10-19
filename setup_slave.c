#include <stdio.h>
#include "slave.h"
#include "setup.h"
#include <assert.h>



void wait_reply(volatile unsigned long *reply, int m) {
	while(*reply != m) {};
}

int slave_tk(ptr_t *ptr){

    int i,j;
	int Ndma[NTh];
	for(i=0;i<NTh;i++)
		Ndma[i]=0;
	if(_MYID==0){
		printf("dfd%d\n",ptr->connInfo_s[1].maxDelay);
	}
	uint8_t minD=1;
	uint8_t maxD=ptr->connInfo_s->maxDelay;
	int preN = ptr->snnInfo_s.preN;
	int MaxN = ptr->snnInfo_s.MaxN;
	int numNReg = ptr->snnInfo_s.numNReg;
	int Ndelay = ptr->snnInfo_s.Ndelay;
	int Nsyn = ptr->snnInfo_s.Nsyn;
	int gStart = ptr->snnInfo_s.gStart;
	int gSize = ptr->snnInfo_s.gSize;
	int Nop = ptr->snnInfo_s.Nop;
	int MIND=minD<<Nop;
	int MAXD=maxD<<Nop;
	int synWt;

	synInfo_t Host;

	for(i=0;i<preN;i++){
		int offset=i*Ndelay*NTh*MaxN;
		if(i==numNReg) 
			ptr->connInfo_s++;

		int iD=0,iTh=0,iN[20];
		// assert(Ndelay<=20);
		for(iD=0;iD<Ndelay;iD++)
			iN[iD]=0;

		int NNN; //select 1,2,3,4
		NNN = numNReg/Nsyn; //40
		
		// assert(NNN>0);

		for(j=i%NNN;j<numNReg;j+=NNN) {
			if(j<gStart||j>=gStart+gSize) continue;
			if(j>=ptr->snnInfo_s.swInfo[iTh].StartN+ptr->snnInfo_s.swInfo[iTh].SizeN){
				iTh++;
				for(iD=0;iD<Ndelay;iD++)
					iN[iD]=0;
			}
			if((0)&&i==j) continue;

			uint8_t dVal;
			for(;;){
				dVal=((i+j/NNN)%Ndelay+1)<<Nop;//for test
				// assert((dVal>=MIND)&&(dVal<MAXD));
				iD = (dVal-MIND)>>Nop;
				// assert(iD<Ndelay && iD>=0);
				break;
			}
			synWt=ptr->connInfo_s->initWt;
			// assert(synWt>=0.);
			int addr=offset+iD*NTh*MaxN+iTh*MaxN+iN[iD];
			int reply=0;
			athread_get(PE_MODE,&ptr->snnInfo_s.sInfoHost[addr],&Host,sizeof(Host),&reply,0,0);
			//athread_put(PE_MODE,&cdma,&dma[_MYID],sizeof(int),&reply,0,0);

			// ptr->snnInfo_s.sInfoHost[addr].postId=j;
			// ptr->snnInfo_s.sInfoHost[addr].wt=synWt;
			// ptr->snnInfo_s.sInfoHost[addr].dl=dVal;
			iN[iD]++;
			if(iN[iD]>Ndma[iTh])
				Ndma[iTh]=iN[iD];
		}
	}

    return 0;
}
