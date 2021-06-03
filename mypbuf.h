#ifndef _MY_PBUF_H_
#define _MY_PBUF_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

typedef int spikegroupid_t ;
//! Type for specifying the delay in time steps
typedef unsigned short int delaystep_t;
//! The size of one allocation chunk size in PropagatedSpikeBuffer
#define PROPAGATED_SPIKE_BUFFER_CHUNK_SIZE 1024

typedef struct StgNode_s
{
	spikegroupid_t stg;
	delaystep_t delay;
	struct StgNode_s* next;
}StgNode;

typedef struct const_iterator
{
	StgNode *node;
}const_iterator;

typedef struct Pbuf_s
{
	int currIdx ;
	StgNode** ringBufferFront;
	StgNode** ringBufferBack;
	StgNode** chunkBuffer;
	int ringBufferFront_size;
	int ringBufferBack_size;
	int chunkBuffer_size;

	StgNode *currentFreeChunk;
	int  nextFreeSrgNodeIdx;
	size_t nextFreeChunkIdx;

	StgNode* recycledNodes;
	int chunkSize;
	int currT;
	double fillitup[32];
}Pbuf_t;
void pbufInit(Pbuf_t *pbuf,int minDelay,int maxDelay);
void pbufFree(Pbuf_t *pbuf);
void scheduleSpikeTargetGroup(Pbuf_t* pbuf,spikegroupid_t stg,delaystep_t delay);

void nextTimeStep(Pbuf_t* pbuf);
void reset(Pbuf_t* pbuf,int minDelay, int maxDelay);
inline size_t length(Pbuf_t* pbuf) { return pbuf->ringBufferFront_size; };

void init(Pbuf_t* pbuf,size_t maxDelaySteps);
StgNode *getFreeNode(Pbuf_t* pbuf);
const_iterator beginSpikeTargetGroups(Pbuf_t* pbuf,int stepOffset);
const_iterator endSpikeTargetGroups();

#endif /*MYPBUF_H_*/
