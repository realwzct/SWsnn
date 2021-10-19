
#include "mypbuf.h" 

void pbufInit(Pbuf_t *pbuf,int minDelay,int maxDelay)
{
	pbuf->chunkSize = PROPAGATED_SPIKE_BUFFER_CHUNK_SIZE;

        pbuf->currIdx=0;
        pbuf->ringBufferFront_size=maxDelay+1;
        pbuf->ringBufferBack_size=maxDelay+1;
        pbuf->chunkBuffer_size=0;
        pbuf->currentFreeChunk=NULL;
        pbuf->nextFreeSrgNodeIdx=0;
        pbuf->nextFreeChunkIdx=0;
        pbuf->recycledNodes=NULL;
        pbuf->ringBufferFront=(StgNode**)malloc((maxDelay+1)*sizeof(StgNode*));
        pbuf->ringBufferBack=(StgNode**)malloc((maxDelay+1)*sizeof(StgNode*));
        reset(pbuf,minDelay,maxDelay);
        pbuf->currT = 0;
}

void pbufFree(Pbuf_t *pbuf)
{
    size_t i;
    for(i=0; i<pbuf->chunkBuffer_size; i++) {
        free((void*)pbuf->chunkBuffer[i]);
    }
}

void init(Pbuf_t* pbuf,size_t maxDelaySteps)
{
    if( pbuf->ringBufferFront_size!=maxDelaySteps+1){
	printf("%s::%d\n",__FILE__,__LINE__);
        fflush(NULL);exit(0);
        pbuf->ringBufferFront=(StgNode**)malloc((maxDelaySteps+1)*sizeof(StgNode*));
        pbuf->ringBufferBack=(StgNode**)malloc((maxDelaySteps+1)*sizeof(StgNode*));
	pbuf->ringBufferFront_size=maxDelaySteps+1;
	pbuf->ringBufferBack_size=maxDelaySteps+1;
    }
    if( pbuf->chunkBuffer_size < 1 ) {
        pbuf->chunkBuffer=(StgNode**)malloc(10*sizeof(StgNode*));
        pbuf->chunkBuffer_size=0;
        pbuf->chunkBuffer[pbuf->chunkBuffer_size]
		=(StgNode*)malloc(pbuf->chunkSize*sizeof(StgNode));
	pbuf->chunkBuffer_size++;
    }
}

void reset(Pbuf_t* pbuf,int minDelay,int maxDelay)
{
    size_t i;
    init(pbuf,maxDelay+minDelay);
    for(i=0;i<pbuf->ringBufferFront_size;i++){
        pbuf->ringBufferFront[i]= NULL;
	pbuf->ringBufferBack[i] = NULL;
    }

    pbuf->currentFreeChunk = pbuf->chunkBuffer[0];
    pbuf->nextFreeChunkIdx = 1;

    pbuf->nextFreeSrgNodeIdx  = 0;

    pbuf->recycledNodes = NULL;

    pbuf->currIdx = 0;
}

const_iterator beginSpikeTargetGroups(Pbuf_t* pbuf,int stepOffset)
{
	const_iterator tmp;
	size_t size;
	size=length(pbuf);
  	tmp.node=pbuf->ringBufferFront[(pbuf->currIdx+stepOffset+size)%size];
        return tmp;
}

    //! End iterator corresponding to beginSynapseGroups
const_iterator endSpikeTargetGroups()
{
	const_iterator tmp;
  	tmp.node=NULL;
        return tmp;
};

StgNode* getFreeNode(Pbuf_t *pbuf)
{
    StgNode *n;
    if (pbuf->recycledNodes != NULL) {
        n = pbuf->recycledNodes;
        pbuf->recycledNodes = pbuf->recycledNodes->next;
    } else if ( pbuf->nextFreeSrgNodeIdx < pbuf->chunkSize ) {
        n = &(pbuf->currentFreeChunk[pbuf->nextFreeSrgNodeIdx++]);
    } else if (pbuf->nextFreeChunkIdx < pbuf->chunkBuffer_size ) {
        pbuf->currentFreeChunk = pbuf->chunkBuffer[pbuf->nextFreeChunkIdx++];
        n = &(pbuf->currentFreeChunk[0]);
        pbuf->nextFreeSrgNodeIdx = 1;
    } else {
        pbuf->currentFreeChunk=(StgNode*)malloc(pbuf->chunkSize*sizeof(StgNode));
        pbuf->chunkBuffer[pbuf->chunkBuffer_size]=pbuf->currentFreeChunk;
	pbuf->chunkBuffer_size++;
	if(pbuf->chunkBuffer_size>10){
		printf("%s::%d\n",__FILE__,__LINE__);
        fflush(NULL);exit(0);
	}
        pbuf->nextFreeChunkIdx++;
        n = &(pbuf->currentFreeChunk[0]);
        pbuf->nextFreeSrgNodeIdx = 1;
    }
    return n;
}

void scheduleSpikeTargetGroup(Pbuf_t *pbuf,spikegroupid_t stg,delaystep_t delay)
{
    StgNode *n = getFreeNode(pbuf);
    int writeIdx=(pbuf->currIdx+delay)%pbuf->ringBufferFront_size;
    n->stg    = stg;
    n->delay  = delay;
    n->next   = NULL;
    if( pbuf->ringBufferFront[writeIdx]==NULL){
        pbuf->ringBufferBack[writeIdx]=n;
        pbuf->ringBufferFront[writeIdx]=n;
    } else {
        pbuf->ringBufferBack[writeIdx]->next=n;
        pbuf->ringBufferBack[writeIdx] = n;
    } 
}

void nextTimeStep(Pbuf_t *pbuf)
{
    if(pbuf->ringBufferFront[pbuf->currIdx]!=NULL)
    {
        pbuf->ringBufferBack[pbuf->currIdx]->next = pbuf->recycledNodes;
        pbuf->recycledNodes = pbuf->ringBufferFront[pbuf->currIdx];
    }
    pbuf->ringBufferBack[pbuf->currIdx] = NULL;
    pbuf->ringBufferFront[pbuf->currIdx] = NULL;
    pbuf->currIdx=(pbuf->currIdx+1)%pbuf->ringBufferFront_size;
    pbuf->currT++;
}


