#ifndef _SETUP_H_
#define _SETUP_H_

#include "mysnn.h"
#include "swstruct.h"

typedef struct para{
    synInfo_t *sInfoHost;
    long len;
}para_t;
#define NTh  64

typedef struct ptr_s
{
	snnInfo_t snnInfo_s;
    connInfo_t *connInfo_s;
}ptr_t;



#endif
