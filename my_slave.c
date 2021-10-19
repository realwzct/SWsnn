#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"
#include <simd.h>
#include <dma.h>
#include "my_slave.h"

intv8 put_get_intv8(intv8 _v,int srcId)
{
    int dst=8;
    int col=COL(srcId),row=ROW(srcId);
    if(COL(_MYID)==col&&ROW(_MYID)==row){
	REG_PUTR(_v,dst);
	REG_PUTC(_v,dst);
    }
    else if(COL(_MYID)!=col&&ROW(_MYID)==row){
        REG_GETR(_v);
	REG_PUTC(_v,dst);
    }
    else if(COL(_MYID)!=col&&ROW(_MYID)!=row){
        REG_GETC(_v);
    }
    else if(COL(_MYID)==col&&ROW(_MYID)!=row){
        REG_GETC(_v);
    }
    return _v;
}

