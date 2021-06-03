#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"
#include <simd.h>
#include <dma.h>

#include "my_slave.h"

int dmaWait(int* reply,int value)
{
	while(reply[0]!=value);	
    return 0;
}
