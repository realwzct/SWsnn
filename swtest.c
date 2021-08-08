#include "mysnn.h"
#include <sys/time.h>    // for gettimeofday() chenged!!
#include <stdlib.h>
int g_timestep;
float currentfactor;
int rank, nproc;

int main(int argc, char *argv[])
{

	int connectnumber = atoi(argv[1]);  			//connect number 1 0000
	int neuronnumber = atoi(argv[2]);				//neuron number 20 0000
	int delay = atoi(argv[3]);          			//max delay
	g_timestep = (int)(1/atof(argv[4]));              //time step
	int simulatetime = atoi(argv[5]);				//simulatetime
	currentfactor = atof(argv[6]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/*神经元种类连接种类*/
	int numGrp=3, numConn=2;
	int randSeed = 42;
	int NE = neuronnumber; //neuron number 20 0000
	int NI=0; int NP=0;

	snnInfo_t snnInfo;//结构体
	grpInfo_t *grpInfo=(grpInfo_t*)malloc(numGrp*sizeof(grpInfo_t));
	connInfo_t *connInfo=(connInfo_t*)malloc(numConn*sizeof(connInfo_t));

	long time0,time1,time2;
	time0=rpcc();

	//initNetwork(&snnInfo,grpInfo,connInfo,numGrp,numConn,randSeed);//初始化网络，但这段没啥用

	NeuronParameter(grpInfo,NE,NI,NP,&snnInfo,delay,connectnumber);//神经元和snn参数设定

	ConnectParameter(connInfo,delay,grpInfo);//连接参数设定

	setConductances(&snnInfo); //电导模型还是电流模型

	buildModel(&snnInfo,grpInfo,connInfo);

	createNetwork(&snnInfo,grpInfo,connInfo);
	MPI_Barrier(MPI_COMM_WORLD);//mpi syn
	
	time1=rpcc();

	runNetwork(&snnInfo,grpInfo,connInfo,snnInfo.nInfoHost,snnInfo.sInfoHost,snnInfo.swInfo,simulatetime,0);
	time2=rpcc();

	if(rank==0){
		printf("running time: %.4lf ms\n", (double)(time2-time1)*1000/CLOCKRATE);
		printf("all time: %.4lf ms\n", (double)(time1-time0)*1000/CLOCKRATE);
	}

    MPI_Finalize();
	return 0;
}
