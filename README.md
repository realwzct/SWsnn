# SWsnn

SWsnn is a new type of high-efficiency parallel spiking neural network simulator based on the SW26010 architecture of TaihuLight supercomputer, which supports large-scale and high-precision brain simulation, and has good acceleration effects

## Hardware dependencies

The programs are based on SW26010 of the Sunway TaihuLight supercomputer.
SW26010 is a many-core accelerator, whose cores are independent of each other, and register-level communication (RLC) mechanism allows data communication between cores with low latency

## Usage

- Download and make source

```bash
    make
```

- Run SWsnn on TaihuLight supercomputer

```bash
    make run
```
  
- This is part of `Makefile`, the program operating parameters can be changed

   ` bsub -b -I -q q_sw_expr -n 56 -np 4 -cgsp 64 -host_stack 256 -share_size 4096 ./swtest 10000 200000 2 0.125 1000 5.5`

   The first parameter: The number of synapses connected to each neuron  
   The second parameter: Number of neurons  
   The third parameter: Max delay  
   The fourth parameter: Time step(The time step can only be the reciprocal of the power of two)  
   The fifth parameter: Simulation time  
   The sixth parameter: Input current parameters 
       
       
       
    
        
- See results in terminal
