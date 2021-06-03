# SWsnn

SWsnn SWsnn is a new type of high-efficiency parallel pulse neural network simulator based on the SW26010 architecture of Taihuzhiguang supercomputer, named, which supports large-scale and high-precision brain simulation, and has good acceleration effects.

## Hardware dependencies

This program based on SW26010 of the Sunway TaihuLight supercomputer.
SW26010 is a many-core accelerator, whose cores are independent of each other, and register-level communication (RLC) mechanism allows data communication between cores with low latency

## Usage

* Download and make source

```bash
    make
```

* Grant shell file permissions 

```bash
    chmod +x run.sh
```

* Run SWsnn on  Taihuzhiguang supercomputer

```bash
    ./run.sh
```


* See results in listening terminal
