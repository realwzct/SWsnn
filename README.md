# SWsnn

SWsnn SWsnn is a new type of high-efficiency parallel pulse neural network simulator based on the SW26010 architecture of Taihuzhiguang supercomputer, named, which supports large-scale and high-precision brain simulation, and has good acceleration effects.

## Hardware dependencies

This program based on SW26010 of the Sunway TaihuLight supercomputer.
SW26010 is a many-core accelerator, whose cores are independent of each other, and register-level communication (RLC) mechanism allows data communication between cores with low latency

## Usage

* Download and make source on each node

```bash
    make
```

* Start SunwayMR resource manager , while specifying `master IP`, `master port`, `shared threads of node`, `shared memory of node`

```bash
    ./sunwaymr -t 192.168.1.85 19113 4 1024
```

* Now, you can run example program **on master** (in a new terminal)

```bash
    ./sunwaymr -a examples/SunwayMRPi.cpp
```

* See results in listening terminal
