bsub -b -I -q q_sw_expr -n 56 -np 4 -cgsp 64 -host_stack 256 -share_size 4096 ./swtest 10000 200000 2 0.125 1000 5.5
#./swtest connectnumber neuronnumber maxdelay timestep simluatetime currentfactor