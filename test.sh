bsub -I -b -q q_sw_share -n 65 -np 1 -cgsp 64 -host_stack 256 -share_size 7000 ./pm /home/export/online1/systest/swsdu/xx/photon/lux/ O 1
bsub -I -b -q q_sw_share -n 65 -np 1 -cgsp 64 -host_stack 256 -share_size 7000 ./pm /home/export/online1/systest/swsdu/xx/photon/water/ O 1
bsub -I -b -q q_sw_share -n 65 -np 1 -cgsp 64 -host_stack 256 -share_size 7000 ./pm /home/export/online1/systest/swsdu/xx/photon/maxtorus/ O 1
