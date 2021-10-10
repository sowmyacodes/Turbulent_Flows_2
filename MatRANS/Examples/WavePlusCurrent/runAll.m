% Run all wave plus current simulations
addpath ~/MatRANS/src;
cd CurrentOnly; unix('rm out_MatRANS.mat'); MatRANS; cd ..;
cd WaveOnly; unix('rm out_MatRANS.mat'); MatRANS; cd ..;
cd Combined; unix('rm out_MatRANS.mat'); MatRANS; cd ..;
