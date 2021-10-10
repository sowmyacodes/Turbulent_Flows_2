% Runs all the Holmedal & Myrhaug (2009) cases
addpath ~/MatRANS/src/;
cd Symmetric; unix('rm out_MatRANS.mat'); MatRANS; cd ..;
cd Asymmetric; unix('rm out_MatRANS.mat'); MatRANS; cd ..;

    
    
