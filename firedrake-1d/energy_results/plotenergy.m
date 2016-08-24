clc
close all


T = dlmread('fv_100_16.txt');
figure 
set(gca,'fontsize', 14)
semilogy( T(:,1), T(:,2))
title('Plot of the drift from initial energy')
xlabel('Wave periods')
ylabel('Error in energy')