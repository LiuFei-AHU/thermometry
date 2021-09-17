% Simulation Results: Permutations
load incoherent % nSlices x permutations

figure, histogram(rmsError(2,:),'BinEdges',[0 0.5 1 1.5 2 2.5 3]) 
% Histogram for 2 slices
set(gca,'fontsize',24), xticks(0:0.5:3), axis([0 3 0 50])
s = ['RMS Error (' char(176) 'C)'];
ylabel('# of permutations'), xlabel(s)

figure, histogram(rmsError(3,:),'BinEdges',[0 0.5 1 1.5 2 2.5 3])
% Histogram for 3 slices
set(gca,'fontsize',24), xticks(0:0.5:3), yticks(0:10:50), axis([0 3 0 50])
ylabel('# of permutations'), xlabel(s)