% Simulation Results: Hotspot Leakage
load Leakage % nSlices x permutations x hotspot size
rmsError = sqz(mean(rmsError,2)); % average RMS error over 50 permutations

figure, p = plot(8:8:32,rmsError(2:3,:),'LineWidth',4);
axis([-inf inf 0 1]), set(gca, 'fontsize',24), xticks([8 16 24 32])
s = ['RMS Error (' char(176) 'C)']; ylabel(s), xlabel('FWHM (Voxels)')
legend('2 Slices', '3 Slices')
legend('boxoff')

load leakageData
minDisplayTemp = 0; maxDisplayTemp = max(tempkacc(:));
for nSlices = 2:3
    figure
    for hotSize = 1:4
        for ii = 1:nSlices
            temp = tempkacc(:,:,ii,hotSize,nSlices-1);
            subplot(5,nSlices,(hotSize-1)*nSlices+ii)
            imagesc(temp,[minDisplayTemp,maxDisplayTemp]); axis image; axis off
            %h = colorbar; ylabel(h,'degrees C');
            if hotSize == 1
                title(sprintf('slice %d',ii));
            end
        end
    end
end
drawnow;
