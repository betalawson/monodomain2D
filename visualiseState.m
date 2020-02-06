function visualiseState( V, Nx, Ny, occ_map, t )
% This function visualises the membrane potential map

% Convert vector V back into a matrix (use transpose for correct
% orientation because reshape works 'columnwise')
V = reshape(V, Nx+1, Ny+1)';

% Calculate *element* Voltage values as the average of node corners
% Some will be NaNs, but we are happy for these to propagate because they
% will only 'infect' sites that are occupied anyway
Ve = ( V(1:Ny,1:Nx) + V(2:Ny+1,1:Nx) + V(1:Ny,2:Nx+1) + V(2:Ny+1,2:Nx+1) ) / 4;

% Make sure all occupied sites won't show up in the plot
Ve(occ_map) = NaN;

load('extra_colormaps.mat','plasma');

Vclr = [ [0.2, 0.2, 0.2]; plasma]; 

% Visualise the remainder
imagesc(flipud(Ve));
colormap(Vclr);
whitebg([0.2 0.2 0.2]);
caxis([-90 40]);
axis equal;


% Add title
title(['Membrane potential map after t = ',num2str(t),' ms'],'FontSize', 24);

end

