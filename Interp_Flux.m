%% Function to find average flux value per panel

function [FLUX_AVG_PANEL] = Interp_Flux(FLUX_SOLARPILOT, N_PANELS)
FLUX_AVG=FLUX_SOLARPILOT(end,2:end);                                        %Average azimuthal flux from the SolarPilot
index= 1:1:length(FLUX_AVG{1,:});                                           %No. of point from SolarPilot
n_panel=1:1:N_PANELS;                                                       %No. of panels
x_panel=length(index)/length(n_panel).*n_panel;                                                     
FLUX_AVG_PANEL=spline(index,FLUX_AVG{1,:},x_panel);
