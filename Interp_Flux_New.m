%% Function to find average flux value per panel

function [FLUX_AVG_PANEL] = Interp_Flux_New(FLUX_SOLARPILOT_AVG, N_PANELS)
index= 1:1:length(FLUX_SOLARPILOT_AVG);                                           %No. of point from SolarPilot
n_panel=1:1:N_PANELS;                                                       %No. of panels
x_panel=length(index)/length(n_panel).*n_panel;                                                     
FLUX_AVG_PANEL=spline(index,FLUX_SOLARPILOT_AVG,x_panel);
