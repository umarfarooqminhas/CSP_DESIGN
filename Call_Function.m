
FLUX_SOLARPILOT = readtable("flux-table.csv");                              %Average azimuthal flux from Solar Pilot
[FLUX_AVG] = Interp_Flux(FLUX_SOLARPILOT, N_PANELS);                        %Interpolation to find 1 flux value per panel
D_RECEIVER=25;                                                            %8.15,
HEIGHT_RECEIVER=30;
TOWER_HEIGHT=250;
[Q_ABS_TOT,Q_INCIDENT_TOT,ETA_THERMAL,M_HTF_TOT] = Receiver_model(FLUX_AVG,D_RECEIVER,HEIGHT_RECEIVER,TOWER_HEIGHT)