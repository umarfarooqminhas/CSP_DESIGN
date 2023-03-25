%% Code for computation of Thermal Effeciency of Receiver
function [Q_ABS_TOT,Q_INCIDENT_TOT,ETA_THERMAL,M_HTF_TOT] = Receiver_model(FLUX_AVG,D_RECEIVER,HEIGHT_RECEIVER,TOWER_HEIGHT)

%% General Parameters

T_HTF_OUT = 550;                                                            %Heat Transfer Fluid Outlet Temperature in °C
T_HTF_IN = 260;                                                             %Heat Transfer Fluid Inlet Temperature in °C
WIND_VELOCITY_AVG = 5.6;                                                    %Average Wind Velocity at Receiver Height in m/s -- Wind Atlas PRIMM NV
P_AMB = 1;                                                                  %Ambient Pressure
ETA_PUMP = 0.65;                                                            %Pump Effeciency
T_AMB = 27;                                                                 %Ambient Temp in °C
T_SKY= 0;                                                                   %Sky Temperature
Q_HTF = 1000;                                                                %Coolent Mass Flowrate in kg/s
SIGMA = 5.678E-8;                                                           %Stefan-Boltzman Constant W/m^2-K^4
EPSILON = 0.88;                                                             %Emissivity

%% Receiver Dimensions

N_PANELS = 16;                                                              %Number of Vertical Panel in Receiver in m
%D_RECEIVER = 8.15;                                                          %Receiver overall Diameter in m
%HEIGHT_RECEIVER = 6.2;                                                      %Paner Height in m
%TOWER_HEIGHT = 120;                                                         %Tower Height from Ground in m
D_OUT_TUBE = 25*10^-3;                                                      %TUBE OUTER DIAMETER in m
THICKNESS_TUBE = 1.25*10^-3;                                                %Tube Thickness in m
D_IN_TUBE = D_OUT_TUBE-2*THICKNESS_TUBE;                                    %Tube Inner Diamter in m
A_TUBE_OUTER = pi*D_OUT_TUBE/2*HEIGHT_RECEIVER;                             %Tube Outer Area
A_RECEIVER = pi*D_RECEIVER*HEIGHT_RECEIVER;                                 %Receiver Outer Diameter
N_TUBE_PAN = floor(pi*D_RECEIVER/(D_OUT_TUBE*N_PANELS));                    %Tubes per panel
N_TUBES_TOT= N_PANELS*N_TUBE_PAN;                                           %Total Number of Tubes
A_PROJ_REC = D_OUT_TUBE*HEIGHT_RECEIVER*N_TUBES_TOT;                        %Projected Area of all Tubes
A_NODE = pi* D_RECEIVER* HEIGHT_RECEIVER/N_PANELS;                          %Each Panel considered as one node

%% FLUX DATA ACQUISITION

%FLUX_SOLARPILOT = readtable("flux-table.csv");                              %Average azimuthal flux from Solar Pilot
%[FLUX_AVG] = Interp_Flux(FLUX_SOLARPILOT, N_PANELS);                        %Interpolation to find 1 flux value per panel

for j=1:N_PANELS
    if j<=N_PANELS/2
    FLUX_AVG_PANEL(j)=FLUX_AVG(j+N_PANELS/2);
    else
    FLUX_AVG_PANEL(j)=FLUX_AVG(j-N_PANELS/2);
    end
end

%% FLOW PATTERN
N_FLOW_PATH= 2;                                                             %Flow paths inside the receiver
FLOW_DIRECTION_1= [1,2,3,4,12,11,10,9];                                     %Crossflow from North to South anticlocwise
FLOW_DIRECTION_2= [16,15,14,13,5,6,7,8];                                    %Crossflow from North to South clockwise                                  

%% INITILIZATION W/ FIRST GUESS VALUE

for j=1:N_PANELS
    T_SURFACE_PANEL(j) = 300;                                               %EACH PANEL SURFACE TEMPERATURE GUESS
    T_PANEL_IN(j) = 600;                                                    %EACH PANEL HTF TEMPERATURE INLET
    T_PANEL_OUT (j) =  600;                                                 %EACH PANEL HTF TEMPERATURE OULET
end

M_HTF=1500/N_FLOW_PATH;                                                      %HTF MASS FLOWRATE GUESS VALUE
T_HTF_HOT=1500;                                                             %HTF OUTLET TEMP GUESS VALUE
T_HTF_AVG_OUT=600;
%% CALCULATIONS
error=1;                                                                    %INITIALIZING ERROR VARIABLE
while error > 10^-4                                                         %CHECK CONVERGENCE
for j=1:N_PANELS
    T_PANEL_AVE(j)=(T_PANEL_OUT(j)+T_PANEL_IN(j))/2;                              %PANEL AVG TEMP FOR T WALL & AVG NEW SURFACE TEMP
    T_FILM(j)=(T_SURFACE_PANEL(j)+T_AMB)/2
end

T_HTF_REF=(T_HTF_AVG_OUT+T_HTF_IN)/2                                      %HTF AVERAGE PROPERTIES EVALUATION
T_SURFACE_AVG = sum(T_SURFACE_PANEL)/length(T_SURFACE_PANEL);               %RECEIVER SURFACE AVG TEMP For Forced Convection calculation
T_FILM_AVG = (T_AMB+T_HTF_OUT)/2;                                           %AVERGAE Film Temp for forced convection co-effient calculation

%% FORCED CONVECTION COEFFICIENT CALCULTION
K_AIR_FILM=5E-12*T_FILM_AVG^3 - 3E-08*T_FILM_AVG^2 + 8E-05*T_FILM_AVG^1 + 0.0236; %Conductivity of Air at FILM TEMP
MU_AIR_FILM= 3E-14*T_FILM_AVG^3 - 4E-11*T_FILM_AVG^2 + 5E-08*T_FILM_AVG^1 + 2E-05;                                                               %Viscosity of Air at FILM TEMP
RHO_AIR_FILM=4E-11*T_FILM_AVG^4 - 5E-08*T_FILM_AVG^3 + 2E-05*T_FILM_AVG^2 - 0.0054*T_FILM_AVG^1 + 1.2947;  %Density of Air at FILM TEMP
CP_AIR_FILM=-6E-07*T_FILM_AVG^3 + 0.0007*T_FILM_AVG^2 - 0.0456*T_FILM_AVG^1 + 1007.1;                                                               %Air Specific Heat
RE_FORCED=RHO_AIR_FILM*WIND_VELOCITY_AVG*D_RECEIVER/MU_AIR_FILM;                        %Reynold Number
KS_D=(D_OUT_TUBE/2)/D_RECEIVER;                                             %Relative roughness
%NUS_FOR=2.57*10^-3*RE_FORCED^0.98;                                           %Nusselt Correlation 
NUS_FOR=0.3+0.488*RE_FORCED^0.5*(1+(RE_FORCED/28200)^0.625)^0.8;
H_FOR_CONV=NUS_FOR*K_AIR_FILM/D_RECEIVER;
%% NATURAL CONVETION COEFFICIENT CALCULTION
BETA = 1/(T_AMB+273.15);                                                             %Volumetric Coefficient
NU_AMB= 8E-11*T_AMB^2 + 9E-08*T_AMB^1 + 1E-05;                              %Kinematic Viscosity at T am
G=9.8;
for j=1:N_PANELS
    GR_NAT(j) = G*BETA*((T_SURFACE_PANEL(j)+273.15)-(T_AMB+273.15))*HEIGHT_RECEIVER^3/NU_AMB^2; %Grashof Number
    NUS_NAT(j)=0.98*GR_NAT(j)^(1/3)*((T_SURFACE_PANEL(j)+273.15)/(T_AMB+273.15))^(-0.14);      %Nusselt Number
    H_NAT_CONV(j)=NUS_NAT(j)*K_AIR_FILM/HEIGHT_RECEIVER;                              %Natural Convection Coefficient
end
%% Mixed Convection 
M=3.2;
for j=1:N_PANELS
    H_MIXED(j)=((H_FOR_CONV^(M)+H_NAT_CONV(j)^(M)))^(1/M);                     %Mixed Convection Coefficient with overall forced coefficient and local Natural Convection coeffient
    Q_CONV_DOT(j)=H_MIXED(j)*A_NODE*((T_SURFACE_PANEL(j)+273.15)-(T_FILM(j)+273.15));         %CONVECTION LOSS OF EACH PANEL
end
%% RADIATION LOSSES CALCULATION
for j=1:N_PANELS
    H_RAD_AMB(j)=SIGMA*EPSILON*((T_SURFACE_PANEL(j)+273.15)^2+(T_AMB+273.15)^2)*((T_SURFACE_PANEL(j)+273.15)+(T_AMB+273.15));
    H_RAD_SKY(j)=SIGMA*EPSILON*((T_SURFACE_PANEL(j)+273.15)^2+(T_SKY^2+273.15))*((T_SURFACE_PANEL(j)+273.15)+(T_SKY+273.15));
    Q_RAD_AMB(j)=1/2*H_RAD_AMB(j)*A_NODE*((T_SURFACE_PANEL(j)+273.15)-(T_AMB+273.15));
    Q_RAD_SKY(j)=1/2*H_RAD_SKY(j)*A_NODE*((T_SURFACE_PANEL(j)+273.15)-(T_SKY+273.15));
end
%% TOTAL SURFACE LOSSES
for j=1:N_PANELS
    Q_RAD_TOT(j)=Q_RAD_SKY(j)+Q_RAD_AMB(j);                                 %Total Radiative Losses per panel (SKY + GROUND)
    Q_LOSS_TOT(j)=Q_CONV_DOT(j)+ Q_RAD_TOT(j);                              %Total Surface Losses per panel (Radiative + Convective)
end

Q_RAD_TOT_SUM = sum(Q_RAD_TOT);                                             %Total Radiative Losses
Q_COV_TOT_SUM = sum(Q_CONV_DOT);                                            %Total Convective Losses
Q_LOSS_TOT_SUM = sum(Q_LOSS_TOT);                                           %Total Surface Losses
%%  NET FLUX ON PANEL
for j=1:N_PANELS
    Q_INCIDENT_PANEL(j)=A_NODE*FLUX_AVG_PANEL(j)*1000;
    Q_ABS_PANEL(j)=Q_INCIDENT_PANEL(j)-Q_LOSS_TOT(j);
end
Q_INCIDENT_TOT= sum(Q_INCIDENT_PANEL);
Q_ABS_TOT = sum(Q_ABS_PANEL);

%% CONDUCTION TUBE WALL

for j=1:N_PANELS
    T_WALL(j)=(T_SURFACE_PANEL(j)+T_PANEL_AVE(j))/2;
    K_TUBE(j)=17;                                                         
    R_TUBE(j)=THICKNESS_TUBE/(K_TUBE(j)*HEIGHT_RECEIVER*D_RECEIVER*pi);     %correction Required?/
end
%% HTF CONVECTIVE HEAT TRANSFER COEFFICIENT CALCS
[CP_HTF,MU_HTF,K_HTF,RHO_HTF] = Molten_Salt_Props(T_HTF_REF);
U_HTF=M_HTF/(N_TUBE_PAN*pi*(D_IN_TUBE/2)^2*RHO_HTF);                        %HTF Velocity
RE_INT_TUBE=RHO_HTF*U_HTF*D_IN_TUBE/MU_HTF;                                 %Renold Number of Tube in flow
F=(0.79*log(RE_INT_TUBE)-1.64)^-2;                                           %Friction Factor for Re (3000 - 5E6)
PR_INT_TUBE=CP_HTF*MU_HTF/K_HTF;
NU_INT=(F/8)*(RE_INT_TUBE-1000)*PR_INT_TUBE/(1+12.7*(F/8)^(0.5)*(PR_INT_TUBE^(2/3)-1));
H_INN_HTF=NU_INT*K_HTF/D_IN_TUBE;
R_CONV_HTF=1/(H_INN_HTF*pi*D_IN_TUBE/2*HEIGHT_RECEIVER*N_TUBE_PAN);         %Convective Heat Transfer Resistance HTF
%% FLOW GRID 

for j=1:4
    if j==1  || j==N_PANELS
        T_PANEL_IN(j)=T_HTF_IN;
    elseif j<=N_PANELS/2 && j~=1 && j~=5
        T_PANEL_IN(j)=T_PANEL_OUT(j-1);
    elseif  j<N_PANELS && j>3*N_PANELS/4
        T_PANEL_IN(j)=T_PANEL_OUT(j+1);
    elseif j==5 
        T_PANEL_IN(j)=T_PANEL_OUT(j+N_PANELS/2);
     elseif j==3*N_PANELS/4
          T_PANEL_IN(j)=T_PANEL_OUT(N_PANELS/4);
    elseif j<3*N_PANELS/4 && j>N_PANELS/2
          T_PANEL_IN(j)=T_PANEL_OUT(j+1);
    end
  T_PANEL_OUT(j)=T_PANEL_IN(j)+Q_ABS_PANEL(j)/(M_HTF*CP_HTF);
    T_SURFACE_PANEL(j)=T_PANEL_AVE(j)+Q_ABS_PANEL(j)*(R_CONV_HTF+R_TUBE(j));

end

for j2=1:4
    j=17-j2;
    if j==1  || j==N_PANELS
        T_PANEL_IN(j)=T_HTF_IN;
    elseif j<=N_PANELS/2 && j~=1 && j~=5
        T_PANEL_IN(j)=T_PANEL_OUT(j-1);
    elseif  j<N_PANELS && j>3*N_PANELS/4
        T_PANEL_IN(j)=T_PANEL_OUT(j+1);
    elseif j==5 
        T_PANEL_IN(j)=T_PANEL_OUT(j+N_PANELS/2);
     elseif j==3*N_PANELS/4
          T_PANEL_IN(j)=T_PANEL_OUT(N_PANELS/4);
    elseif j<3*N_PANELS/4 && j>N_PANELS/2
          T_PANEL_IN(j)=T_PANEL_OUT(j+1);
    end
  T_PANEL_OUT(j)=T_PANEL_IN(j)+Q_ABS_PANEL(j)/(M_HTF*CP_HTF);
    T_SURFACE_PANEL(j)=T_PANEL_AVE(j)+Q_ABS_PANEL(j)*(R_CONV_HTF+R_TUBE(j));

end
for j=5:8
   
    if j==1  || j==N_PANELS
        T_PANEL_IN(j)=T_HTF_IN;
    elseif j<=N_PANELS/2 && j~=1 && j~=5
        T_PANEL_IN(j)=T_PANEL_OUT(j-1);
    elseif  j<N_PANELS && j>3*N_PANELS/4
        T_PANEL_IN(j)=T_PANEL_OUT(j+1);
    elseif j==5 
        T_PANEL_IN(j)=T_PANEL_OUT(j+N_PANELS/2);
     elseif j==3*N_PANELS/4
          T_PANEL_IN(j)=T_PANEL_OUT(N_PANELS/4);
    elseif j<3*N_PANELS/4 && j>N_PANELS/2
          T_PANEL_IN(j)=T_PANEL_OUT(j+1);
    end
  T_PANEL_OUT(j)=T_PANEL_IN(j)+Q_ABS_PANEL(j)/(M_HTF*CP_HTF);
    T_SURFACE_PANEL(j)=T_PANEL_AVE(j)+Q_ABS_PANEL(j)*(R_CONV_HTF+R_TUBE(j));

end

for j2=1:4
   j=13-j2;
    if j==1  || j==N_PANELS
        T_PANEL_IN(j)=T_HTF_IN;
    elseif j<=N_PANELS/2 && j~=1 && j~=5
        T_PANEL_IN(j)=T_PANEL_OUT(j-1);
    elseif  j<N_PANELS && j>3*N_PANELS/4
        T_PANEL_IN(j)=T_PANEL_OUT(j+1);
    elseif j==5 
        T_PANEL_IN(j)=T_PANEL_OUT(j+N_PANELS/2);
     elseif j==3*N_PANELS/4
          T_PANEL_IN(j)=T_PANEL_OUT(N_PANELS/4);
    elseif j<3*N_PANELS/4 && j>N_PANELS/2
          T_PANEL_IN(j)=T_PANEL_OUT(j+1);
    end
  T_PANEL_OUT(j)=T_PANEL_IN(j)+Q_ABS_PANEL(j)/(M_HTF*CP_HTF);
    T_SURFACE_PANEL(j)=T_PANEL_AVE(j)+Q_ABS_PANEL(j)*(R_CONV_HTF+R_TUBE(j));

end

T_HTF_AVG_OUT = (T_PANEL_OUT(N_PANELS/2)+T_PANEL_OUT(N_PANELS/2+1))/2;
%% NEW MASS FLOWRATE

M_HTF=Q_ABS_TOT/(N_FLOW_PATH*CP_HTF*(T_HTF_OUT-T_HTF_IN));

M_HTF_TOT=M_HTF*N_FLOW_PATH;

error= abs(T_HTF_AVG_OUT-T_HTF_OUT)/T_HTF_OUT;
end

ETA_THERMAL=Q_ABS_TOT/Q_INCIDENT_TOT;