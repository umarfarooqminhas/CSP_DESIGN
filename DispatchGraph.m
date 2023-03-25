%% Clear Workspace
clc;
close all;
clear;
count=1;
i=1;
%% FLUX FROM SOLAR PILOT
Tamb_v=AmbientTemp();                                                       %NEOM Data
for j=4596:4596
    omit=[1417:1:1488,	2953:1:2976 ,4441	:1:4464	, 6673	:1:6696, 8161 :1:8184];


    if any(omit(:) == j)
    else

    
    avg_flux=0;
    filename=['param_flux-data_Receiver 1_',num2str(j),'.csv'];
    %filename=['m_flux-data_Receiver 1_',num2str(j),'.csv'];
    try
        FLUX_SOLARPILOT = readtable(filename);                              %Average azimuthal flux from Solar Pilot
        mmm=FLUX_SOLARPILOT(2:end,2:end);
    lens = length(mmm{1,:});
        for l=1:lens
            avg_flux(l)=mean(mmm{:,l});
        end
    avg_flux=avg_flux';
    %avg_flux_Tab(j-6,:)=avg_flux
    N_PANELS=16;
    [FLUX_AVG_PANEL] = Interp_Flux_New(avg_flux, N_PANELS);
    %FLUX_AVG_PANEL_Tab(j-6,:)=FLUX_AVG_PANEL
    catch
         Q_ABS_TOT_Tab(i)=0;
         M_HTF_Tab(i)=0;
         FLUX_AVG_PANEL=0;
    end

    if mean(FLUX_AVG_PANEL)>150
        count=count+1;
    %Call Receiver Model
        [Q_ABS_TOT,Q_INCIDENT_TOT,ETA_THERMAL,M_HTF] = Receiver_model(FLUX_AVG_PANEL,23,27.6,250,25);
        Q_ABS_TOT_Tab(i)=Q_ABS_TOT/10^6;                                    %Megawatt
        M_HTF_Tab(i)=M_HTF;
        ETA_THERMAL_Tab(i)=ETA_THERMAL;
    else
        Q_ABS_TOT_Tab(i)=0;
        M_HTF_Tab(i)=0; 
    end
    i=i+1
    end
end

%% Dispatching
 Hours_storage=14;
for i=1:24
    Q_NOM=345;
    Q_ABS_TOT_Tab_NOR(i)=Q_ABS_TOT_Tab(i)/Q_NOM;
    Q_TO_TANK(i)=0;
    hour_stor(1)=1;
    if Q_ABS_TOT_Tab_NOR(i)>1
        Q_PB(i)=1;
        Q_TO_TANK=Q_ABS_TOT_Tab_NOR(i)-Q_PB(i);
        Q_DEF_Check= hour_stor(i-1)-Hours_storage;
                if Q_DEF_Check>1
                    Q_DEF(i) = Q_DEF_Check;
                    hour_stor(i)=Hours_storage;
                else
                     hour_stor(i)=Q_TO_TANK+hour_stor(i-1);
                     if   hour_stor(i)> Hours_storage
                         Q_DEF(i)=hour_stor(i)-Hours_storage;
                         hour_stor(i)=Hours_storage;
                     else
                       Q_DEF(i)=0;

                     end
                end
          
    elseif i>1
        if Q_ABS_TOT_Tab_NOR(i)<1 && hour_stor(i-1)>0 
         Q_DEF(i)=0;
            if  Q_ABS_TOT_Tab_NOR(i)+hour_stor(i-1)>1
             Q_PB(i)=1;
             hour_stor(i)=hour_stor(i-1)-1;
            else
                Q_PB(i)=Q_ABS_TOT_Tab_NOR(i)+hour_stor(i-1);
                hour_stor(i)=0;
            end
        else
          Q_PB(i)=0;
          hour_stor(i)=0;
        end
    elseif i==1 && hour_stor(1)>0
        Q_PB(i)=1;
        hour_stor(i)=hour_stor(i)-1;
    else
         Q_PB(i)=0;
    end
end

%% Ploting

xaxis=0:1:23;
yyaxis left;


p=plot(xaxis,Q_ABS_TOT_Tab_NOR(),'o-',LineWidth=1.5);
c = p.Color;
p.Color = "#EDB120";

hold on
plot(xaxis,Q_PB,'k+-',LineWidth=1.5);

hold on
plot(xaxis,Q_DEF,'k+-',LineWidth=1.5) ;
ylabel('Q ABSORDED Q PB QDEF x 345 MW',LineWidth=1.5);
hold on
yyaxis right

plot(xaxis,hour_stor,'r*-',LineWidth=1.5);
ylabel('State of Charge [h]',LineWidth=1.5);
xlim([0 23]);

legend('Q ABS [MW]','Q PB [MW]','Q DEFS[MW]','TES SOC [h]',Location='northwest') 

