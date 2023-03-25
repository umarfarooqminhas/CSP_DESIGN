%% Clear Workspace
clc;
close all;
clear;
count=1
%% FLUX FROM SOLAR PILOT
for j=7:18
filename=['param_flux-data_Receiver 1_',num2str(j),'.csv']
FLUX_SOLARPILOT = readtable(filename);                              %Average azimuthal flux from Solar Pilot
mmm=FLUX_SOLARPILOT(2:end,2:end);
lens = length(mmm{1,:});
for i=1:lens
    avg_flux(i)=mean(mmm{:,i});
end
avg_flux=avg_flux'
%avg_flux_Tab(j-6,:)=avg_flux
N_PANELS=16;
[FLUX_AVG_PANEL] = Interp_Flux_New(avg_flux, N_PANELS);
%FLUX_AVG_PANEL_Tab(j-6,:)=FLUX_AVG_PANEL

if mean(FLUX_AVG_PANEL)>100
    count=count+1;
    %% Call Receiver Model
    [Q_ABS_TOT,Q_INCIDENT_TOT,ETA_THERMAL,M_HTF] = Receiver_model(FLUX_AVG_PANEL);
    Q_ABS_TOT_Tab(j-6)=Q_ABS_TOT;
    M_HTF_Tab(j-6)=M_HTF;
end

end
%% 
for i=1:24
    if i>=7 && i<=18
        Q_ABS_TOT_Tabb(i)=Q_ABS_TOT_Tab(i-6)/10^6;
    else
        Q_ABS_TOT_Tabb(i)=0;
    end
end
Hours_storage=14;
%% Dispatch Strategy
for i=1:24
    
    Q_NOM=40;
    Q_ABS_TOT_Tab_NOR(i)=Q_ABS_TOT_Tabb(i)/Q_NOM;
    Q_TO_TANK(i)=0;
    hour_stor(1)=3.7;
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
    end
end
yyaxis left;

plot(Q_ABS_TOT_Tab_NOR)

hold on
plot(Q_PB)

hold on
plot(Q_DEF) 

hold on
yyaxis right
plot(hour_stor)

legend('Hour','Q_ABS','Q_PB','Q_DEF')


