%% Clear Workspace
clc;
close all;
clear;
count=1;
i=1;
%% FLUX FROM SOLAR PILOT
Tamb_v=AmbientTemp();                                                       %NEOM Data
for j=1:6672
    omit=[1417:1:1488,	2953:1:2976 ,4441	:1:4464	, 6673	:1:6696, 8161 :1:8184];


    if any(omit(:) == j)
    else

    
    avg_flux=0;
    filename=['param_flux-data_Receiver 1_',num2str(j),'.csv'];
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
        [Q_ABS_TOT,Q_INCIDENT_TOT,ETA_THERMAL,M_HTF] = Receiver_model(FLUX_AVG_PANEL,23,27.6,250,Tamb_v(i));
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
for j=1:2232
    omit=[1417 :1:1440];


    if any(omit(:) == j)
    else

    
    avg_flux=0;
    filename=['m_flux-data_Receiver 1_',num2str(j),'.csv'];
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
        [Q_ABS_TOT,Q_INCIDENT_TOT,ETA_THERMAL,M_HTF] = Receiver_model(FLUX_AVG_PANEL,23,27.6,250,Tamb_v(i));
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
%% %% Dispatch Strategy

for h=14
 Hours_storage=h;

i=1;
for i=1:8760
    Q_NOM=345;
    Q_ABS_TOT_Tab_NOR(i)=Q_ABS_TOT_Tab(i)/Q_NOM;
    Q_TO_TANK(i)=0;
    hour_stor(h,1)=0;
    if Q_ABS_TOT_Tab_NOR(i)>1
        Q_PB(h,i)=1;
        Q_TO_TANK=Q_ABS_TOT_Tab_NOR(i)-Q_PB(h,i);
        Q_DEF_Check= hour_stor(h,i-1)-Hours_storage;
                if Q_DEF_Check>1
                    Q_DEF(h,i) = Q_DEF_Check;
                    hour_stor(h,i)=Hours_storage;
                else
                     hour_stor(h,i)=Q_TO_TANK+hour_stor(h,i-1);
                     if   hour_stor(h,i)> Hours_storage
                         Q_DEF(h,i)=hour_stor(h,i)-Hours_storage;
                         hour_stor(h,i)=Hours_storage;
                     else
                       Q_DEF(h,i)=0;

                     end
                end
          
    elseif i>1
        if Q_ABS_TOT_Tab_NOR(i)<1 && hour_stor(h,i-1)>0 
         Q_DEF(h,i)=0;
            if  Q_ABS_TOT_Tab_NOR(i)+hour_stor(h,i-1)>1
             Q_PB(h,i)=1;
             hour_stor(h,i)=hour_stor(h,i-1)-1;
            else
                Q_PB(h,i)=Q_ABS_TOT_Tab_NOR(i)+hour_stor(h,i-1);
                hour_stor(h,i)=0;
            end
        else
          Q_PB(h,i)=0;
          hour_stor(h,i)=0;
        end
    elseif i==1 && hour_stor(1)>0
        Q_PB(h,i)=1;
        hour_stor(h,i)=hour_stor(h,i)-1;
    else
         Q_PB(h,i)=0;
    end
end


%% Ploting
%{

 yyaxis left;

plot(Q_ABS_TOT_Tab_NOR(,'r--')

hold on
plot(Q_PB)

hold on
plot(Q_DEF,'k+--') 

hold on
yyaxis right
plot(hour_stor,'b*-')

legend('Q ABS','Q PB','Q_DEF','Hour') 

%}
%% 
for t=1:12

Tower_Height=195+t*5;
[CP_HTF,MU_HTF,K_HTF,RHO_HTF] = Molten_Salt_Props((565+260)/2);
eta_pump=0.80;
eta_EM=0.98;
ETA_AVG_PB=0.4;
DP_SF = RHO_HTF*9.81*Tower_Height;
for j=1:8760
    if Q_PB(h,j)>0
        M_HTF_IN=Q_PB(h,j)*Q_NOM*10^6/(CP_HTF*(560-260))
        M_HTF_IN_PB_V(j)=M_HTF_IN;
        %[W_NET,ETA_PLANT] = Power_Block_OFF(M_HTF_IN,Tamb_v(i));                        %T_AMB(j)
        W_NET=ETA_AVG_PB*Q_PB(h,j)*Q_NOM;
        W_Gross_Tab(j)=W_NET;
       % ETA_PLANT_Tab(j)=ETA_PLANT;
        if Q_ABS_TOT_Tab(j)>0;
            P_Aux_pump(j) = (M_HTF_IN*DP_SF)/RHO_HTF/eta_pump/eta_EM+( 345*15774)/RHO_HTF/eta_pump/eta_EM;
        else
            P_Aux_pump(j) = ( 345*15774)/RHO_HTF/eta_pump/eta_EM;      %% F=0.015, V=3.41, D=16", Mu=0.001835     
        end
        W_NET_Tab(j)=W_Gross_Tab(j)-P_Aux_pump(j)/10^6;
    else
         W_NET_Tab(j)=0;
         ETA_PLANT_Tab(j)=0;  
  
    end

   
end
AEP(h)= sum(W_NET_Tab);
%%


A_field=1843666;

Reciever_Dia=23;
Receiver_Height=27.6;
PB_nom=150;                                                                 %MW

%% Hours Storage Sensitivity Check

    Salt_Inventory=(h+1)*3600*741/1000;                                     %tonne
    V_Storage=Salt_Inventory*1000/1785;                                     %Density 1785 kg/m3
    H_Tank=12.2;
    D_tank_guess=12.2;
    err_tank_D=1;
    while err_tank_D>10^-4
        D_tank=sqrt(4*V_Storage/(pi()*H_Tank));
        err_tank_D=abs((D_tank-D_tank_guess)/D_tank);
        D_tank_guess=D_tank;
    end
     D_tank= D_tank_guess;                           
    [LCOE] = LCOE_Cacls(A_field,Tower_Height,Reciever_Dia,Receiver_Height,PB_nom,D_tank,H_Tank,Salt_Inventory,AEP(h));
    LCOE_h(h)=LCOE;
end
end

%{

ETA_OVERALL=AEP./5167690;                                                      %5167690 (DNI*A_Field/10^6) sum along year

%% Defocus at eaach our

for h=1:18
    Q_DEF_h(h)=sum(Q_DEF(h,:));
end
%%
xaxis=[10 11 12 13 14 15 16 17 18];
yyaxis left;
plot(xaxis,Q_DEF_h(10:18)./(4336),'-o',LineWidth=1.5);
ylabel('Defocus [%]')
hold on
yyaxis right;
plot(xaxis,LCOE_h(10:18),'-o',LineWidth=1.5)
ylabel('LCOE [$/MWh]')
xlabel('Stoarge Hours [h]');
ylim([100 150]);
legend('Defocus %','LCOE [$/MWh]')

%% LCOE 2 Graph with SAM
figure;
xaxis=[1:1:16];
%load 'LCOE_h_SM3.mat';
plot(xaxis,LCOE_h(1:16),'-o',LineWidth=1.5)
ylabel('LCOE [$/MWh]')
xlabel('Stoarge Hours [h]');
%ylim([100 150]);

hold on;
SAM_LCOE=readtable('LCOE storage Hours.csv');
SAM_LCOE_v=SAM_LCOE{:,2}*10';

plot(xaxis,SAM_LCOE_v(1:16),'-o',LineWidth=1.5)
ylabel('LCOE [$/MWh]')
xlabel('Stoarge Hours [h]');
ylim([80 260]);
legend('LCOE [$/MWh] OUR MODEL', 'LCOE [$/MWh] SAM MODEL')
%}

