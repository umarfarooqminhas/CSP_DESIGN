%% Economic Analysis
% Input = A_Field, Tower_Height, Reciever_Dia, Receiver_Height, PB_nom,
% Input: D_Tank, H_Tank, Salt_Inventory, AEP (in MWH)
function [LCOE] = LCOE_Cacls(A_field,Tower_Height,Reciever_Dia,Receiver_Height,PB_nom,D_Tank,H_Tank,Salt_Inventory,AEP)
%% Site and Site Improvement Cost
C_site_purchase = 0;                                                             %$/m2
C_site = A_field*C_site_purchase;
C_improvement_ref = 16;                                                          % SAM Figure $/m2
C_site_improvement = C_improvement_ref * A_field;
C_Site_Total = C_site + C_site_improvement;

%% Heliostats Cost

C_Heliostat_Ref = 96;                                                           % $/m2 Reference 2021 cost SunRing NREL Report for 40,000 Heliostats
A_ref_helio = 1078560;                                                          % 26.9 * 40000 Ref nrel.gov/docs/fy22osti/80482.pdf
SF_helio = 1;
Pi_Heliostat=1.08;                                                              %Average CPI
C_Heliostat=C_Heliostat_Ref*A_ref_helio*(A_field/A_ref_helio)^SF_helio*Pi_Heliostat;

%% Tower Cost

C_fixed_tower = 3*10^6;                                                         % [$]Fixed Ref: SAM 
SF_Exponent_Receiver = 0.0113;                                                  % SAM
pi_tower=1.2;                                                                   %Avergae CPI concrete
C_Tower = C_fixed_tower*exp(SF_Exponent_Receiver*(Tower_Height-0.5*Reciever_Dia+0.5*12.2))*pi_tower;  % [$] SAM Model

%% Receiver Cost

C_Rec_ref = 103*10^6;                                                            % [$] Ref SAM
Ref_Area_Rec = 1571;                                                               % m2
SF_rec = 0.7;
pi_rec= 1.35;                                                                    %Average  CPI Stainless Steel Alloy
C_Rec=C_Rec_ref*((pi*(Reciever_Dia)*(Receiver_Height))/Ref_Area_Rec)^SF_rec*pi_rec;

%% PB Cost

C_PB_ref = 122077000;                                                                  % $ Steam Turbine, Steam Cycle and Reated System (Sargent & Lundy Report Case SL-014940- 100 MW 8hr)
PB_ref_size = 100;                                                                  %MWe
SF_PB=0.68;                                                                         %Avg CPI Metal
Pi_PB=1.54;
C_BOF_ref=36413000;                                                                        % $ BOP Foundation Instrumentation and cONTROL Steam Turbine, Steam Cycle and Reated System, Support Structures (Sargent & Lundy Report Case SL-014940- 100 MW 8hr)
SF_BOP=0.68;

C_PB=C_PB_ref*(PB_nom/PB_ref_size)^0.68*Pi_PB+C_BOF_ref*(PB_nom/PB_ref_size)^SF_BOP;

%% TES + Salt Cost                                                                 %Abengoa G018149 2017 Report


CF_HT_REF=10016;
V_HT_REF=pi*(42.4/2)^2*12.2;
CF_CT_REF=4361;
V_HT=pi*(D_Tank/2)^2*H_Tank;

C_TES_BOP=8480;
SF_TES=0.8;
pi_TES=1.35;                                                                    %Avg CPI
C_TES_TOT=(CF_HT_REF*(V_HT/V_HT_REF)^SF_TES+CF_CT_REF*(V_HT/V_HT_REF)^SF_TES+C_TES_BOP*(V_HT/V_HT_REF)^SF_PB)*pi_TES;

C_Molten_Salt=1100;                                                              %$/tonne
C_Salt=Salt_Inventory*C_Molten_Salt;

Contigency = 0.07;

TDC = (C_Site_Total+C_Heliostat+C_Tower+C_Rec+C_PB + C_TES_TOT + C_Salt)*(1+Contigency);
C_Contigency= TDC - (C_Site_Total+C_Heliostat+C_Tower+C_Rec+C_PB + C_TES_TOT + C_Salt);

%% In Direct Cost

EPC = 0.13;                                                                    %SAM
Tax_rate= 0.05;                                                                     %SAM at 80% total Cost
Tax_Base=0.8;
TST= TDC*Tax_Base*Tax_rate;

TIDC= TDC*EPC+TST;

CAPEX = TDC+TIDC;

TIC_SP=TIDC/PB_nom;

%% O&M

OM_fix = 66 * 10^3;                                                           % SAM Fix share of O&M costs [$/(MW*y)]
OM_var = 3.5;                                                                 % SAM Variable share of O&M costs [$/MWh]

figure;
p=pie([C_Site_Total C_Heliostat C_Tower C_Rec C_TES_TOT+C_Salt C_PB  C_Contigency TIDC]);
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Site: ';'Heliostat: ';'Tower: ';'Reciever: ';'TES: ';'PB, AUxiliaries and BOP: ';'Contigency: ';'Indirect: '};  
combinedtxt = strcat(txt,percentValues); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
pText(4).String = combinedtxt(4);
pText(5).String = combinedtxt(5);
pText(6).String = combinedtxt(6);
pText(7).String = combinedtxt(7);
pText(8).String = combinedtxt(8);

%% LCOE Calcs
Discount_Rate=0.1;
Tax_rate=0.25;
Inflation_rate=0.03;
Lifetime = 30;
LCOE=100;
CDCF=10001;
error=10001;
while CDCF>10000 || CDCF<-10000
%while error>10000
      
for i=1:1:Lifetime+1

    C_Depriciation=CAPEX/Lifetime;
    if i==1
       C_OM(i)=0;
       C_TOT(i)=CAPEX+C_OM(i);
       Prod(i)=0;
       Revenue(i)=Prod(i)*LCOE;
       Profit(i)=Revenue(i)-C_TOT(i);
       %Tax_y(j-1)=Profit(j-1)*Tax_rate;
       Tax_y(i)=0;
       DCF(i)=(Revenue(i)-C_TOT(i)-Tax_y(i))*(1+Discount_Rate)^(-(i-1));
    else
       C_OM(i)=(PB_nom*OM_fix+OM_var*AEP)*(1+Inflation_rate)^(i-1);
       C_TOT(i)=C_OM(i);
       Prod(i)=AEP;
       Revenue(i)=Prod(i)*LCOE;
       Profit(i)=Revenue(i)-C_TOT(i)-C_Depriciation;
       Tax_y(i)=Profit(i)*Tax_rate;
       DCF(i)=(Revenue(i)-C_TOT(i)-Tax_y(i))*(1+Discount_Rate)^(-(i-1));
    end
end  
 CDCF=sum(DCF);
    if CDCF>1
        LCOE=LCOE-0.001;
    elseif CDCF <1
        LCOE=LCOE+0.001;
    end
end