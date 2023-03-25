%% Fluid properties
function [CP_HTF,MU_HTF,K_HTF,RHO_HTF] = Molten_Salt_Props(T_HTF_REF)

load 'CP MOLTEN SALT.txt';
load 'Density Kg-m-3 Molten Salt.txt';
load 'K W-m-k Molten Salt.txt';
load 'mu Pa-s Molten Salt.txt';


for j=1:length(CP_MOLTEN_SALT(:,1))
    if CP_MOLTEN_SALT(j,1)>T_HTF_REF
        CP_HTF=CP_MOLTEN_SALT(j,2);
        break
    end
    j=j+1;
end

for j=1:length(Density_Kg_m_3_Molten_Salt(:,1))
    if Density_Kg_m_3_Molten_Salt(j,1)>T_HTF_REF
        RHO_HTF=Density_Kg_m_3_Molten_Salt(j,2);
        break
    end
    j=j+1;
end

for j=1:length(K_W_m_k_Molten_Salt(:,1))
    if K_W_m_k_Molten_Salt(j,1)>T_HTF_REF
        K_HTF=K_W_m_k_Molten_Salt(j,2);
        break
    end
    j=j+1;
end

for j=1:length(mu_Pa_s_Molten_Salt(:,1))
    if mu_Pa_s_Molten_Salt(j,1)>T_HTF_REF
        MU_HTF=mu_Pa_s_Molten_Salt(j,2);
        break
    end
    j=j+1;
end