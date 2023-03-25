%% Fluid properties
function [CP_HTF] = CP_Salt_Props(T_HTF_REF)

load 'CP MOLTEN SALT.txt';


for j=1:length(CP_MOLTEN_SALT(:,1))
    if CP_MOLTEN_SALT(j,1)>T_HTF_REF
        CP_HTF=CP_MOLTEN_SALT(j,2);
        break
    else 
      CP_HTF=CP_MOLTEN_SALT(length(CP_MOLTEN_SALT(:,1)),2);  
    end
   
    j=j+1;
end

end