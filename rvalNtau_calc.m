function [Rval_pv,Rval_sst,Rval_vip,tauvec] = rvalNtau_calc(outstruct,matdir)

cellnames = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'}; %laminar glutamatergic cell names
cellinds = 9:15; %laminar glutamatergic cell indices in Bvals or corrB
ranks = [1 2 3 3 4 4 4]; %true layer ordering of laminar glutamatergic cells
tord = 1; %setting correlations to be based on cell density per region, '0' sets it to be totals

Rval_pv = zeros(length(outstruct),3); %preallocating rval and tau vectors across nG values in outstruct
Rval_sst = zeros(length(outstruct),3);
Rval_vip = zeros(length(outstruct),3);
tauvec = zeros(length(outstruct),1);
for i = 1:length(outstruct)
    
    [~,PearsonStruct] = CorrelationsCalc_Density(outstruct,i,matdir,tord); %calculating rvals per nG in outstruct
    Pnames = fieldnames(PearsonStruct);
    for j = 1:length(Pnames)
        curparam_Rval(j,:) = PearsonStruct.(Pnames{j});
    end
    Rval_pv(i,:) = curparam_Rval(1,:);
    Rval_sst(i,:) = curparam_Rval(2,:);
    Rval_vip(i,:) = curparam_Rval(3,:);
    
    taustruct = TauCalc_mod(outstruct,i,cellnames,cellinds,ranks,matdir); %tau value calculation per nG in outstruct
    tauvec(i) = taustruct.tau;
        
end

end