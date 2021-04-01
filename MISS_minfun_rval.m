matdir = '/Users/christophermezias/Documents/MISS-MatFiles';
load([matdir filesep 'resnorm_MRx3Prefilter_l90.mat'],'outstruct','Rval_pv','Rval_sst','Rval_vip');
load([matdir filesep 'MRx3_L90_inds.mat'],'geneinds');
preloadinds = geneinds;
load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','classkey','gene_names','ct_labvec','C_indivcells','genevct');
method = 'MRx3';
lambda = 90;

rvals = zeros(length(outstruct),1);
ng_param_list = zeros(1,length(outstruct));
for i = 1:length(outstruct)
    param = outstruct(i).nGen;
    [E_red,C_red,nGen] = GeneSelector(genevct,voxvgene,gene_names,param,lambda,method,preloadinds);
    Bvals = outstruct(i).Bvals;
    Bvals = Bvals.';
    CB_red = C_red * Bvals;
    rvals(i) = corr(E_red(:),CB_red(:));
    ng_param_list(i) = param;
end

save([matdir filesep 'rval_MRx3Prefilter_l90.mat'],'rvals','Rval_pv','Rval_sst','Rval_vip','ng_param_list','lambda');    
    