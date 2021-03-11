% cd /Users/christophermezias/Documents/MISS_General/MISS-Pipeline;
matdir = '/Users/christophermezias/Documents/MISS_General/MatFiles';
% 
load([matdir filesep 'meanexprmats.mat'],'meanexprmat_ct');
load([matdir filesep 'PresetInputs.mat'],'regvgene','classkey','gene_names','C_indivcells','ct_labvec'); %load for L1 norm E
genevct = meanexprmat_ct;
method = 'MRx3';
% testnG = [368 388];
testnG = 408;
voxvgene = regvgene;
sigmas = 100000;
lambdas = 92;
lnames = {'l92'};
% lambdas = [98 99 95];
% lnames = {'l98','l99','l95'};
% lambdas = [250 250];
% lnames = {'posterior','error'};
k = 5;
crossval = 1;
minfun = 'dist2origin';
% errnames = {'error','posterior'};

%tonG_inputmat_norm_inversionalgo
savename = 'majtypes_elbows_dist2origin_prefilter_l92';

% preloadinds = MRx3_Selector(genevct,voxvgene,max(ng_param_list),lambda);
% preloadinds = MRx3_Selector_PerCTInit(genevct,voxvgene,max(ng_param_list),lambda);

for i = 1:length(testnG)
    
    lambda = lambdas(i);
    ng_param_list = testnG(i);

    [fitstruct_,outstruct_] = nG_ParameterFitter(voxvgene, genevct,...
                                    gene_names, method, C_indivcells,...
                                    ct_labvec, ng_param_list, lambda,...
                                    sigmas, k, crossval, matdir);
                                
    outstruct.(lnames{i}) = outstruct_;
    fitstruct.(lnames{i}) = fitstruct_;
                                
end

save([matdir filesep savename '.mat'],...
    'fitstruct','lambda','matdir','method','ng_param_list','outstruct');

% Figure_4ab_taulayerslice(outstruct.error,1,[25 32 48],0,matdir)
% Figure_4ab_taulayerslice(outstruct.posterior,1,[25 32 48],0,matdir)
