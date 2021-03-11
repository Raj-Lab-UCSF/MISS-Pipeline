% cd /Users/christophermezias/Documents/MISS_General/MISS-Pipeline;
matdir = '/Users/christophermezias/Documents/MISS_General/MatFiles';
% 
load([matdir filesep 'meanexprmats.mat'],'meanexprmat_ct');
load([matdir filesep 'PresetInputs.mat'],'regvgene','classkey','gene_names','C_indivcells','ct_labvec'); %load for L1 norm E
% load([matdir filesep 'subtype_labels.mat'],'subt_labvec');
% load([matdir filesep 'PresetInputs_rawE.mat'],'regvgene','classkey','gene_names'); %load for raw E
genevct = meanexprmat_ct;
method = 'MRx3';
% testnG = [390:550 750:850];
% testnG = [100:10:1200];
testnG = 50:1000;
% testnG = 810:10:1000;
% testnG = 100;
voxvgene = regvgene;
lambda = 92;
ng_param_list = testnG;
sigmas = 100000;
k = 5;
crossval = 1;
minfun = 'dist2origin';
errtype = 'error';

%tonG_inputmat_norm_inversionalgo
savename = 'superfine_majtypes_reclassify_prefilter_l92';

% preloadinds = MRx3_Selector(genevct,voxvgene,max(ng_param_list),lambda);
% preloadinds = MRx3_Selector_PerCTInit(genevct,voxvgene,max(ng_param_list),lambda);

[fitstruct,~] = nG_ParameterFitter(voxvgene, genevct,...
                                    gene_names, method, C_indivcells,...
                                    ct_labvec, ng_param_list, lambda,...
                                    sigmas, k, crossval, matdir);

% for i = 1:length(ng_param_list)
%     
%     param = ng_param_list(i);
%     classifystruct(i).nGen = param;
%     
%     fprintf('QC metrics calculation, nG parameter value %d/%d\n',i,length(ng_param_list))
%     
%     fprintf('Determining subset, nG parameter value %d/%d\n',i,length(ng_param_list))
%     % Create reduced versions of voxvgene and genevct according to the
%     % method and parameter specified by the user
%     [E_red,C_red,nGen] = GeneSelector(C_indivcells,voxvgene,gene_names,param,lambda,method,preloadinds);
%     
%     [error,accuracy,std_err,~,~,~] = reclassifyCT(C_red,subt_labvec,'recluster',[],5);
% %     [error,accuracy,std_err,~,~,posteriors] = reclassifyCT(C_red,subt_labvec,'GMMPosteriors',0,[]);
% %     [error,accuracy,std_err,~,~,posteriors] = reclassifyCT(C_red,subt_labvec,'GMMPosteriors',1,5);
% %     [error,accuracy,~,~,~,~] = reclassifyCT(C_red,subt_labvec,'SVM',[],5);
%     classifystruct(i).error = error;
%     classifystruct(i).accuracy = accuracy;
%     classifystruct(i).std_err = std_err;
%     
% end

%getting vectors from structs
for i = 1:length(fitstruct)
    negloglike(i) = fitstruct(i).negloglikelihood;
    accuracy(i) = fitstruct(i).gmmstruct.accuracy;
    posteriors(i) = fitstruct(i).gmmstruct.gmmpost;
end

if strcmp(errtype,'loglike')
    error = negloglike;
elseif strcmp(errtype,'error')
    error = 1 - accuracy;
elseif strcmp(errtype,'posts')
    error = 1 - posteriors;
end

%Minmax norms
% norm_error = (negloglike - min(negloglike)) / (max(negloglike) - min(negloglike));
norm_error = (error - min(error)) / (max(error) - min(error));
% post_notright = 1 - posteriors;
% norm_post = (post_notright - min(post_notright)) / (max(post_notright) - min(post_notright));
% normnG = ((ng_param_list - min(ng_param_list)) / (length(gene_names) - min(ng_param_list))).';
normnG = (ng_param_list / (3855)).';
normnG = normnG.';

if strcmp(minfun,'dist2origin')
    dist2origin = sqrt((normnG-0).^2 + (norm_error-0).^2);
    [~,elbowind] = min(dist2origin);
elseif strcmp(minfun,'end+1std')
    elbowind = find(norm_error<norm_error(end)+std(norm_error),'first');
elseif strcmp(minfun,'end+2std')
    elbowind = find(norm_error<norm_error(end)+2*std(norm_error),'first');
elseif strcmp(minfun,'globalmin')
    [~,elbowind] = min(norm_error);
elseif strcmp(minfun,'changept')
    elbowind = findchangepts(norm_error,'Statistic','linear');
elseif strcmp(minfun,'kneept')
    [~,elbowind] = knee_pt(norm_error,normnG);
end

save([matdir filesep savename '.mat'],...
    'fitstruct','lambda','matdir','method','ng_param_list',...
    'elbowind','norm_error','normnG','negloglike',...
    'ng_param_list')

% save([matdir filesep savename '.mat'],...
%     'classkey','lambda','matdir','method','ng_param_list',...
%     'outstruct','preloadinds','LinR_pv','LinR_sst','LinR_vip',...
%     'LinR_micro','tauvec','sumfit_vec','sumfit_fb','sumfit_nctx',...
%     'prodmax','summax','-v7.3')
