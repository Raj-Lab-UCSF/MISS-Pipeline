% cd /Users/christophermezias/Documents/MISS_General/MISS-Pipeline;
matdir = '/Users/christophermezias/Documents/MISS-MatFiles';

% load([matdir filesep 'meanexprmats.mat'],'meanexprmat_ct');
% load([matdir filesep 'PresetInputs.mat'],'regvgene','classkey','gene_names'); %load for L1 norm E
% load([matdir filesep 'PresetInputs_rawE.mat'],'regvgene','classkey','gene_names'); %load for raw E
subtype_extract = 0;
if subtype_extract
    load([matdir filesep 'TasicSub_Inputs.mat'],'voxvgene','classkey','gene_names','ct_labvec','C_indivcells','genevct','ct_group');
    uct = unique(ct_group);
else
    load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','classkey','gene_names','ct_labvec','C_indivcells','genevct');
end
% voxvgene = regvgene;

% genevct = meanexprmat_ct;
tord = 1;
method = 'MRx3';
% testnG = 50:1000;
testnG = [1005:5:1500 1510:10:2000 2050:50:3800];
lambda = 90;
ng_param_list = testnG;

%tonG_inputmat_norm_inversionalgo
savename = 'highernG_resnorms_MRx3Prefilter_l90';

% preloadinds = MRx3_Selector(genevct,voxvgene,max(ng_param_list),lambda);
% preloadinds = MRx3_Selector_PerCTInit(genevct,voxvgene,max(ng_param_list),lambda);
% preloadinds = MRx3_Selector_Prefilter(genevct,voxvgene,max(ng_param_list),lambda);
load([matdir filesep 'MRx3_L90_inds.mat'],'geneinds');
preloadinds = geneinds;

Rval_pv = zeros(length(ng_param_list),3);
Rval_sst = zeros(length(ng_param_list),3);
Rval_vip = zeros(length(ng_param_list),3);
Rval_micro = zeros(length(ng_param_list),3);

for i = 1:length(ng_param_list)
    
    param = ng_param_list(i);
    outstruct(i).nGen = param;
    
    fprintf('QC metrics calculation, nG parameter value %d/%d\n',i,length(ng_param_list))
    
    fprintf('Determining subset, nG parameter value %d/%d\n',i,length(ng_param_list))
    param = ng_param_list(i);
    % Create reduced versions of voxvgene and genevct according to the
    % method and parameter specified by the user
    [E_red,C_red,nGen] = GeneSelector(genevct,voxvgene,gene_names,param,lambda,method,preloadinds);
    
    % Infer cell density per voxel in arbitrary units
    fprintf('Nonnegative matrix inversion, nG parameter value %d/%d\n',i,length(ng_param_list))
%     B = CellDensityInference(E_red,C_red); %default MI, E-L1 ext, C-L1 in
    [B,mnresnorm,fronorm] = CellDensityInference_resnorm(E_red,C_red);
    outstruct(i).Bvals = B;
    outstruct(i).nGen = nGen;
%     outstruct(i).rval = rval;
    outstruct(i).mnresnorm = mnresnorm;
    outstruct(i).fronorm = fronorm;
    if strcmp(method,'MRx3')
        outstruct(i).lambda = lambda;
    end
    
%     Bcorrected = Density_to_Counts_nnnr(B,{1:21,22:24,25},matdir);
%     
%     if subtype_extract
%         Bcorr = zeros(size(Bcorrected,1),length(unique(ct_group)));
%         outstruct(i).Bcorrected = Bcorrected;
%         for j = 1:length(uct)
%             curinds = (ct_group==uct(j));
%             Bcorr(:,j) = sum(Bcorrected(:,curinds),2);
%         end
%         outstruct(i).corrB = Bcorr;
%     else
%         outstruct(i).corrB = Bcorrected;
%     end

    % Sum and average over CCF regions
    [sumB,meanB] = Voxel_To_Region(B,matdir);
    if subtype_extract
        outstruct(i).sumB = sumB;
        outstruct(i).meanB = meanB;
        for j = 1:length(uct)
            curinds = (ct_group==uct(j));
            outstruct(i).Bsums(:,j) = sum(sumB(:,curinds),2);
            outstruct(i).Bmeans(:,j) = sum(meanB(:,curinds),2);
        end
    else
        outstruct(i).Bsums = sumB; % total cells per region
        outstruct(i).Bmeans = meanB; % mean cell count per region
    end

    % Calculate Pearson and Lin R
    [~,PearsonStruct] = CorrelationsCalc_Density(outstruct,i,matdir,tord);
    outstruct(i).Pearson = PearsonStruct;
    Pnames = fieldnames(PearsonStruct);
    for j = 1:length(Pnames)
        curparam_Rval(j,:) = PearsonStruct.(Pnames{j});
    end
    Rval_pv(i,:) = curparam_Rval(1,:);
    Rval_sst(i,:) = curparam_Rval(2,:);
    Rval_vip(i,:) = curparam_Rval(3,:);
    
    fprintf('Done, nG parameter value %d/%d\n',i,length(ng_param_list))
end

save([matdir filesep savename '.mat'],...
    'classkey','lambda','matdir','method','ng_param_list',...
    'outstruct','preloadinds','Rval_pv','Rval_sst','Rval_vip','-v7.3')

