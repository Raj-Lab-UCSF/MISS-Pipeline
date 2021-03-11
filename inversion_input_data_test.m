% cd /Users/christophermezias/Documents/MISS_General/MISS-Pipeline;
matdir = '/Users/christophermezias/Documents/MISS_General/MatFiles';

% load([matdir filesep 'meanexprmats.mat'],'meanexprmat_ct');
% load([matdir filesep 'PresetInputs.mat'],'regvgene','classkey','gene_names'); %load for L1 norm E
% load([matdir filesep 'PresetInputs_rawE.mat'],'regvgene','classkey','gene_names'); %load for raw E
subtype_extract = 0;
if subtype_extract
    load([matdir filesep 'meanexprmats.mat'],'meanexprmat_subt','ct_group');
    load([matdir filesep 'PresetInputs.mat'],'regvgene','classkey','gene_names','ct_labvec','C_indivcells');
    genevct = meanexprmat_subt;
    uct = unique(ct_group);
else
    load([matdir filesep 'PresetInputs.mat'],'regvgene','classkey','gene_names','ct_labvec','C_indivcells');
    load([matdir filesep 'meanexprmats.mat'],'meanexprmat_ct');
    genevct = meanexprmat_ct;
end
% voxvgene = regvgene;

% genevct = meanexprmat_ct;
tord = 1;
method = 'MRx3';
% testnG = [390:550 750:850];
% testnG = [399 408 414];
testnG = 399;
% testnG = 810:10:1000;
% testnG = [350 375 400 425 450];
% voxvgene = regvgene; %This line for E as is from file
normtype = 'L2';
excl_thalamus = 0;
voxvgene = regvgene;
% gnrinds = {1:21,22:25};
% if excl_thalamus
%     exclregs = 17:21;
%     exclinds = exclregs;
%     classkey(exclinds) = [];
%     genevct(:,exclinds) = [];
%     voxvgene = E2C_Erescale_nothal(regvgene,genevct,normtype,1,matdir); %This line for E2C rescale
% else
%     voxvgene = E2C_Erescale(regvgene,genevct,normtype,1,matdir); %This line for E2C rescale
% end 
lambda = 98;
ng_param_list = testnG;

%tonG_inputmat_norm_inversionalgo
savename = 'reclassify_elbowinds_MRx3Prefilter_l98';

% preloadinds = MRx3_Selector(genevct,voxvgene,max(ng_param_list),lambda);
% preloadinds = MRx3_Selector_PerCTInit(genevct,voxvgene,max(ng_param_list),lambda);
preloadinds = MRx3_Selector_Prefilter(genevct,voxvgene,max(ng_param_list),lambda);

tauvec = zeros(length(ng_param_list),1);
LinR_pv = zeros(length(ng_param_list),3);
LinR_sst = zeros(length(ng_param_list),3);
LinR_vip = zeros(length(ng_param_list),3);
LinR_micro = zeros(length(ng_param_list),3);
Rval_pv = zeros(length(ng_param_list),3);
Rval_sst = zeros(length(ng_param_list),3);
Rval_vip = zeros(length(ng_param_list),3);
Rval_micro = zeros(length(ng_param_list),3);
sumfit_vec = zeros(length(ng_param_list),1);

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
    [B,resnorm,fronorm] = CellDensityInference_resnorm(E_red,C_red); %default MI, E-L1 ext, C-L1 in
%     B = CellDensityInference_L2norm(E_red,C_red); %lsqnonneg, no norms
%     B = CellDensityInference_nonorm(E_red,C_red); %lsqnonneg, no norms
%       B = CellDensityInference_lsqlin(E_red,C_red); %lsqlin, L1s
%       B = CellDensityInference_lsqlin_nonorm(E_red,C_red); %lsqlin, none
%     B = CellDensityInference_lsqlin_L2norm(E_red,C_red); %lsqlin, L2s
%     [B,intcpts] = CellDensityInference_intcpt_nonorm(E_red,C_red); %lsqnonneg, w/intcp, none
%     [B,intcpt] = CellDensityInference_lsqlin_intcpt(E_red,C_red);
%     B = CellDensityInference_lsqlin_nnnr(E_red,C_red);
%     B = CellDensity_Corr(E_red,C_red);
%     outstruct(i).intcpt = intcpt;
%     outstruct(i).resnorm = resnorm;
%     outstruct(i).fronorm = fronorm;
    outstruct(i).Bvals = B;
    outstruct(i).nGen = nGen;
    if strcmp(method,'MRx3')
        outstruct(i).lambda = lambda;
    end
    
    Bcorrected = Density_to_Counts_nnnr(B,{1:21,22:24,25},matdir);
%     outstruct(i).Bvals = B;
%     Bcorrected = Density_to_Counts(B,matdir);
%     Bcorrected = Density_to_Counts_nnnr(B,{find(ismember(ct_group,1:21)),find(ismember(ct_group,22:24)),find(ismember(ct_group,25))},matdir);
%     [Bcorrected,Binterneurons,Bfactor] = Density_to_Counts_fit2data(B,0,'neo',matdir);
%     [Bcorrected,Bfactor] = Density_to_Counts_GNR(B,gnrinds,matdir);
 
%     outstruct(i).corrB = Binterneurons;
%     outstruct(i).corrB = Bcorrected;
    
    if subtype_extract
        Bcorr = zeros(size(Bcorrected,1),length(unique(ct_group)));
        outstruct(i).Bcorrected = Bcorrected;
        for j = 1:length(uct)
            curinds = (ct_group==uct(j));
            Bcorr(:,j) = sum(Bcorrected(:,curinds),2);
        end
        outstruct(i).corrB = Bcorr;
    else
        outstruct(i).corrB = Bcorrected;
    end

    % Sum and average over CCF regions
%     [sumB,meanB] = Voxel_To_Region(Binterneurons,matdir);
    [sumB,meanB] = Voxel_To_Region(Bcorrected,matdir);
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
    
%     outstruct(i).Bsums = sumB; % total cells per region
%     outstruct(i).Bmeans = meanB; % mean cell count per region

    % Calculate Pearson and Lin R
    [LinRstruct,PearsonStruct] = CorrelationsCalc_Density(outstruct,i,matdir,tord);
    outstruct(i).LinR = LinRstruct;
    outstruct(i).Pearson = PearsonStruct;
    LinRnames = fieldnames(LinRstruct);
    Pnames = fieldnames(PearsonStruct);
    for j = 1:length(LinRnames)
        curparam_LinR(j,:) = LinRstruct.(LinRnames{j});
        curparam_Rval(j,:) = PearsonStruct.(Pnames{j});
    end
    LinR_pv(i,:) = curparam_LinR(1,:);
    LinR_sst(i,:) = curparam_LinR(2,:);
    LinR_vip(i,:) = curparam_LinR(3,:);
%     LinR_all(i,:) = curparam_LinR(4,:);
    LinR_micro(i,:) = curparam_LinR(4,:);
%     LinR_neuron(i,:) = curparam_LinR(6,:);
    Rval_pv(i,:) = curparam_Rval(1,:);
    Rval_sst(i,:) = curparam_Rval(2,:);
    Rval_vip(i,:) = curparam_Rval(3,:);
%     Rval_all(i,:) = curparam_Rval(4,:);
    Rval_micro(i,:) = curparam_Rval(4,:);
%     Rval_neuron(i,:) = curparam_Rval(6,:);

    % Calculate adjusted Kendall's tau for layer-type glutamatergic neurons
    fprintf('Calculating tau, nG parameter value %d/%d\n',i,length(ng_param_list))
    
    majinds = 9:15;
    cnames = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
    
    ranks = [1 2 3 3 4 4 4];
    cell_inds = majinds;
    cell_names = cnames;
    
    taustruct = TauCalc_mod(outstruct,i,cell_names,cell_inds,ranks,matdir);
    tauvec(i) = taustruct.tau;
    outstruct(i).tau = taustruct.tau;

    % Calculate SumFit criterion
    fprintf('Calculating sum fit, nG parameter value %d/%d\n',i,length(ng_param_list))
    sumfit = taustruct.tau + LinR_pv(i,3) + LinR_sst(i,3) + LinR_vip(i,3) + LinR_micro(i,3);
    sumfit_vec(i) = sumfit;
    outstruct(i).sumfit = sumfit;
    
    fprintf('Done, nG parameter value %d/%d\n',i,length(ng_param_list))
end

sumfit_fb = LinR_pv(:,2) + LinR_sst(:,2) + LinR_vip(:,2) + LinR_micro(:,2) + tauvec;
sumfit_nctx = LinR_pv(:,1) + LinR_sst(:,1) + LinR_vip(:,1) + LinR_micro(:,1) + tauvec;

nctxprod = LinR_pv(:,1) .* LinR_sst(:,1) .* LinR_vip(:,1) .* tauvec;
nctxsum = LinR_pv(:,1) + LinR_sst(:,1) + LinR_vip(:,1) + tauvec;
[~,prodmax] = max(nctxprod);
[~,summax] = max(nctxsum);

save([matdir filesep savename '.mat'],...
    'classkey','lambda','matdir','method','ng_param_list',...
    'outstruct','preloadinds','LinR_pv','LinR_sst','LinR_vip',...
    'LinR_micro','tauvec','sumfit_vec','sumfit_fb','sumfit_nctx',...
    'prodmax','summax','-v7.3')

% save([matdir filesep savename '.mat'],...
%     'classkey','lambda','matdir','method','ng_param_list',...
%     'outstruct','preloadinds','LinR_pv','LinR_sst','LinR_vip',...
%     'LinR_micro','tauvec','sumfit_vec','sumfit_fb','sumfit_nctx',...
%     'prodmax','summax','-v7.3')
