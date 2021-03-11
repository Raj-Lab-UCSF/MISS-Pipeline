function [meanexprmat_ct, meanexprmat_subt, classkey, classkey_subt, ...
    C_indivcells, ct_labvec, subt_labvec, ct_group, entrez_names] = ...
    Cell_Type_Data_Extract_Tasic(classstruct,excl_names,directory)

if nargin < 2
    excl_names = 'VLMC';
    if nargin < 3
        directory = [cd filesep 'MatFiles'];
    end
end

load([directory filesep 'ISH_data.mat'],'GENname','wherecoronal');
load([directory filesep 'mouse_VISp_2018-06-14_genes-rows.mat'],'filestruct'); 
genesrows = filestruct;

uni_wb_genenames = unique(GENname(wherecoronal));
genebool = ismember(genesrows.gene_symbol,uni_wb_genenames);
ct_entrez = genesrows.gene_entrez_id(genebool);
entrezbool = ismember(classstruct.VISp.EntrezID,ct_entrez);
entrez_names = genesrows.gene_symbol(entrezbool);

C_indivcells = [];
corrs_indivcells = [];
regnames = fieldnames(classstruct);
sublabs = {};
typelabs = {};

for i = 1:length(regnames)
    classnames = fieldnames(classstruct.(regnames{i}));
    exclinds = ismember(classnames,{'EntrezID','Outlier','Low_Quality'});
    classnames(exclinds) = [];
    for j = 1:length(classnames)
        typenames = fieldnames(classstruct.(regnames{i}).(classnames{j}));
        typenames([1 end-1:end]) = [];
        for k = 1:length(typenames)
            curtype = typenames{k};
            subnames = fieldnames(classstruct.(regnames{i}).(classnames{j}).(typenames{k}));
            subnames([1 end-1:end]) = [];
            for h = 1:length(subnames)
                C_indivcells = [C_indivcells classstruct.(regnames{i}).(classnames{j}).(typenames{k}).(subnames{h}).normalizedexpr(entrezbool,:)];
                corrs_indivcells = [corrs_indivcells; classstruct.(regnames{i}).(classnames{j}).(typenames{k}).(subnames{h}).corrs];
                cursublabs = cell(1,classstruct.(regnames{i}).(classnames{j}).(typenames{k}).(subnames{h}).cellcount);
                cursub = subnames{h};
                cursublabs(:) = {cursub};
                sublabs = [sublabs cursublabs];
                curtypelabs = cell(size(cursublabs));
                curtypelabs(:) = {curtype};
                typelabs = [typelabs curtypelabs];
            end
        end
    end
end

for i = 1:length(typelabs)
    curtype = typelabs{i};
    if strcmp(curtype,'Olig1')
        typelabs{i} = 'Oligo';
    elseif strcmp(curtype,'OPC')
        typelabs{i} = 'Oligo';
    elseif strcmp(curtype,'Macrophage')
        typelabs{i} = 'Macro';
    elseif strcmp(curtype,'Micro')
        typelabs{i} = 'Macro';
    elseif strcmp(curtype,'Peri')
        typelabs{i} = 'Endo';
    elseif strcmp(curtype,'SMC')
        typelabs{i} = 'Endo';
    end
end

classkey = unique(typelabs(~ismember(typelabs,excl_names)),'stable');
nonneurontypes = {'Astro','Macro','Oligo','Endo'};
classkey = setdiff(classkey,nonneurontypes,'stable');
classkey = [classkey, nonneurontypes];
ct_labvec = zeros(size(typelabs));
subt_labvec = zeros(size(sublabs));
dex = 1;
for i = 1:length(classkey)
    membs = ismember(typelabs,classkey{i});
    ct_labvec(membs) = i;
    cursubkeys = unique(sublabs(membs));
    for j = 1:length(cursubkeys)
        submembs = ismember(sublabs,cursubkeys{j});
        subt_labvec(submembs) = dex;
        dex = dex + 1;
    end 
end

excltype_inds = (ct_labvec==0);
C_indivcells(:,excltype_inds) = [];
corrs_indivcells(excltype_inds) = [];
ct_labvec(excltype_inds) = [];
subt_labvec(excltype_inds) = [];
typelabs(excltype_inds) = [];
sublabs(excltype_inds) = [];

meanexprmat_subt = zeros(size(C_indivcells,1),length(unique(subt_labvec)));
ncells = zeros(length(unique(subt_labvec)),1);
ct_group = zeros(1,length(unique(subt_labvec)));
classkey_subt = cell(1,length(ct_group));
for i = 1:length(unique(subt_labvec))
    curinds = (subt_labvec==i);
    curcells = C_indivcells(:,curinds);
    curcorrs = corrs_indivcells(curinds);
    if ~all(isnan(curcorrs))
        curcorrs(isnan(curcorrs)) = nanmean(curcorrs);
    else
        curcorrs(isnan(curcorrs)) = 1;
    end
    meanexprmat_subt(:,i) = curcells * curcorrs / sum(curcorrs);
    ncells(i) = sum(curinds);
    ct_group(i) = unique(ct_labvec(curinds));
    classkey_subt{i} = unique(sublabs(curinds));
end

meanexprmat_ct = zeros(size(C_indivcells,1),length(unique(ct_labvec)));
for i = 1:length(unique(ct_labvec))
    curinds = (ct_group==i);
    cursubs = meanexprmat_subt(:,curinds);
    curncells = ncells(curinds);
    meanexprmat_ct(:,i) = cursubs * curncells / sum(curncells);
end
end
    


    
    
    