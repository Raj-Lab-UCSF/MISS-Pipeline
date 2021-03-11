directory = '/Users/christophermezias/Documents/MISS_General/MatFiles';
classstruct = scRNAseq_Data_Extract(directory);

load([directory filesep 'ISH_data.mat'],'GENname','wherecoronal');
load([directory filesep 'mouse_VISp_2018-06-14_genes-rows.mat'],'filestruct'); 
genesrows = filestruct;

uni_wb_genenames = unique(GENname(wherecoronal));
genebool = ismember(genesrows.gene_symbol,uni_wb_genenames);
ct_entrez = genesrows.gene_entrez_id(genebool);
entrezbool = ismember(classstruct.VISp.EntrezID,ct_entrez);
entrez_names = genesrows.gene_symbol(entrezbool);
members = ismember(filestruct.gene_symbol,entrez_names);

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
                C_indivcells = [C_indivcells classstruct.(regnames{i}).(classnames{j}).(typenames{k}).(subnames{h}).normalizedexpr(members,:)];
                corrs_indivcells = [corrs_indivcells;classstruct.(regnames{i}).(classnames{j}).(typenames{k}).(subnames{h}).corrs];
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

load([directory filesep 'PresetInputs.mat'],'classkey');
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

excl_type = 0;

excltype_inds = (ct_labvec==excl_type);
C_indivcells(:,excltype_inds) = [];
corrs_indivcells(excltype_inds) = [];
ct_labvec(excltype_inds) = [];
subt_labvec(excltype_inds) = [];
typelabs(excltype_inds) = [];
if sum(excl_type) > 0
    classkey(excl_type) = [];
end

corrs_indivcells(isnan(corrs_indivcells)) = 1;

meanexprmat_subt = zeros(size(C_indivcells,1),length(unique(subt_labvec)));
ncells = zeros(length(unique(subt_labvec)),1);
ct_group = zeros(1,length(unique(subt_labvec)));
for i = 1:length(unique(subt_labvec))
    curinds = (subt_labvec==i);
    curcells = C_indivcells(:,curinds);
    curcorrs = corrs_indivcells(curinds);
    meanexprmat_subt(:,i) = curcells * curcorrs / sum(curcorrs);
    ncells(i) = length(find(curinds));
    majct = ct_labvec(curinds);
    ct_group(i) = unique(majct);
end

meanexprmat_ct = zeros(size(C_indivcells,1),length(unique(ct_labvec)));
for i = 1:length(unique(ct_labvec))
    curinds = (ct_group==i);
    cursubs = meanexprmat_subt(:,curinds);
    curncells = ncells(curinds);
    meanexprmat_ct(:,i) = cursubs * curncells / sum(curncells);
end

save([directory filesep 'meanexprmats.mat'],'meanexprmat_subt','meanexprmat_ct','ct_group');
    
    
    


    
    
    