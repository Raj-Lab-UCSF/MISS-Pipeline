function [Zeisel_gene_names,repcells,repcellinds] = Zeisel_geneset(filedir)
%This function generates the unity gene set from the differential
%expression subsetting from Zeisel, et al., 2018 for their correlation
%mapping procedure, and the coronal AGEA ISH atlas.

metacols = h5read([filedir filesep 'l5_all.agg.loom'],'/col_attrs/ClusterName');
load([filedir filesep 'Zeisel_cellIDs.mat'],'cellIncl');
metacols = cellfun(@deblank,metacols,'UniformOutput',false);
[newcellnames,metainds] = sort(metacols);
cellIncl = logical(cellIncl);
inclcols = newcellnames(cellIncl);
exclcols = newcellnames(~cellIncl);
classkey = inclcols;

metacells = h5read([filedir filesep 'l5_all.agg.loom'],'/matrix');
cellcols = h5read([filedir filesep 'l5_all.loom'],'/col_attrs/ClusterName');
cellrows = h5read([filedir filesep 'l5_all.loom'],'/row_attrs/Gene');
newgennames = cellfun(@deblank,cellrows,'UniformOutput',false);
% newcellnames = cellfun(@deblank,cellcols,'UniformOutput',false);
[newgennames,geninds] = sort(newgennames);
metacells = metacells(metainds,geninds);
metacells = metacells(cellIncl,:);
[~,maxgen_percell] = max(metacells,[],2);
unimax = unique(maxgen_percell);
maxgenanmes = newgennames(unimax);

% filedir = '/Users/christophermezias/Documents/MISS_General';
rawdata = readcell([filedir filesep 'mmc4.csv']);
cnames = rawdata(2:4:1058,1);
cnames{39} = 'OBDOP2';
genname_rows = 2:4:1058;
genname_raw = rawdata(genname_rows,3:end);
inclgen = ismember(cnames,classkey);
genname_incl = genname_raw(inclgen,:);
genname_vec = reshape(genname_incl,size(genname_incl,1)*size(genname_incl,2),1);
unigen_spex = unique(genname_vec);
genname_vec = [genname_vec;maxgenanmes];
unigen = unique(genname_vec);

load([filedir filesep 'ISH_gene_names.mat'],'ISH_gene_names');
entrez_inds = ismember(unigen,ISH_gene_names);
entrez_names = unigen(entrez_inds);

criteria = 1;
celltest = ismember(genname_incl,entrez_names);
gensum = sum(celltest,2);
reprows = (gensum>=criteria);
% reprows = unique(reprows);
newcnames = cnames(inclgen);
newcnames = newcnames(reprows);
repcellinds = ismember(classkey,newcnames);
repcells = classkey(repcellinds);

Zeisel_gene_names = entrez_names;

end


