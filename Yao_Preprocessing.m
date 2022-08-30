clear; clc;

%% Metadata processing for grouping clusters into subclasses
% Define directories
glbl = '/Users/justintorok/Documents/MATLAB/MISS';
lcl_sc = 'RawYao';
lcl_ish = 'MISS-MatFiles';

% Average trimmed means gene expression data per cluster over subclasses
file_yao_metadata = [glbl filesep lcl_sc filesep 'yao_metadata.csv'];
yao_metadata = readtable(file_yao_metadata);
subclasses = yao_metadata.subclass_label;
classkey = unique(subclasses);
clusters = yao_metadata.cluster_label;
classkey_clust = unique(clusters);
subclass_id_vec = zeros(1,length(classkey_clust));
numcells_cluster = subclass_id_vec;

file_yao_trimmed_means= [glbl filesep lcl_sc filesep 'yao_trimmed_means.csv'];
yao_trimmed_means = readtable(file_yao_trimmed_means,'VariableNamingRule','preserve');
gene_names_yao = yao_trimmed_means.feature;
[gene_names_yao,gene_names_sortinds] = sort(gene_names_yao);

for i = 1:length(classkey_clust)
    cluster_i = classkey_clust{i};
    cluster_i_inds = ismember(clusters,cluster_i);
    numcells_cluster(i) = sum(cluster_i_inds);
    subclass_i = unique(subclasses(cluster_i_inds));
    subclass_id_vec(i) = find(ismember(classkey,subclass_i));
end

genevct_allgenes = zeros(length(gene_names_yao),length(classkey));
for i = 1:length(classkey)
    clusters_i = classkey_clust(subclass_id_vec == i);
    numcells_i = numcells_cluster(subclass_id_vec == i);
    trimmed_sum_i = zeros(length(gene_names_yao),1);
    for j = 1:length(clusters_i)
        trimmed_mean_cluster_j = yao_trimmed_means.(clusters_i{j});
        trimmed_sum_i = trimmed_sum_i + (numcells_i(j)*trimmed_mean_cluster_j);
    end
    genevct_allgenes(:,i) = trimmed_sum_i / sum(numcells_i);
end
genevct_allgenes = genevct_allgenes(gene_names_sortinds,:);

%% Consensus gene set with AGEA coronal series, get final C and E matrices
load([glbl filesep lcl_ish filesep 'ISH_input_data.mat'],'GENname',...
    'wherecoronal','GENsetid','GENsetid_coronal','V');
uni_wb_genenames = unique(GENname(wherecoronal));
genebool = ismember(gene_names_yao,uni_wb_genenames);
gene_names = gene_names_yao(genebool);

corinds = ismember(GENsetid,GENsetid_coronal);
genlist = GENname(corinds);
voxvgene = zeros(size(V,1), length(gene_names));
for i = 1:length(gene_names)
    gene_i = gene_names{i};
    gene_i_inds = ismember(genlist,gene_i);
    voxvgene(:,i) = mean(V(:,gene_i_inds),2);
end

genevct = genevct_allgenes(genebool,:);

%% Save Yao_Inputs.mat
load([glbl filesep lcl_ish filesep 'Tasic_Inputs.mat'],'GENGDmod','listBmap',...
    'nonzerovox','structIndex','structList');
save([glbl filesep lcl_ish filesep 'Yao_Inputs.mat'],'classkey','gene_names',...
    'gene_names_yao','genevct','genevct_allgenes','listBmap','nonzerovox',...
    'structIndex','structList','voxvgene','-v7.3')
