function [fitstruct, outstruct] = nG_ParameterFitter(voxvgene_, genevct_,...
                                    gene_names_, method_, C_indivcells_,...
                                    ct_labvec_, ng_param_list_, lambda_,...
                                    k_, crossval_, matdir_,preloadinds_)
% This function performs a brute-force parameter sweep of optionally
% supplied nG values to find one that produces a minimum loss under a
% GMM-based nearest centroid classification algorithm, as described in
% Mezias et al, 2020. The user must supply voxvgene_, genevct_,
% gene_names_, C_indivcells_, and ct_labvec_, all of which are outputs of
% upstream preprocessing functions, as well as the method_ character array
% that indicates which gene subset selection method to use. This function
% outputs two structs: fitstruct, which contains all of the information
% required to construct the loss function curves in the manuscript and
% determine the nG* value for a given lambda and sigma hyperparameter pair,
% and outstruct, which contains the matrix inversion outputs at nG*.
%
% This function calls upon various dependent functions that perform
% different aspects of the data processing and quantitative assessment.

% Default nG parameter lists
if nargin < 12
    preloadinds = [];
    if nargin < 11
        matdir_ = [cd filesep 'MatFiles'];
        if nargin < 10
            crossval_ = 1;
            if nargin < 9
                k_ = 5;
                if nargin < 8
                    lambda_ = 98;
                    if nargin < 7
                        if nargin < 6
                            error(sprintf('Error. \nUser must supply a voxel x gene matrix, a gene x cell type matrix, a gene names cell array, and a subsetting method'))
                        end
                        if ismember(method_, {'MRx3','mRMR','colAMD'})
                            ng_param_list_ = [100:70:1500, 1750:250:3500, 3855]; 
                        elseif strcmp(method_, 'DBSCAN')
                            ng_param_list_ = [0.0001, 0.0004:0.0001:0.001, 0.002:0.001:0.005, 0.015:0.01:0.195]; 
                        elseif strcmp(method_, 'Entropy')
                            ng_param_list_ = [0.25,0.5:0.125:2.75,2.8:0.05:3.25];
                        elseif strcmp(method_, 'Zeisel')
                            ng_param_list_ = [2:5, 6:2:50]; % made up something here
                        else
                            error('Error. \n%s is an incorrect subsetting method identifier',method_)
                        end
                    end    
                end           
            end
        end
    end
end

fitstruct = struct;
outstruct = struct;
nG_total = length(gene_names_);

fprintf('Initializing preloaded gene indices\n')
if isempty(preloadinds_)
    if strcmp(method_,'MRx3')
        preloadinds_ = MRx3_Selector_Prefilter(genevct_,voxvgene_,nG_total,lambda_);
    else
        preloadinds_ = [];
    end
end

for i = 1:length(ng_param_list_)
    fprintf('Determining subset, nG parameter value %d/%d\n',i,length(ng_param_list_))
    param = ng_param_list_(i);
    % Create reduced versions of voxvgene and genevct according to the
    % method and parameter specified by the user
    [~,~,nGen,~,C_ind_red] = GeneSelector_IndivCells(genevct_,voxvgene_,C_indivcells_,gene_names_,param,lambda_,method_,preloadinds_);
    
    % Calculate classification error
    tic
    fprintf('Determining GMM classification error, nG parameter value %d/%d\n',i,length(ng_param_list_))
    savegroups = 1;
    gmmstruct = GMM_Nearest_Neighbor_Posterior(C_ind_red, ct_labvec_, k_, crossval_, savegroups);
    fitstruct(i).gmmstruct = gmmstruct;
    fitstruct(i).lambda = lambda_;
    fitstruct(i).crossval = crossval_;
    fitstruct(i).nG_param = param;
    fitstruct(i).nGen = nGen;
    fitstruct(i).error = gmmstruct.error;
    fitstruct(i).accurary = gmmstruct.accuracy;
    fitstruct(i).posteriors = gmmstruct.gmmpost;
    toc
    fprintf('Done, GMM fitting, nG parameter value %d/%d\n',i,length(ng_param_list_))
end

% Elbow determination for nG range supplied
fprintf('Determining optimal nG value\n');
errors = zeros(1,length(fitstruct));
for i = 1:length(fitstruct)
   errors(i) = fitstruct(i).error;
end
norm_error = (errors - min(errors)) / (max(errors) - min(errors));
normnG = ng_param_list_/nG_total;
dist2origin = sqrt((normnG).^2 + (norm_error).^2);
[~,elbowind] = min(dist2origin);
nG_opt = fitstruct(elbowind).nGen;
nG_param_opt = fitstruct(elbowind).nG_param;
outstruct.nGen = nG_opt;

% Infer cell density per voxel in arbitrary units
tic
fprintf('Nonnegative matrix inversion at optimal nG value \n')
[E_red,C_red] = GeneSelector(genevct_,voxvgene_,gene_names_,nG_param_opt,lambda_,method_,preloadinds_);
B = CellDensityInference(E_red,C_red);
outstruct.corrB = B;
toc

if strcmp(method_,'MRx3')
    outstruct.lambda = lambda_;
end

% Sum and average over CCF regions
[sumB,meanB] = Voxel_To_Region(B,matdir_);
outstruct.Bsums = sumB; % total cells per region
outstruct.Bmeans = meanB; % mean cell count per region

end