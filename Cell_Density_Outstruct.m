function outstruct = Cell_Density_Outstruct(genevct_,voxvgene_,...
                                    gene_names_,nG_param_list_,lambda_,...
                                    missmethod_,infmethod_,preloadinds_,matdir_)
if nargin < 7
    infmethod_ = 'inversion';
    if nargin < 8
        preloadinds_ = [];
        if nargin < 9
            matdir_ = [cd filesep 'MISS-MatFiles'];
        end
    end
end

outstruct = struct;
for i = 1:length(nG_param_list_)
    tic;
    [E_red,C_red,nGen] = GeneSelector(genevct_,voxvgene_,gene_names_,nG_param_list_(i),lambda_,missmethod_,preloadinds_);
    if strcmp(infmethod_,'inversion')
        B = CellDensityInference(E_red,C_red);
    elseif strcmp(infmethod_,'corr')
        B = CellDensity_Corr(E_red,C_red);
    elseif strcmp(infmethod_,'inv+res')
        [B,mnresnorm_,fronorm_] = CellDensityInference_resnorm(E_red,C_red);
        outstruct(i).resnorm = mnresnorm_;
        outstruct(i).fronorm = fronorm_;
    end
    outstruct(i).corrB = B;
%     outstruct(i).nG_param = nG_param_list_(i);
    outstruct(i).nGen = nGen;
    if strcmp(missmethod_,'MRx3')
        outstruct(i).lambda = lambda_;
    end
    [sumB,meanB] = Voxel_To_Region(B,matdir_);
    outstruct(i).Bsums = sumB; % total cells per region
    outstruct(i).Bmeans = meanB; % mean cell count per region
    toc;
end
end
