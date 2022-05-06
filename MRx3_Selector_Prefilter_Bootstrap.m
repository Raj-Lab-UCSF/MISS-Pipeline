function geneinds_mat = MRx3_Selector_Prefilter_Bootstrap(C_raw,E,n,lambda,...
                        resample_rate,niters,useParallel,n_cores)
% This function calculates the best n genes (rows) of matrix C according to
% the mRMR criterion with projection-error-based prefiltering of genes,
% using bootstrapping following the prefiltering step

if nargin < 8
    n_cores = 4;
    if nargin < 7
        useParallel = 0;
    end
end

rng(0);
geneinds_mat = zeros(niters,n);
remaininds = 1:size(C_raw,1);
if n > length(remaininds)
    n = length(remaininds);
end

% Create a column-normalized (i.e. by cell type) C matrix for the purpose
% of calculating the approximate F-statistic and pairwise correlations for
% the MRx3 criterion
ctmean = mean(C_raw,1);
ctmean = repmat(ctmean,size(C_raw,1),1);
ctnorm = C_raw ./ ctmean;
ctnorm(ctnorm<(0.1*std(nonzeros(ctnorm)))) = 0; % filter out excessively small signal
C_cellnorm = ctnorm;

% Create a row-normalized (i.e. by gene) C matrix for the purpose of
% calculating the rank-1 approximation to the projection error for the mRMR
% criterion; also indicate which genes have zero expression across all cell
% types
gsum = sum(C_raw,2);
zerosuminds = find(gsum == 0);
remaininds(zerosuminds) = []; %#ok<*FNDSB>

gsum = repmat(gsum,1,size(C_raw,2));
gnorm = C_raw ./ gsum;
gnorm(isnan(gnorm)) = 0;
C_genenorm = gnorm;

% Normalize E (should be redundant)
E = E.';
esum = sum(E,2);
esum = repmat(esum,1,size(E,2));
E = E ./ esum;

% Prefilter genes based on projection error and lambda
fprintf('Prefiltering based on projection error\n');
projerr = Inf(1,size(C_raw,1));
if ~useParallel
    for i = 1:size(C_raw,1)
    %     fprintf('Calculating projection error for gene %d/%d \n',i,size(C_raw,1))
        if ismember(i,remaininds)
            g_g = C_genenorm(i,:);
            E_i = E(i,:);
            noti_inds = setdiff(remaininds,i);
            pinv_C_noti = pinv(C_genenorm(noti_inds,:));
            E_g_noti = E(noti_inds,:);
            D_noti = pinv_C_noti * E_g_noti; 
            err1_i = sum((E_i - g_g*D_noti).^2);
            err2_i = sum((g_g*pinv_C_noti).^2)*dot(E_i,E_i);
            projerr(i) = err1_i + err2_i;
        end
    end
else
    parpool(n_cores);
    parfor i = 1:size(C_raw,1)
        if ismember(i,remaininds)
            g_g = C_genenorm(i,:);
            E_i = E(i,:);
            noti_inds = setdiff(remaininds,i);
            pinv_C_noti = pinv(C_genenorm(noti_inds,:));
            E_g_noti = E(noti_inds,:);
            D_noti = pinv_C_noti * E_g_noti; 
            err1_i = sum((E_i - g_g*D_noti).^2);
            err2_i = sum((g_g*pinv_C_noti).^2)*dot(E_i,E_i);
            projerr(i) = err1_i + err2_i;
        end
    end 
    delete(gcp('nocreate'));
end

error_inds = find(projerr > prctile(projerr,lambda));
remaininds = setdiff(remaininds,error_inds);
excludeinds = setdiff(1:size(C_raw,1),remaininds);
nonzeroexcludeinds = setdiff(excludeinds,zerosuminds);

Fcalc = @(g) sum((g - mean(g)).^2)/(length(g) - 1);
Rmat = corr(C_cellnorm.');

for j = 1:niters
    fprintf('Bootstrap Iteration %d/%d\n',j,niters);
    remaininds_res = remaininds;
    nonzeroexcludeinds_res = nonzeroexcludeinds;
    zerosuminds_res = zerosuminds;
    bootstrapexclude_size = round((1-resample_rate)*length(remaininds_res));
    bootstrapexcludeinds = randperm(length(remaininds_res));
    bootstrapexcludeinds = sort(remaininds_res(bootstrapexcludeinds(1:bootstrapexclude_size)));
%     display(bootstrapexcludeinds)
%     display(remaininds_res)
    remaininds_res = remaininds_res(~ismember(remaininds_res,bootstrapexcludeinds));
    % mRMR Initialization
    Vinit = zeros(1,length(remaininds_res));
    C_init = C_cellnorm(remaininds_res,:);
    for i = 1:length(Vinit)
        g_ct = C_init(i,:);
        Vinit(i) = Fcalc(g_ct);
    end
    [~,maxind] = max(Vinit);
    % geneinds = [geneinds,remaininds(maxind)];
    geneinds = remaininds_res(maxind);
    remaininds_res(maxind) = [];
    % S = C(geneinds,:);
    T = C_cellnorm(remaininds_res,:);

    % mRMR Propagation/Termination
    while length(geneinds) < n
        if isempty(remaininds_res)
            if ~isempty(bootstrapexcludeinds)
                geneinds = [geneinds, bootstrapexcludeinds(1)];
                bootstrapexcludeinds(1) = [];
            elseif ~isempty(nonzeroexcludeinds_res)
                geneinds = [geneinds, nonzeroexcludeinds_res(1)];
                nonzeroexcludeinds_res(1) = [];
            elseif ~isempty(zerosuminds_res)
                geneinds = [geneinds, zerosuminds_res(1)];
                zerosuminds_res(1) = [];
            end
        end
        szT = size(T,1);
        Vtest = zeros(1,szT);
        for i = 1:szT
            g = T(i,:);
            ind = remaininds_res(i);
            Rs = abs(Rmat(ind,geneinds));
            wc = sum(Rs)/length(geneinds);
            Vtest(i) = Fcalc(g) / wc;
        end
        [~,maxind] = max(Vtest);
    %     Vmax = [Vmax,Vm];
        geneinds = [geneinds,remaininds_res(maxind)];
        remaininds_res(maxind) = [];
    %     S = C(geneinds,:);
        T = C_cellnorm(remaininds_res,:);
    end
    geneinds_mat(j,:) = geneinds;
    clear geneinds
end
end