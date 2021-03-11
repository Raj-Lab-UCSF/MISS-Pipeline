function B = CellDensity_Corr(E_red,C_red)

% unitygenes = ismember(entrez_names,Zeisel_gene_names);
E_corr = E_red;
C_corr = C_red;
% E_corr = E_corr.'; C_corr = C_corr.';
nvox = size(E_corr,2);
logE = log2(E_corr + 1); % Zeisel et al initial log2 normalization
mn_logE = mean(logE,2); std_logE = std(logE,[],2);
Z_logE = (logE - repmat(mn_logE,1,nvox)) ./ repmat(std_logE,1,nvox);
Enans = isnan(Z_logE);
ntypes = size(C_corr,2);
mn_C = mean(C_corr,2); 
Z_C = (C_corr - repmat(mn_C,1,ntypes)) ./ repmat(mn_C,1,ntypes);
Cnans = isnan(Z_C);
allnans = logical(Cnans(:,1) + Enans(:,1)); 
Bcorr = corr(Z_C(~allnans,:),Z_logE(~allnans,:));
Bcorr(Bcorr<0) = 0;
B = Bcorr.';

end