function [LinRvals,PeaRvals] = CorrelationsCalc_Density(outstruct,idx,directory,tord)
% Calculates Lin and Pearson correlations from a pre-created outstruct,
% with a specified index (idx) that indicates which row of the outstruct to
% calculate these metrics for.

if nargin < 4
    tord = 0;
    if nargin < 3
        directory = [cd filesep 'MatFiles'];
    end
end

% Define gross anatomical regions of interest from AIBS CCF
neoinds = 57:94; % neocortical regions in the CCF atlas
foreinds = [1:11,23:94,141:156,170:212]; % forebrain regions in the CCF atlas
wbinds = 1:212; % all 212 (bilateral) regions in the CCF atlas
indcell = {neoinds,foreinds,wbinds};

% The following function calculates the Lin's Concordance Correlation
% Coefficient.
LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

for i = 1:length(indcell)
    testinds = indcell{i};
    
    % Extract Kim et al, 2017 interneuron data from pre-created file
    if tord
        load([directory filesep 'kim_density_listB_order.mat']);
        kim_totals_reorder = kim_dense_reorder;
    else
        load([directory filesep 'kim_totals_reorder_m.mat']);
    end
    kim_totals_reorder = kim_totals_reorder(testinds,:);
    nonnaninds = zeros(size(kim_totals_reorder));
    for j = 1:size(kim_totals_reorder,2)
        nonnaninds(:,j) = ~isnan(kim_totals_reorder(:,j));
    end
    nonnaninds = logical(nonnaninds);
    kim_pv_wb = kim_totals_reorder(nonnaninds(:,1),1);
    kim_sst_wb = kim_totals_reorder(nonnaninds(:,2),2);
    kim_vip_wb = kim_totals_reorder(nonnaninds(:,3),3);

    %Extract MISS simulated interneuron density across AIBS CCF regions
    if tord
        sumcts_wb = (outstruct(idx).Bmeans(1:213,:) + outstruct(idx).Bmeans(214:end,:)) / 2;
        sumcts_wb = sumcts_wb * 125;
    else
        sumcts_wb = outstruct(idx).Bsums(1:213,:) + outstruct(idx).Bsums(214:end,:);
    end
    sumcts_wb(12,:) = [];
    sumcts_wb_mod = sumcts_wb(testinds,:);
    sumcts_pv_wb = sumcts_wb_mod(nonnaninds(:,1),3);
    sumcts_sst_wb = sumcts_wb_mod(nonnaninds(:,2),6);
    sumcts_vip_wb = sumcts_wb_mod(nonnaninds(:,3),7);
    
    %Lin R values
    LinRvals.pv(i) = LinRcalc(kim_pv_wb,sumcts_pv_wb);
    LinRvals.sst(i) = LinRcalc(kim_sst_wb,sumcts_sst_wb);
    LinRvals.vip(i) = LinRcalc(kim_vip_wb,sumcts_vip_wb);
    
    %Pearson R values
    PeaRvals.pv(i) = corr(kim_pv_wb,sumcts_pv_wb);
    PeaRvals.sst(i) = corr(kim_sst_wb,sumcts_sst_wb);
    PeaRvals.vip(i) = corr(kim_vip_wb,sumcts_vip_wb);
end

end