function SFigure2c_Generator(elbowind,outstruct,ng_param_list,naround)

%This function creates panels b & c for Supplemental Figure 1. Panel
%(b) is recreated first, then panel (c). Panel (b) is the
%elbow curve for parameter selection using major cell type classes from
%Tasic, et al., 2018. Panel (c) demonstrates that nG parameter choices
%close to the chosen elbow value produce fairly similar results. 

%INPUTS: 
%fitstruct: structure with the error values per nG, for (b)
%outstruct: structure with nG values and per-voxel cell type maps, for (c)
%ng_param_list: gene set size vector, necessary for (b) & (c)
%naround: distance around the elbow to correlate with elbow maps, for (c)
%classkey: the string array of major cell type name labels from Tasic, et
%al., 2018, necessary for (c)

%THIS IS PANEL C
%calculating correlations between elbow and nearby nG cell type maps
corr_range = (elbowind-naround):(elbowind+naround);
meanrs = zeros(length(corr_range),1);
ngcomp = meanrs;
rvals = zeros(length(corr_range),size(outstruct(1).corrB,2));
for i = 1:length(corr_range)
    compmap = outstruct(i).corrB;
    rs = corr(compmap,outstruct(length(corr_range)-naround).corrB);
    rs(~(1:size(rs,1)+1:size(rs,1)*size(rs,2))) = 0;
    rs = diag(rs);
    rvals(i,:) = rs.';
    meanrs(i) = mean(squeeze(rvals(i,:)));
    ngcomp(i) = ng_param_list(corr_range(i));
end

% plotlabs = [{'Mean'} classkey];
% cmap = hsv(size(outstruct(1).corrB,2));
figure('Position',[0 0 600 400]); hold on;
plot(ngcomp,meanrs,'k-','LineWidth',3); hold on;
plot(ngcomp,rvals,'LineWidth',1); hold on; 
colormap(hsv(size(outstruct(1).corrB,2)));
% legend(plotlabs,'Location','eastoutside','FontSize',14);
set(gca,'FontSize',18);
xlim([ngcomp(1) ngcomp(end)]);
ylim([min(rvals(:))-0.05 max(rvals(:))+0.05]);
xticks([ngcomp(1) median(ngcomp) ngcomp(end)]);
yticks([min(rvals(:)) (max(rvals(:))-min(rvals(:)))/2+min(rvals(:)) max(rvals(:))]);
xlabel('{\itn}_G','FontSize',20)
ylabel('R-Val.','FontSize',20);
title('Correlation Between Elbow & Nearby {\itn}_G','FontSize',20);

end


