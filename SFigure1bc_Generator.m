function SFigure1bc_Generator(fitstruct,outstruct,naround,classkey)

%This function creates panels b & c for Supplemental Figure 1. Panel
%(b) is recreate first, then panel (c), then panel (d). Panel (b) is the
%elbow curve for parameter selection using major cell type classes from
%Tasic, et al., 2018. Panel (c) demonstrates that nG parameter choices
%close to the chosen elbow value produce fairly similar results. Panel (d)
%demonstrates that there is no bias in residual or square error between E
%and C*B to being lower in scrnaSeq sampled regions.

%INPUTS: 
%fitstruct: structure with the error values per nG, for (b)
%outstruct: structure with nG values and per-voxel cell type maps, for (b)
%& (c)
%naround: distance around the elbow to correlate with elbow maps, for (c)
%classkey: the string array of major cell type name labels from Tasic, et
%al., 2018, necessary for (c)

%THIS IS PANEL B
%Generating an error vs. nG curve for parameter elbow value selection
%These vectors are normalized & 3855 is the total n-genes in the set usead
for i = 1:length(fitstruct)
    error(i) = fitstruct(i).gmmstruct.error;
    ng_param_list(i) = outstruct(i).nGen;
end

norm_error = (error - min(error)) / (max(error) - min(error));
normnG = ng_param_list / 3855;
dist2origin = sqrt((normnG-0).^2 + (norm_error-0).^2);
[~,elbowind] = min(dist2origin);

figure('Position',[0 0 400 600]); hold on;
plot(normnG,norm_error,'k-','LineWidth',2.5); hold on;
plot([outstruct(elbowind).nGen outstruct(elbowind).nGen],...
    [0 max(norm_error)],'r--','LineWidth',2.5); hold on;
xlim([0 max(normnG)]);
ylim([0 max(norm_error)]);
set(gca,'FontSize',18);
xticks([0 500/3855 1000/3855]);
xticklabels({'0','500','1000'});
yticks([0 0.025 0.05 0.075 0.1]);
xlabel('\it{n}_G','FontSize',20)
ylabel('Error','FontSize',20);
title('Norm. Error Vs. \it{n}_G','FontSize',20);
text(outstruct(elbowind).nGen+2,0.5*max(norm_error),...
    ['\it{n}_G = ' num2str(outstruct(elbowind).nGen)],'FontSize',16);

%THIS IS PANEL C
%calculating correlations between elbow and nearby nG cell type maps
corr_range = (elbowind-naround):(elbowind+naround);
meanrs = zeros(length(corr_range),1);
ngcomp = meanrs;
rvals = zeros(length(corr_range),size(outstruct(1).Bvals,2));
for i = corr_range
    compmap = outstruct(i).Bvals;
    rvals(i,:) = corr(compmap,outstruct(elbowind).Bvals);
    meanrs(i) = mean(squeeze(rvals(i,:)));
    ngcomp(i) = outstruct(corr_range(i)).nGen;
end

plotlabs = ['Mean';classkey];
cmap = hsv(size(outstruct(1).Bvals,2));
figure('Position',[0 0 400 600]); hold on;
plot(ngcomp,meanrs,'k-','LineWidth',3); hold on;
plot(ngcomp,rvals,'Color',cmap,'LineWidth',1); hold on;
legend(plotlabs,'Location','eastoutside','FontSize',14);
xlim([ngcomp(1) ngcomp(end)]);
ylim([min(rvals(:))-0.05 max(rvals(:))+0.05]);
xticks([ngcomp(1) median(ngcomp) ngcomp(end)]);
yticks([min(rvals(:)) (min(rvals(:))+max(rvals(:)))/2 max(rvals(:))]);
xlabel('\it{n}_G','FontSize',20)
ylabel('R-Val.','FontSize',20);
title('Correlation Between Elbow & Nearby \it{n}_G','FontSize',20);

end


