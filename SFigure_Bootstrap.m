function SFigure_Bootstrap(outstruct_elbow,idx,outstruct_bs,savenclose,matdir)

if nargin < 5
    matdir = [cd filesep 'MatFiles'];
    if nargin < 4
        savenclose = 0;
        if nargin < 3
            load([matdir filesep 'outstruct_bootstrap.mat'],'outstruct_bootstrap');
            outstruct_bs = outstruct_bootstrap;
            if nargin < 2
                idx = 1;
                if nargin < 1
                    load([matdir filesep 'MISSMaps_2BSaved_ng606.mat'],'outstruct');
                    outstruct_elbow = outstruct;
                end
            end
        end
    end
end

elbowngs = zeros(length(outstruct_bs),1);
elbowcorrs = elbowngs;
kimcorrs_bs = zeros(length(outstruct_bs),3);
kimrhos_bs = kimcorrs_bs;
taus_bs = elbowngs;
corrB_full = outstruct_elbow(idx).corrB;
elbow_full = outstruct_elbow(idx).nGen;

for i = 1:length(outstruct_bs)
    elbowngs(i) = outstruct_bs(i).nGen;
    elbowcorrs(i) = corr(outstruct_bs(i).corrB(:),corrB_full(:));
    [~,kimRs,kimrhos] = CorrelationsCalc_Density(outstruct_bs,i,matdir,1);
    kimcorrs_bs(i,1) = kimRs.pv(1);
    kimcorrs_bs(i,2) = kimRs.sst(1);
    kimcorrs_bs(i,3) = kimRs.vip(1);
    kimrhos_bs(i,1) = kimrhos.pv(1);
    kimrhos_bs(i,2) = kimrhos.sst(1);
    kimrhos_bs(i,3) = kimrhos.vip(1);
    if ismember('tau',fieldnames(outstruct_bs))
        taus_bs(i) = outstruct_bs(i).tau;
    else
        ts = TauCalc_mod(outstruct_bs,i,{'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'},...
            9:15,[1 2 3 3 4 4 4],matdir);
        taus_bs(i) = ts.tau;
    end
end
kimcorr_full = zeros(1,3); kimrho_full = kimcorr_full;
[~,kimRs,kimrhos] = CorrelationsCalc_Density(outstruct_elbow,idx,matdir,1);
kimcorr_full(1) = kimRs.pv(1);
kimcorr_full(2) = kimRs.sst(1);
kimcorr_full(3) = kimRs.vip(1);
kimrho_full(1) = kimrhos.pv(1);
kimrho_full(2) = kimrhos.sst(1);
kimrho_full(3) = kimrhos.vip(1);
ts = TauCalc_mod(outstruct_elbow,idx,{'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'},...
    9:15,[1 2 3 3 4 4 4],matdir);
tau_full = ts.tau;
vals_full = [kimcorr_full(1),kimrho_full(1),kimcorr_full(2),kimrho_full(2),...
    kimcorr_full(3),kimrho_full(3),tau_full];

figure('Units','inches','Position',[0 0 20 6]);
subplot(1,2,1);
histogram(elbowngs,10,'FaceColor',[0 0 1],'EdgeColor',[0 0 1]); hold on;
h1 = plot([elbow_full elbow_full],[0 30],'r--','LineWidth',3);
xlabel('$n_{G}^{*}$','Interpreter','latex'); 
ylabel('\# of Simulations','Interpreter','latex');
title('Bootstrapped Elbow $n_{G}^{*}$','Interpreter','latex');
legend(h1,{sprintf('Full Dataset $n_{G}^{*} = %d$',elbow_full)},'location',...
    'northeast','Interpreter','latex');
set(gca,'FontSize',24,'TickLabelInterpreter','latex'); hold off;
subplot(1,2,2);
histogram(elbowcorrs,10,'FaceColor',[0 0.5 1],'EdgeColor',[0 0.5 1]);
xlabel('$R$','Interpreter','latex'); 
title({'Correlations to $D$ at $n_{G}^{*} = 606$','Across All Cell Types'},...
    'Interpreter','latex');
set(gca,'FontSize',24,'TickLabelInterpreter','latex');

if savenclose
    print('bootstrap_sim_histograms','-dtiffn');
    close
end

rng(0);
xlabs = {'$R_{Pvalb+}$','$\rho_{Pvalb+}$','$R_{Sst+}$','$\rho_{Sst+}$',...
    '$R_{Vip+}$','$\rho_{Vip+}$','$\tau_{adj}$'};
cmap = [[0 1 0];[0.4 1 0];[0 0 1];[0 0.4 1];[1 0 1];[1 0.4 1];[1 0.75 0.4]];
vals = [kimcorrs_bs(:,1),kimrhos_bs(:,1),kimcorrs_bs(:,2),kimrhos_bs(:,2),...
    kimcorrs_bs(:,3),kimrhos_bs(:,3),taus_bs];
xpos = zeros(size(vals,2),1);
g = xpos;
xposscatter = @(x,y) 0.2 * (2*rand(length(x),1) - 1) + y;
for j = 1:size(vals,2)
    curinds = (1:size(vals,1)) + (j-1)*size(vals,1);
    g(curinds) = j;
    xpos(curinds) = xposscatter(curinds,j);
end

figure('Units','inches','Position',[0 0 18 6]); hold on; box on;
gscatter(xpos,vals(:),g,cmap(1:size(vals,2),:),[],15,'off');
b = boxplot(vals,'Colors',cmap(1:size(vals,2),:),'Symbol','');
for k = 1:length(vals_full)
   h1 = plot([k-0.4, k+0.4],[vals_full(k),vals_full(k)],'r:','LineWidth',3); 
end
set(b,{'linew'},{3})
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gca, 'XTick', 1:size(vals,2), 'XTickLabel', xlabs);
legend(h1,{'Full Dataset'},'Interpreter','latex','Location','northeast');
title('Bootstrapped Quantitative Validation','Interpreter','latex'); 
ylabel('Validation Metric Value','Interpreter','latex');
xlabel('');
set(gca,'TickLength',[0 0],'FontSize',24,'LineWidth',0.75);

if savenclose
    print('bootstrap_sim_barplot','-dtiffn');
    close
end

% 
% % testgeneinds = geneinds_bootstrap(:,1:606);
% testoverlap = zeros(100);
% for i = 1:100
%     for j = 1:100
% %         testoverlap(i,j) = length(intersect(testgeneinds(i,:),testgeneinds(j,:)))/606;
%         testoverlap(i,j) = corr(outstruct_bs(i).corrB(:),outstruct_bs(j).corrB(:));
%     end
% end
% x = triu(ones(100),1); x = logical(x(:));
% testoverlap = testoverlap(x);
% histogram(testoverlap); xlabel('Pairwise Correlations, corrB')

% figure;
% subplot(1,3,1);
% histogram(kimcorrs_bs(:,1),'FaceColor',[0 1 0],'EdgeColor',[0 0 0]);
% xlabel('Bootstrapped Pvalb Corrs');
% subplot(1,3,2);
% histogram(kimcorrs_bs(:,2),'FaceColor',[0 0 1],'EdgeColor',[0 0 0]);
% xlabel('Bootstrapped Sst Corrs');
% subplot(1,3,3);
% histogram(kimcorrs_bs(:,3),'FaceColor',[1 0 1],'EdgeColor',[0 0 0]);
% xlabel('Bootstrapped Vip Corrs');
% figure;
% histogram(taus,'FaceColor',[1 0 0],'EdgeColor',[0 0 0]);
% xlabel('Bootstrapped Tau Reorder');




end