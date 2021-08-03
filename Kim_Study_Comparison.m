function Kim_Study_Comparison(outstruct,idx,region,savenclose,directory)

if nargin < 5
    directory = [cd filesep 'MatFiles'];
    if nargin < 4
        savenclose = 0;
        if nargin < 3
            region = 'neo';
        end
    end
end

sampleinds = [70,92,186]; % Secondary motor area, primary visual area, dorsal lateral geniculate
switch region
    case 'fore'
        reginds = [1:11,23:56,141:156,170:212];
    case 'hind'
        reginds = [12:22 95:140 157:169];
    case 'neo'
        reginds = 57:94;
%         reginds = [57:65,67,69:94];
    case 'neo+fore'
        reginds = [1:11,23:94,141:156,170:212];
    otherwise
        reginds = 1:212;
end
sampleinds = intersect(sampleinds,reginds);
typenames = {'Pvalb','Sst','Vip'};
typecolors = {'g','b','m'};
typesymbol = {'o','s','d'};
kimdata = struct;
testdata = struct;

load([directory filesep 'kim_density_listB_order.mat'],'kim_dense_reorder');
load([directory filesep 'Tasic_Inputs.mat'],'classkey');
kim_density_reg = kim_dense_reorder(reginds,:);
kimmaxes = nanmax(kim_density_reg);
nonnaninds_reg = zeros(size(kim_density_reg));
for j = 1:size(kim_density_reg,2)
    nonnaninds_reg(:,j) = ~isnan(kim_density_reg(:,j));
end
nonnaninds_reg = logical(nonnaninds_reg);
for j = 1:size(kim_density_reg,2)
    kimdata.region.(typenames{j}) = kim_density_reg(nonnaninds_reg(:,j),j)/kimmaxes(j);
end

kim_density_sample = kim_dense_reorder(sampleinds,:);
nonnaninds_sample = zeros(size(kim_density_sample));
for j = 1:size(kim_density_sample,2)
    nonnaninds_sample(:,j) = logical(~isnan(kim_density_sample(:,j)));
end
nonnaninds_sample = logical(nonnaninds_sample);
for j = 1:size(kim_density_reg,2)
    kimdata.sample.(typenames{j}) = kim_density_sample(nonnaninds_sample(:,j),j)/kimmaxes(j);
end

meancts_wb = (outstruct(idx).Bmeans(1:213,:) + outstruct(idx).Bmeans(214:end,:))/2;
meancts_wb(12,:) = []; % not in ISH parcellation but in CCF
meancts_reg = meancts_wb(reginds,:);
testmaxes = zeros(1,length(typenames));
for j = 1:length(typenames)
    testdensity = meancts_reg(nonnaninds_reg(:,j),ismember(classkey,typenames{j}));
    testmaxes(j) = max(testdensity);
    testdata.region.(typenames{j}) = testdensity/testmaxes(j);
end

meancts_sample = meancts_wb(sampleinds,:);
for j = 1:length(typenames)
    testdensity = meancts_sample(nonnaninds_sample(:,j),ismember(classkey,typenames{j}));
    testdata.sample.(typenames{j}) = testdensity/testmaxes(j);
end

figure('Units','inch','Position',[0 0 16 4]);
for j = 1:length(typenames)
    subplot(1,3,j); hold on;
    set(gca,'FontSize',12);
    scatter(testdata.region.(typenames{j}),kimdata.region.(typenames{j}),60,[typecolors{j} typesymbol{j}],'filled'); hold on;
    scatter(testdata.sample.(typenames{j}),kimdata.sample.(typenames{j}),200,'rp','filled'); hold on;
    [pearR,pearp] = corr(testdata.region.(typenames{j}),kimdata.region.(typenames{j}));
    [spearR,spearp] = corr(testdata.region.(typenames{j}),kimdata.region.(typenames{j}),'type','Spearman');
    title(typenames{j},'FontSize',20)
    if pearp < 0.001
        txt_R = sprintf('R = %.2f***',pearR);
    elseif pearp < 0.01
        txt_R = sprintf('R = %.2f**',pearR);
    elseif pearp < 0.05
        txt_R = sprintf('R = %.2f*',pearR);
    else
        txt_R = sprintf('R = %.2f',pearR);
    end
    if spearp < 0.001
        txt_rho = sprintf([char(961) ' = %.2f***'],spearR);
    elseif spearp < 0.01
        txt_rho = sprintf([char(961) ' = %.2f**'],spearR);
    elseif spearp < 0.05
        txt_rho = sprintf([char(961) ' = %.2f*'],spearR);
    else
        txt_rho = sprintf([char(961) ' = %.2f'],spearR);
    end
    txt = {txt_R, txt_rho};
    text(0.625*(max(testdata.region.(typenames{j}))-min(testdata.region.(typenames{j})))+min(testdata.region.(typenames{j})),...
        0.15*(max(kimdata.region.(typenames{j}))-min(kimdata.region.(typenames{j})))+min(kimdata.region.(typenames{j})),...
        txt,'FontSize',18);
    if j == 1
        ylabel('Empirical Density','FontSize',20);
    end
    xlabel('Inferred Density','FontSize',20);
    set(gca,'YTick',[min(kimdata.region.(typenames{j})), (1+min(kimdata.region.(typenames{j})))/2, 1]);
    ytickformat('%.2f'); xtickformat('%.2f'); 
    set(gca,'XTick',[min(testdata.region.(typenames{j})), (1+min(testdata.region.(typenames{j})))/2, 1]);
    set(gca,'FontSize',18);
    p = polyfit(testdata.region.(typenames{j}),kimdata.region.(typenames{j}),1);
    x_maxy = (1 - p(2))/p(1);
    if x_maxy < 1
        plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
    else
        plot([0,1],[p(2),(p(2) + p(1))],'Color',[0.9 0.3 0.1]); hold on;
    end
    xlim([min(testdata.region.(typenames{j})) 1]);
    ylim([min(kimdata.region.(typenames{j})) 1]);
    box on
end

if savenclose
    print('kim_study_comparison','-dtiffn');
    close
end
end