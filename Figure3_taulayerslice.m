function Figure3_taulayerslice(outstruct,idx,mapmethod,slicelocs,savenclose,directory)

if nargin < 6
    directory = [cd filesep 'MatFiles'];
    if nargin < 5
        savenclose = 0;
        if nargin < 4
            slicelocs = [24,32];
            if nargin < 3
                mapmethod = 'Inv & SubS';
            end
        end
    end
end
load([directory filesep 'tau_calc_dependencies.mat'],'GENGDmod',...
    'nonzerovox','structList','structIndex','neomask3d','neobounds3d');
load([directory filesep 'input_struct_voxelrender.mat'],'input_struct');

ranks = [1 2 3 3 4 4 4];
cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
cell_inds = 9:15;
titlecellnames = {'L2/3-IT', 'L4','L5-IT', 'L5-PT', 'L6-CT','L6-IT','L6b'};
taustruct = TauCalc_mod(outstruct,idx,cell_names,cell_inds,ranks,directory);
B = outstruct(idx).corrB;
% B = outstruct(idx).Bvals;
ngen = outstruct(idx).nGen;
if ismember('method',fieldnames(outstruct(idx)))
    method = [outstruct(idx).method ' '];
else
    method = '';
end

newVoxMap = zeros(size(GENGDmod));
pervoxglut = struct;
for i = 1:length(cell_names)
    pervoxglut.(cell_names{i}) = newVoxMap;
end
newVoxMap(nonzerovox) = 1;
glutnames = cell_names;

for i = 1:length(structList)
    [~,voxinds] = ismember(structIndex{i},nonzerovox);
    allvox = nonzerovox(voxinds);
    curvox = B(voxinds,cell_inds);
    for k = 1:length(glutnames)
        curglut = squeeze(curvox(:,k));
        pervoxglut.(glutnames{k})(allvox) = curglut;
        clear curglut
    end
    clear allvox curvox voxinds
end

[d1,~,~] = find(neomask3d);
d1inds = unique(d1);
maplocs = d1inds;
for k = 1:length(glutnames)
    curmap = pervoxglut.(glutnames{k});
    for j = 10:length(maplocs)
        curloc = maplocs(j);
        ctxsurf = squeeze(neobounds3d(curloc,:,:));
        curneomask = squeeze(neomask3d(curloc,:,:));
        [ctxsurfx,ctxsurfy] = find(ctxsurf);
        curslc = squeeze(curmap(curloc,:,:));
        mnzs = min(nonzeros(curslc(:)));
        if isempty(mnzs)
            curslc(curslc==0) = 0;
        else
            curslc(curslc==0) = mnzs;
        end
        grayslc = mat2gray(curslc);
        binslc_whole = logical(grayslc);
        paleoinds = find(~curneomask);
        binslc = binslc_whole;
        binslc(paleoinds) = 0;
        skeleton = bwskel(binslc);
        [skelx,skely] = find(skeleton);
        for s = 1:length(find(skeleton))
            curx = skelx(s);
            cury = skely(s);
            dists(s) = min(sqrt((ctxsurfx-curx).^2 + (ctxsurfy-cury).^2));
        end
        if isempty(skelx) || isempty(curx) || isempty(dists)
            mndists(j-9,k) = -1;
        else
            mndists(j-9,k) = mean(dists(:));
        end
        clear dists cury curx
    end
end

lranks = 0;
trueranks = 0;
for j = 1:size(mndists,1)
    curmns = mndists(j,:);
    nosig = (curmns==-1);
    n_nosig(j) = length(find(nosig));
    curmns(nosig) = [];
    [~,rankinds] = sort(curmns,'ascend');
    curranks = ranks;
    curranks(nosig) = [];
    sliceranks = curranks;
    trueranks = [trueranks curranks];
    curranks = curranks(rankinds);
    lranks = [lranks curranks];
    if isempty(curranks)
        tau_per = 0;
    else
        [tau_per,~] = corr(curranks.', sliceranks.', 'type', 'Kendall');
    end
    tau_slice(j) = (1-(n_nosig(j)/7))*tau_per;
    sqerrvec(j) = (1 - tau_per);
end

f1 = figure; hold on;
set(f1,'Position',[0 0 250 850]); hold on; %this to set the size
subplot(15,2,3,'Parent',f1);
sgtitle([mapmethod ', \tau = ' sprintf('%.2f',taustruct.tau)],'FontSize',18);

redcellnames = glutnames;
slice_name = {'Rostral', 'Caudal'};
maplocs = slicelocs*2-1;

newVoxMap = input_struct.brain_atlas;
newVoxMap = double(logical(newVoxMap));

index = 3; index2 = 5;
for k = 1:length(redcellnames)
    curmap = pervoxglut.(redcellnames{k});
    curmap = imresize3(curmap,[133 81 115]);
    curmap(curmap<0) = 0;
%     curmax = max(max(max(curmap)));
    curpctile = prctile(nonzeros(curmap),97.5);
    for j = 1:length(maplocs)
        curloc = maplocs(j);
        slice_raw = squeeze(curmap(curloc,:,:));
        im = slice_raw;
        se = offsetstrel('ball',3,1,4);
        im = imdilate(im,se);
        im = imerode(im,se);
        im = interpn(im,2,'spline');
        bw = squeeze(newVoxMap(curloc,:,:));
        im_ = imdilate(bw,se);
        im_ = imerode(im_,se);
        bim_ = interpn(im_,2,'spline');
        bim_(bim_ < 0.5) = 0;
        bim_(bim_ >= 0.5) = 1;
        slice_final = im .* bim_;
        bwbounds = bwboundaries(bim_);
        plotmaxes = zeros(length(bwbounds),2); plotmins = plotmaxes;
        for m = 1:length(bwbounds)
            plotmaxes(m,:) = max(bwbounds{m});
            plotmins(m,:) = min(bwbounds{m});
        end
        plotmaxes = max(plotmaxes,[],1);
        plotmins = min(plotmins,[],1);
    
        figure(f1); subplot(15,2,[index index2]); hold on;
%         imagesc(slice_final,[0 0.65*curmax+eps]); hold on;
        imagesc(slice_final,[0 curpctile]); hold on;
        colormap(flipud(pink)); hold on;
        for m = 1:length(bwbounds)
            boundary = bwbounds{m};
            plot(boundary(:,2),boundary(:,1),'k','LineWidth',0.5); hold on;
        end
        ylim([plotmins(1)-10, plotmaxes(1)+10]); xlim([plotmins(2)-10, plotmaxes(2)+10]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'Ydir','reverse')
        box on;
        set(gca,'BoxStyle','full');
        set(gca,'XAxisLocation','origin')
        if index == 3
            ylabel(titlecellnames{k},'FontSize',14);
            title(slice_name{j},'FontSize',14);
        elseif index == 4
            title(slice_name{j},'FontSize',14);
        elseif index == 7 || index == 11 || index == 15
            ylabel(titlecellnames{k},'FontSize',14);
        elseif index == 19 || index == 23 || index == 27
            ylabel(titlecellnames{k},'FontSize',14);
        end
        index = index + 1;
        index2 = index2 + 1;
        clear slice_raw slice_final
    end
    index = index + 2;
    index2 = index2 + 2;
end

if savenclose
    print([mapmethod '_' 'Figure3_layertype'],'-dtiff');
    close
end
end