function MISS_Brainframe(outstruct,types,datsetname,elbowind,xfac,savenclose,voxthresh,cmap_range,img_name,matdir)

if nargin < 9
    matdir = cd;
    if nargin < 8
        img_name = 'invNsubs';
        if nargin < 7
            cmap_range = [0 0 0;1 1 1];
            if nargin < 6
                voxthresh = 1;
                if nargin < 5
                    savenclose = 0;
                    if nargin < 4
                        xfac = 1;
                        if nargin < 3
                            elbowind = 582;
                        end
                    end
                end
            end
        end
    end
end

load([matdir filesep 'default_mouse.mat'],'input_struct');
if strcmp(datsetname,'Tasic')
    load([matdir filesep 'Tasic_Inputs.mat'],'GENGDmod','structList','structIndex','nonzerovox','classkey');
elseif strcmp(datsetname,'Zeisel')
    load([matdir filesep 'Zeisel_Inputs.mat'],'GENGDmod','structList','structIndex','nonzerovox','classkey');
elseif strcmp(datsetname,'Yao')
    load([matdir filesep 'Yao_Inputs.mat'],'GENGDmod','structList','structIndex','nonzerovox','classkey');
end

input_struct.nbin = 10;
input_struct.voxUreg = 0;
input_struct.xfac = xfac;
input_struct.pointsize = 0.001;
input_struct.savenclose = savenclose;
input_struct.img_format = 'tiff';
input_struct.regsUbins = 0;

ng_plotted = outstruct(elbowind).nGen;

for i = 1:length(types)
    
    cmap = twocolor(cmap_range{i}(1,:),cmap_range{i}(2,:),input_struct.nbin);
    input_struct.cmap = cmap;
    newVoxMap = zeros(size(GENGDmod));
    
    curinput = outstruct(elbowind).corrB(:,types(i));
    
    for k = 1:length(structList)
        [~,voxinds] = ismember(structIndex{k},nonzerovox);
        allvox = nonzerovox(voxinds);
        curvox = curinput(voxinds);
        newVoxMap(allvox) = curvox;
    end
    
    datinput = newVoxMap;
    datinput = imresize3(datinput,[133 81 115]);
    datinput(datinput < 0) = 0;
    input_struct.data = datinput;
    input_struct.img_labels = [img_name classkey{types(i)}];
    input_struct.voxthresh = voxthresh(i);
    brainframe(input_struct);
    
end



        

