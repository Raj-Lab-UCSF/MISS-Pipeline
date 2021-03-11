load('/Users/christophermezias/Documents/MISS_General/MatFiles/kimvols_raw.mat');
load('/Users/christophermezias/Documents/MISS_General/MatFiles/listB.mat');
listB(12,:) = [];

naninds = isnan(rawinds);
newinds = rawinds(~(naninds));
newmeta = rawmeta(~(naninds));
newvols = rawvols(~(naninds));

notreginds = strcmp('AGG',newmeta);
newinds(notreginds) = [];
newmeta(notreginds) = [];
newvols(notreginds) = [];

newinds = [newinds;46;196;211;212];
newmeta = [newmeta;listB(46);listB(196);listB(211);listB(212)];
newvols = [newvols;NaN;NaN;NaN;NaN];

[~,sortinds] = sort(newinds);

kim_vols_reorder = newvols(sortinds);
kim_labs_reorder = newmeta(sortinds);
save('/Users/christophermezias/Documents/MISS_General/MatFiles/kim_vols_listB_order.mat','kim_vols_reorder','kim_labs_reorder');

load('/Users/christophermezias/Documents/MISS_General/MatFiles/kim_totals_reorder_m.mat');
kim_dense_reorder = kim_totals_reorder ./ repmat(kim_vols_reorder,1,3);
save('/Users/christophermezias/Documents/MISS_General/MatFiles/kim_density_listB_order.mat','kim_dense_reorder');
