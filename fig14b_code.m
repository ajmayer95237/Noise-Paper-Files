M = load('../data/raw/quantile/numbersteps50.011162.mat');
N = load('../data/raw/otsu/meanConcenmeanSpots.mat');
S = load('../data/raw/quantile/meanConcenmeanSpots.mat');


dataSet = S;
mRNA = dataSet.N_all;

nonzero_ind = [];


cln3_822 = [];
cln3_822_it = 1;
cln3_822pk = [];
cln3_822pk_it = 1;
cln3_882 = [];
cln3_882_it = 1;
cln3_882pk = [];
cln3_882pk_it = 1;
cln3_909 = [];
cln3_909_it = 1;
cln3_909pk = [];
cln3_909pk_it = 1;

bni1_822 = [];
bni1_822_it = 1;
bni1_822pk = [];
bni1_822pk_it = 1;
bni1_882 = [];
bni1_882_it = 1;
bni1_882pk = [];
bni1_882pk_it = 1;
bni1_909 = [];
bni1_909_it = 1;
bni1_909pk = [];
bni1_909pk_it = 1;

whi3_822 = [];
whi3_822_it = 1;
whi3_822pk = [];
whi3_822pk_it = 1;
whi3_882 = [];
whi3_882_it = 1;
whi3_882pk = [];
whi3_882pk_it = 1;
whi3_909 = [];
whi3_909_it = 1;
whi3_909pk = [];
whi3_909pk_it = 1;

for i = 1:length(mRNA)
    sphere_i = M.which_condition_sphere(i,:);
    probe_i = sphere_i(3);
    strain_i = sphere_i(4);
    pk_i = sphere_i(5);
    hypha_i = sphere_i(6);
    if probe_i == 1
       if strain_i == 822
          if pk_i == 1
             cln3_822(cln3_822_it,:) = [i,hypha_i];
             cln3_822_it = cln3_822_it+1;
          else
             cln3_822pk(cln3_822pk_it,:) = [i,hypha_i];
             cln3_822pk_it = cln3_822pk_it+1;
          end
       elseif strain_i == 882
          if pk_i == 1
             cln3_882(cln3_882_it,:) = [i,hypha_i];
             cln3_882_it = cln3_882_it+1;
          else
             cln3_882pk(cln3_882pk_it,:) = [i,hypha_i];
             cln3_882pk_it = cln3_882pk_it+1;
          end  
       else
          if pk_i == 1
             cln3_909(cln3_909_it,:) = [i,hypha_i];
             cln3_909_it = cln3_909_it+1;
          else
             cln3_909pk(cln3_909pk_it,:) = [i,hypha_i];
             cln3_909pk_it = cln3_909pk_it+1;
          end            
       end
    elseif probe_i == 2
       if strain_i == 822
          if pk_i == 1
             bni1_822(bni1_822_it,:) = [i,hypha_i];
             bni1_822_it = bni1_822_it+1;
          else
             bni1_822pk(bni1_822pk_it,:) = [i,hypha_i];
             bni1_822pk_it = bni1_822pk_it+1;
          end
       elseif strain_i == 882
          if pk_i == 1
             bni1_882(bni1_882_it,:) = [i,hypha_i];
             bni1_882_it = bni1_882_it+1;
          else
             bni1_882pk(bni1_882pk_it,:) = [i,hypha_i];
             bni1_882pk_it = bni1_882pk_it+1;
          end  
       else
          if pk_i == 1
             bni1_909(bni1_909_it,:) = [i,hypha_i];
             bni1_909_it = bni1_909_it+1;
          else
             bni1_909pk(bni1_909pk_it,:) = [i,hypha_i];
             bni1_909pk_it = bni1_909pk_it+1;
          end            
       end 
    elseif probe_i == 3
       if strain_i == 822
          if pk_i == 1
             whi3_822(whi3_822_it,:) = [i,hypha_i];
             whi3_822_it = whi3_822_it+1;
          else
             whi3_822pk(whi3_822pk_it,:) = [i,hypha_i];
             whi3_822pk_it = whi3_822pk_it+1;
          end
       elseif strain_i == 882
          if pk_i == 1
             whi3_882(whi3_882_it,:) = [i,hypha_i];
             whi3_882_it = whi3_882_it+1;
          else
             whi3_882pk(whi3_882pk_it,:) = [i,hypha_i];
             whi3_882pk_it = whi3_882pk_it+1;
          end  
       else
          if pk_i == 1
             whi3_909(whi3_909_it,:) = [i,hypha_i];
             whi3_909_it = whi3_909_it+1;
          else
             whi3_909pk(whi3_909pk_it,:) = [i,hypha_i];
             whi3_909pk_it = whi3_909pk_it+1;
          end            
       end 
    end

end
data_cln3 = cln3_822;
data_cln3_all = [cln3_822; cln3_822pk; cln3_882pk; cln3_882; cln3_909pk; cln3_909];
data_bni1 = bni1_882;
data_bni1_882_all = [bni1_882pk; bni1_882];
data_bni1_all = [bni1_822; bni1_822pk; bni1_882pk; bni1_882; bni1_909pk; bni1_909];
data_whi3 = whi3_822pk;
data_whi3_all = [whi3_822; whi3_822pk; whi3_882; whi3_882pk; whi3_909];

data = data_cln3;
hypha = data(:,2);
numhypha = length(unique(hypha));
dataNall = dataSet.N_all(data);
dataVol = dataSet.spherevol(data);
dataSpot = dataSet.meanSpotWeights(data);
dataSpotExpected = dataSet.meanSpotWeightsExpected;
minNall = min(dataNall);
maxNall = max(dataNall);

notNAind = ~isnan(dataNall(:,1));
dataNall=dataNall(notNAind,:);
dataSpot = dataSpot(notNAind,:);
dataSpotExpected = dataSpotExpected(notNAind,:);
dataVol = dataVol(notNAind,:);

nonzeroVols = find(dataVol(:,1) > 0);
concens = nonzeros(dataNall(nonzeroVols,1)./dataVol(nonzeroVols,1));

nonzeromeanSpots = find(dataSpot(nonzeroVols,1) > 0); 
meanSpot = nonzeros(dataSpot(nonzeroVols,1));
meanSpotExpected = dataSpotExpected(nonzeromeanSpots,1);


figure(4)
set(groot,'defaultAxesTickLabelInterpreter','latex');  

bins=logspace(-2,3,40); [bincounts,whichbin] = histc(concens,bins);
for i_bin = 1:length(bins)-1
    is_inbin = whichbin == i_bin; 
    medianspot(i_bin) = median(meanSpot(is_inbin));
    medianspot2(i_bin) = median(meanSpotExpected(is_inbin));
end
binsc = 0.5*(bins(1:end-1)+bins(2:end));
semilogx(binsc,medianspot,'LineWidth',2);
hold on
semilogx(binsc,medianspot2,'LineWidth',2);
scatter(concens,meanSpot,3)

title("CLN3 RNA counts in 2.5 $\mu m$ nuclear neighborhoods",'interpreter','latex')
xlabel("Number of RNA in neighborhood / neighborhood volume",'interpreter','latex')
ylabel({"Mean number of RNA detected";  "per spot in neighborhood "},'interpreter','latex')
 
