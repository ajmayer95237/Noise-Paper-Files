M = load('../data/raw/quantile/numbersteps50.011162.mat');
N = load('../data/raw/quantile/meanConcenmeanSpots.mat');
S = load('../data/raw/quantile/meanConcenmeanSpots.mat');
mRNA = M.N_all;

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

data = data_bni1_all;
hypha = data(:,2);
numhypha = length(unique(hypha));
%dataCV = M.CVsquared(data);
dataNall = S.N_all(data);
%dataEntropy = M.entropyballs(data);
dataVol = M.spherevol(data);
%dataconcen = M.meanconcens(data);
dataSpot = S.meanSpotWeights(data);
minNall = min(dataNall);
maxNall = max(dataNall);
meanCV = [];
meanEntropy = [];

figure(3)
notNAind = ~isnan(dataNall(:,1));
dataNall=dataNall(notNAind,:);
dataSpot = dataSpot(notNAind,:);
dataVol = dataVol(notNAind,:);

nonzeroVols = find(dataVol(:,1) > 0);
concens = nonzeros(dataNall(nonzeroVols,1)./dataVol(nonzeroVols,1));
meanSpot = nonzeros(dataSpot(nonzeroVols,1));

scatter(concens,meanSpot)
title("Neighborhoods across all hyphae")
xlabel("Number of RNA / Neighborhood volume")
ylabel("Mean spot size")
set(gca,'xscale','log')
hold on



nthresh = 1000;
threshs = linspace(min(concens),max(concens),nthresh);
slopes = (.01:.01:10);
f = zeros(length(meanSpot),1);
fdiff = norm(f-meanSpot,1);
slope_opt = 0;
thresh_opt = 0;
[x,xorder] = sort(concens);
y = meanSpot(xorder);

for thresh = threshs
    for slope = slopes
        ftemp = double(x < thresh);
        thresh_ind = sum(ftemp);
        
        for i = (thresh_ind+1):length(x)
            ftemp(i) = 1+slope*log(x(i)/thresh);          
            %ftemp(i) = 1+slope*(x(i)-thresh);

        end
        %plot(x,ftemp)
        %hold on
        
        diff = norm(ftemp-y,1);
        if diff < fdiff
           slope_opt = slope;
           thresh_opt = thresh;
           fdiff = diff;
           f = ftemp;
        end
        
            
    end    
end
plot(x,f)
title(sprintf('Optimal threshold is %f mRNA', mean(dataVol(nonzeroVols,1))*thresh_opt))

figure(4)
set(groot,'defaultAxesTickLabelInterpreter','latex');  

bins=logspace(-2,3,40); [bincounts,whichbin] = histc(concens,bins);
for i_bin = 1:length(bins)-1, is_inbin = whichbin == i_bin; medianspot(i_bin) = median(meanSpot(is_inbin)); end
binsc = 0.5*(bins(1:end-1)+bins(2:end));
semilogx(binsc,medianspot,'LineWidth',2);
hold on
scatter(concens,meanSpot,3)
title("BNI1 spot counts in 2.5 $\mu m$ nuclear neighborhoods",'interpreter','latex')
xlabel("Number of RNA / Neighborhood volume",'interpreter','latex')
ylabel("Mean spot size",'interpreter','latex')

% figure(4)
% xtest = sort(concens);
% ytest = ones(length(xtest),1);
% slopeTest = (max(meanSpot)-1)/(xtest(end)-xtest(30));
% for i = 30:82
%     ytest(i) = 1+slopeTest*(xtest(i)-xtest(30));
% end
% plot(xtest,ytest)
% hold on
% 
% title("Synthetic Data")
% xlabel("Mean concentration")
% ylabel("Mean spot size")
% set(gca,'xscale','log')
% 
% meanSpotTest = ytest+(rand(length(xtest),1)-.5)/2;
% nthresh = 100;
% threshs = linspace(min(concens),max(concens),nthresh);
% slopes = (1:.01:3);
% f = zeros(length(meanSpotTest),1);
% fdiff = norm(f-meanSpotTest,1);
% slope_opt = 0;
% thresh_opt = 0;
% [x,xorder] = sort(concens);
% %y = meanSpotTest(xorder);
% 
% for thresh = threshs
%     for slope = slopes
%         ftemp = double(x < thresh);
%         thresh_ind = sum(ftemp);
%         
%         for i = (thresh_ind+1):length(x)
%             ftemp(i) = 1+slope*(x(i)-thresh);
%         end
%         %plot(x,ftemp)
%         %hold on
%         
%         diff = norm(ftemp-meanSpotTest,1);
%         if diff < fdiff
%            slope_opt = slope;
%            thresh_opt = thresh;
%            fdiff = diff;
%            f = ftemp;
%         end
%         
%             
%     end    
% end
% scatter(x,meanSpotTest)
% 
% plot(x,f)
% title(sprintf('Optimal threshold is %f mRNA', mean(dataVol(:,1))*thresh_opt))
% legend("Original Data","Perturbed Data", "Fit")