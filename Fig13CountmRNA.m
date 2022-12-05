M = load('../data/raw/quantile/BurstSize.mat');
keyboard
mRNA = M.N_mRNAinnucl;

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
    nucleus_i = M.which_condition_nucl(i,:);
    probe_i = nucleus_i(3);
    strain_i = nucleus_i(4);
    pk_i = nucleus_i(5);
    hypha_i = nucleus_i(6);
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
dataN_mRNA = M.N_mRNAinnucl(data);
%dataEntropy = M.entropyballs(data);
%dataconcen = M.meanconcens(data);
minN_mRNA = min(dataN_mRNA);
maxN_mRNA = max(dataN_mRNA);


figure(3)
notNAind = ~isnan(dataN_mRNA(:,1));
dataN_mRNA=dataN_mRNA(notNAind,:);

nonzeroBurstsinds = find(dataN_mRNA(:,1) > 0);
nonzeroBursts = dataN_mRNA(nonzeroBurstsinds,1);
histogram(nonzeroBursts)

