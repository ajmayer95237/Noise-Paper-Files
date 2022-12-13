%%
tic
load('../data/raw/spotinfo_all.mat'); % load file with all spot coordinates and intensity values
N_allhyph = size(spotinfo_all,1); % count total number of hyphae analyzed across all conditions
lastfilename = ''; % placeholder that helps us to decide whether we have changed image in going from one hypha to the next one
maskdir = '../data/3Dmasks';

maxN_samples = 10; % max no. spheres analyzed per hypha

maxN_spheres = 1e5; % ball park for max number of spheres across all hyphae
% maxN_spheres = maxN_samples*N_allhyph; % better bound

N_all = zeros(maxN_spheres,1); % number of all mRNAs across all spheres
N_incluster = zeros(maxN_spheres,1); % number of mRNAs that were detected within a cluster
CVsquared = zeros(maxN_spheres,1);
entropyballs = zeros(maxN_spheres,1);
meanconcens = zeros(maxN_spheres,1);
meanSpotWeights = zeros(maxN_spheres,1);
which_condition_sphere = zeros(maxN_spheres,6); % a set of variables telling us what kind of hypha the sphere was located in
spherevol = zeros(maxN_spheres,1); % volume of the sphere
analyzeddir = '../data/raw';

is_implot = false; % do we make plots
is_0counts = true; % do we show mRNAs assigned weight 0 in the plots?

clustersize = 0.5; % separation (in um) of spots that are considered to be in a single cluster together
sphererad = 2.5; % radius of the sphere we are measuring
subsphererad = 1.0; %radius of sphere within nuclear neighborhood
xsphere = -ceil(sphererad/.11):ceil(sphererad/.11); % set up the sphere we will do measurements in [3D mask]
zsphere = -ceil(sphererad/.3):ceil(sphererad/.3); 
xsubsphere = -ceil(subsphererad/.11):ceil(subsphererad/.11); % set up the sphere we will do measurements in [3D mask]
zsubsphere = -ceil(subsphererad/.3):ceil(subsphererad/.3); 
[xsphere,ysphere,zsphere] = meshgrid(xsphere,xsphere,zsphere);
[xsubsphere,ysubsphere,zsubsphere] = meshgrid(xsubsphere,xsubsphere,zsubsphere);
is_sphere = sqrt(0.11^2*xsphere.^2+.11^2*ysphere.^2+0.3^2*zsphere.^2)<sphererad;
is_subsphere = sqrt(0.11^2*xsubsphere.^2+.11^2*ysubsphere.^2+0.3^2*zsubsphere.^2)<subsphererad;


i_sampleall = 1; % which sample we are on
is_otsuweight = false; % are we thresholding using otsu or using quantile scores. Otsu's method looked better for whi3
is_quantileweight = true; % by eye, quantile weights looked better for cln3 and bni1

maxN_nucl = 10000; % store data about nuclei
%maxSubsamples = 5;
N_mRNAinnucl = zeros(maxN_nucl,1); % how many mRNAs in each nucleus
which_condition_nucl = zeros(maxN_nucl,6); % what kind of hypha was the nucleus drawn
i_allnucl = 1;
%%
for i_probe = 1:3% loop over probes: cln3, bni1, whi 3
    
    for i_strain = [822 882 909] % loop over strains
        
        for i_pk = 1:2 % loop over whether PK or not PK
            switch i_pk
                case 1
                    is_pk = false;
                    pkstr = '';
                case 2
                    is_pk = true;
                    pkstr = ' pk';
            end
            
            switch i_probe
                case 1
                    probename = 'cln3';
                    
                case 2
                    probename = 'bni1';
                    
                case 3
                    probename = 'whi3';
                    
            end
            
            
            spotinfofile = [analyzeddir filesep 'spotinfo_all.mat'];
            if is_otsuweight
                concenfielddir = [analyzeddir filesep 'otsu'];
            elseif is_quantileweight
                concenfielddir = [analyzeddir filesep 'quantile']; % create the directory for storing the data
            else
                error('Either is_otsuweight or is_quantileweight must be set equal to true');
            end
            
            if ~exist(concenfielddir,'dir'), mkdir(concenfielddir); end
            NumberSteps = 50;

            concenfielddatafile = [concenfielddir filesep 'meanConcenmeanSpots.mat'];                                            
                        
            lastfilename = '';            
            allhyph = find(is_goodhyph_all);
            for i = 1:numel(allhyph) % only analyze the hyphae that have good intensities
                i
                i_allhyph = allhyph(i);
                imfileprefix = spotinfo_all{i_allhyph,1}; % find the image that this hypha is contained in
                %                     figfile = [figdir filesep imfileprefix '-RNAcheck.fig']
                if spotinfo_all{i_allhyph,2} == i_strain ...
                        && spotinfo_all{i_allhyph,5} == is_pk ...
                        && (strcmpi(spotinfo_all{i_allhyph,3},probename) ...
                        || strcmpi(spotinfo_all{i_allhyph,4},probename))
                    
                    
                    whichhyph = spotinfo_all{i_allhyph,6};
                    if strcmpi(spotinfo_all{i_allhyph,3},'cln3'), i_probe1=1; end
                    if strcmpi(spotinfo_all{i_allhyph,3},'bni1'), i_probe1=2; end
                    
                    i_probe2 = 3;
                    
                    if ~strcmpi(imfileprefix,lastfilename) % load the files containing the hyphal and nuclear masks
                        maskfile = [maskdir filesep imfileprefix '-mask.mat'];
                        lastfilename = imfileprefix;
                        load(maskfile,'hyph3D_label_uncompress');
                        
                        nuclfile = [analyzeddir filesep imfileprefix '-nucl.mat'];
                        load(nuclfile,'nucl3D_label_uncompress');
                    end
                    
                    is_hyph3D = hyph3D_label_uncompress==whichhyph; % find the specific hypha of interest
                    label_nucl3D = nucl3D_label_uncompress.*is_hyph3D;  % labeled nuclei in the hypha of interest                                      
                    
                    [Ny,Nx,Nz] = size(hyph3D_label_uncompress);
                    
                    if strcmpi(spotinfo_all{i_allhyph,3},probename) % read out the data about the locations and intensities of spots
                        spot_intintens = spotinfo_all{i_allhyph,8};
                        spot_snr = spotinfo_all{i_allhyph,10};
                        spot_xyz = spotinfo_all{i_allhyph,11};
                    elseif strcmpi(spotinfo_all{i_allhyph,4},probename)
                        probe1name = spotinfo_all{i_allhyph,3};
                        spot_intintens = spotinfo_all{i_allhyph,12};
                        spot_snr = spotinfo_all{i_allhyph,14};
                        spot_xyz = spotinfo_all{i_allhyph,15};
                    else
                        spot_snr = [];
                        spot_intintens = [];
                    end
                    
                    is_spotweightflag = false; 
                    
                    if is_otsuweight % we are using otsu's method to assign weights to spots
                        
                        if ~isempty(spot_intintens) % use the integrated intensity
                                                        
                            [intensthresh1,metric1] = multithresh(spot_intintens,1); % threshold is between not-spots and spots
                            [intensthresh2,metric2] = multithresh(spot_intintens(spot_intintens>intensthresh1),1); % gives the threshold between single mRNAs and multiple mRNAs
                            
                            if metric1>0.1 && metric2>0.1
                                
                                dthresh = intensthresh2 - intensthresh1;
                                
                                spotweight = floor((spot_intintens-intensthresh1)/dthresh); % assigns each spot a weight, showing our estimate of how many mRNAs it contains
                                %                         if is_0counts
                                %                             spotweight = spotweight+1;
                                %                         end
                                spotweight = max(spotweight+1,0);  
                                
                                hyphvol = spotinfo_all{i_allhyph,7};
                                
                                is_0weight = spotweight==0; % indicate which spots contain no mRNAs
                                
                                spotweight(is_0weight) = []; % remove spots that don't contain any mRNAs
                                spot_xyz(is_0weight,:) = [];                                                                
                                
                                is_spotweightflag = true;
                            end
                        end
                    elseif is_quantileweight
                        is_spot = spot_snr>1.4; % identify every spot with signal-to-noise>1.4.
                        min_intens = median(spot_intintens(spot_intintens>quantile(spot_intintens(is_spot),0.1) ...
                            & spot_intintens<quantile(spot_intintens(is_spot),0.25))); % integrated intensity should be between 10% and 25%-ile. Median
                        % of all spots in this range is the minimum
                        % intensity.
                        
                        spotweight = round(spot_intintens./min_intens);
                        is_0weight = spotweight<=0;
                        spotweight(is_0weight) = [];
                        spot_xyz(is_0weight,:) = [];
                        is_spotweightflag = true;
                        
                        
                    end
                    %%
                    if is_spotweightflag
                        
                        spot_inds = sub2ind([Ny,Nx,Nz],round(spot_xyz(:,2)),round(spot_xyz(:,1)),...
                                    round(spot_xyz(:,3))); % find linear indices for all spots
                        whichnucl = nucl3D_label_uncompress(spot_inds); % which spots are in nuclei

                        N_nucl = max(label_nucl3D(:));
                        nuc_centroids = zeros(3,N_nucl);
                        for i_nucl = 1:N_nucl
                            nuci = label_nucl3D == i_nucl;
                            if max(nuci(:)) == 0
                               nuc_centroids(:,i_nucl) = [0; 0; 0];
                               i_allnucl = i_allnucl+1;
                            else
                                nuci_c = regionprops3(nuci,'Centroid').Centroid;
                                if length(nuci_c(:)) > 3
                                   nuc_centroids(:,i_nucl) = [0; 0; 0];
                                   i_allnucl = i_allnucl+1;
                                else
                                    nuc_centroids(:,i_nucl) =[nuci_c(1); nuci_c(2); nuci_c(3)]; %nuci_c in form [cx,cy,cz]
                                    N_mRNAinnucl(i_allnucl) = sum(spotweight(whichnucl==N_nucl));
                                    which_condition_nucl(i_allnucl,:) = [i_probe1 i_probe2 i_probe i_strain i_pk i_allhyph];
                                    i_allnucl = i_allnucl+1;
                                end
                            end
                        end
                        
                        spotweight(whichnucl>0) = []; % remove any spot that is contained in a nucleus
                        spot_xyz(whichnucl>0,:) = [];
                      
                        
                        is_nucl = is_hyph3D & label_nucl3D > 0;
                        is_hyph3D = is_hyph3D & label_nucl3D==0; % create a modified mask of the hyphae, deleting nuclei
                        is_hyph3D_tosample = is_hyph3D; % reference to the mask - remove voxels as they get sampled
                        i_sample = 1;      
                        
                        %%                                                           
                        while i_sample<=N_nucl && any(is_hyph3D_tosample(:)) 
                              sphere_c = round(nuc_centroids(:,i_sample)); %[cx,cy,cz]
                              if sphere_c(1) == 0 %Nucleus not contained in hypha
                                 i_sample = i_sample+1;
                                 continue
                              end
                              
                              sample_mask = false(Ny,Nx,Nz);
                              sample_mask(sphere_c(2),sphere_c(1),sphere_c(3)) = 1;
                              sample_mask = imdilate(sample_mask,is_sphere) & is_hyph3D_tosample; % mask of the spherical volume I am sampling
                              maskinds = find(sample_mask);
                              spot_insphere = (0.11^2*(spot_xyz(:,2)-sphere_c(2)).^2 ...
                                         +0.11^2*(spot_xyz(:,1)-sphere_c(1)).^2 ...
                                         +0.3^2*(spot_xyz(:,3)-sphere_c(3)).^2) < sphererad^2;
                              spherevol(i_sampleall) = sum(sample_mask(:))*0.11^2*0.3; % volume of the sphere
                              N_all(i_sampleall) = sum(spotweight(spot_insphere)); % number of mRNAs in sphere                           
                              if ~isempty(spotweight(spot_insphere))
                                meanSpotWeights(i_sampleall) = mean(spotweight(spot_insphere));    
                              end
                              
                              is_hyph3D_tosample = is_hyph3D_tosample & (1-sample_mask);
                              which_condition_sphere(i_sampleall,:) = [i_probe1 i_probe2 i_probe i_strain i_pk i_allhyph];
                              
                              i_sample = i_sample+1;
                              i_sampleall = i_sampleall+1;
                              
                        end
                    
                end
            end    
            
        end
        
        
    end
    
end
end


N_all(i_sampleall:end) = [];
meanSpotWeights(i_sampleall:end) = [];
which_condition_sphere(i_sampleall:end,:) = [];
spherevol(i_sampleall:end) = [];

save(concenfielddatafile,'N_all','which_condition_sphere',"meanSpotWeights");
