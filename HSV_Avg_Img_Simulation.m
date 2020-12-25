
%% Load Gaussian width in nm, based on experimental data

clear;close all;clc;

data_file_name = 'Grouped_Data_example.mat';

[all_accepted_radii_um_preEx, ...
    particle_avg_accepted_radii_um_preEx,...
    all_normalized_radii,...
    particle_avg_normalized_radii,...
    particle_std_normalized_radii ,...
    particle_std_radii ,...
    all_accepted_gaufit2_sigma_um_preEx,...
    particle_avg_accepted_gaufit2_sigma_um_preEx,...
    median_PAA_gaufit2_sigma_preEx_nm,...
    median_TG_gaufit2_sigma_preEx_nm,...
    median_pooled_gaufit2_sigma_preEx_nm,...
    sample_types] = load_grouped_data_and_compile_statistics(data_file_name);

%%

% Flip sample-type order of loaded variables
all_accepted_radii_um_preEx = fliplr(all_accepted_radii_um_preEx);
all_accepted_gaufit2_sigma_um_preEx = fliplr(all_accepted_gaufit2_sigma_um_preEx);
sample_types = fliplr(sample_types);
particle_avg_accepted_radii_um_preEx = fliplr(particle_avg_accepted_radii_um_preEx);
particle_std_radii = fliplr(particle_std_radii);

% ====== num of particle to include in the simulated average image ======
num_particles = [362 396];

% Set input value, for std of radii
override_circ_std_with_reported_value = 1;
if override_circ_std_with_reported_value
    circ_std = [9.2, 14.3]; % overriding with reported values
else
    circ_std = cellfun(@median,particle_std_radii).*1000; % nm
end

% % Model parameters
% mu_mode = 'randomly_selected_from_mean_particle_radii_of_gel_and_plus_or_minus_median_std_radii';
% sigma_mode = 'median_of_respective_gel';

profile_mag_normalization = 'peak_to_1';
simu_pixel_size = 0.01; % per pixel is 0.01 nm

%% Generate radial intensity profiles with Gaussians

xvals = 0:simu_pixel_size:400; % to 400 nm, to be consistent with actual sum images

gau_mu = cell(1,2);
gau_sigma = cell(1,2);
sum_int_profile = cell(2,8); % summed over particles
simu_std_of_radii = cell(1,2);

% Generate Gaussian parameters
for s = 1:length(sample_types)
    
    gau_mu{s} = zeros( num_particles(s) , 8 );
    gau_sigma{s} = zeros( num_particles(s) , 8 );
    
    simu_std_of_radii{s} = zeros( num_particles(s) , 1);
    
    for r = 1:8
        sum_int_profile{s,r} = zeros(size(xvals));
    end
    
    % Sample-wide constant parameters
    mu_median_of_pooled = median( [all_accepted_radii_um_preEx{1} ; all_accepted_radii_um_preEx{2}] ) .* 1000; %nm
    mu_median_of_respective_gel = median( [all_accepted_radii_um_preEx{s}] ) .* 1000; % nm
    mu_mean_of_respective_gel = mean( [all_accepted_radii_um_preEx{s}] ) .* 1000; % nm
    sigma_median_of_pooled = median( [all_accepted_gaufit2_sigma_um_preEx{1} ; all_accepted_gaufit2_sigma_um_preEx{2}] )  .* 1000; %nm
    sigma_median_of_respective_gel = median( [all_accepted_gaufit2_sigma_um_preEx{s}] ) .* 1000; % nm
    
    for p = 1 : num_particles(s)
        
        % mu_mode = randomly_selected_from_mean_particle_radii_of_gel_and_plus_or_minus_median_std_radii
        gau_mu{s}(p,:) = ones(1,8) .* particle_avg_accepted_radii_um_preEx{s}( ...
            ceil( rand(1,1).* length(particle_avg_accepted_radii_um_preEx{s}) ) ) ...
            .* 1000; % nm
        r_with_pos_sigma = datasample(1:8,4,'Replace',false);
        r_with_neg_sigma = 1:8;
        r_with_neg_sigma(r_with_pos_sigma) = [];
        gau_mu{s}(p,r_with_pos_sigma) = ...
            gau_mu{s}(p,r_with_pos_sigma)  ...
            + circ_std(s); % nm
        gau_mu{s}(p,r_with_neg_sigma) = ...
            gau_mu{s}(p,r_with_neg_sigma)  ...
            - circ_std(s); % nm
        
        simu_std_of_radii{s}(p) = std(gau_mu{s}(p,:),1);
        
        
        % sigma_mode = 'median_of_respective_gel'
        gau_sigma{s}(p,:) = sigma_median_of_respective_gel .* ones( 1 , 8 );
        
        for r = 1:8
            % Make Gaussian
            curr_int_profile = normpdf( xvals , ...
                gau_mu{s}(p,r) , gau_sigma{s}(p,r)  );
            
            switch profile_mag_normalization
                case 'peak_to_1'
                    curr_int_profile = ...
                        curr_int_profile .* (1 ./ max(curr_int_profile));
            end
            
            sum_int_profile{s,r} = sum_int_profile{s,r} + curr_int_profile;
            
        end
    end
end

% ======= Analysis of simulated results ========

sample_color = {[0 0 1],[0.91 0.41 0.17]};

BKG_intensity_for_simu_radii_profiles = 0;

% Extract radii and FWHM from simulated sum images
for s = 1:length(sample_types)
    for r = 1:8
        
        % Get radii from gaufit
        gaufit2_outputs(s,r) = perform_gaufit2_analysis( xvals , sum_int_profile{s,r} , BKG_intensity_for_simu_radii_profiles, 1, 0, 0);
        % Get FWHM
        profile_peak_indices(s,r) = find(sum_int_profile{s,r} == max(sum_int_profile{s,r}), 1);
        FWHM_outputs(s,r) = perform_FWHM_analysis_v2( sum_int_profile{s,r} ,...
            profile_peak_indices(s,r) , BKG_intensity_for_simu_radii_profiles, 0.5);
        
    end
end

% Final metrics from simulated sum images
for s = 1:length(sample_types)
    simu_particle_std_of_radii(s) = ...
        std( [gaufit2_outputs(s,:).mu] );
    simu_repre_FWHM_of_radii_profile(s) = ...
        median([FWHM_outputs(s,:).FWHM] .* simu_pixel_size);
end




%% Graphic display

% =======  Simulated radial line profiles  ========

figure(300);clf;

row_num = 2;
col_num = 2;
simu_radial_profile_subplot_index = [0 col_num] + 1;
simu_img_subplot_index = [0 col_num] + 2;

plot_other_gel = 0;

for s = 1:length(sample_types)
    
    subplot(row_num,col_num, simu_radial_profile_subplot_index(s) );
    hold on;
    
    
    title(sprintf([...
        'Simulated AVG-img intensity profiles' '\n' ...
        'Std of radii = ' num2str(simu_particle_std_of_radii(s),'%.1f') ' nm\n'...
        'Median FWHM = ' num2str(simu_repre_FWHM_of_radii_profile(s),'%.1f') ' nm'...
        ]));
    
    if s == 2
        xlabel('distance from center (nm)')
    end
    
    set(gca,'YTickLabel','')
    
    for r = 1:8
        
        plot(xvals,sum_int_profile{s,r},'Color',sample_color{s});
        
        if plot_other_gel
            s_other = 2 - (s~=1);
            h_other = plot(xvals,sum_int_profile{s_other,r},'Color',[sample_color{s_other},0.15]);
        end
        
        xlim([0 250]);
        
    end
end

% Generate 400nmx400nm simulated imaged
% Pixel size = 1 nm; same as actual image

% All 8 profiles are essentially identical, using index 1 as
% representative

repre_profile_index = 1;

simu_sum_image = cell(1,2);

for s = 1:length(sample_types)
    coordinate_range = -200:1:200;
    simu_sum_image{s} = zeros(length(coordinate_range));
    
    profile_to_use = sum_int_profile{s,repre_profile_index};
    
    for i = 1:length(coordinate_range)
        for j = 1:length(coordinate_range)
            if coordinate_range(i) == 0 & coordinate_range(j) == 0
                continue
            end
            dist_to_origin = sqrt(coordinate_range(i).^2+coordinate_range(j).^2);
            simu_sum_image{s}(i,j) = ...
                profile_to_use(ceil(dist_to_origin.*(1/simu_pixel_size)));
        end
    end
    
    figure(300);
    subplot(row_num,col_num,simu_img_subplot_index(s));
    
    axis_line_color = [1 1 1].*0.5;
    
    img_center_index = (size(simu_sum_image{s},1)+1)/2;
    img_full_length = size(simu_sum_image{s},1);
    
    imagesc(simu_sum_image{s}); colormap gray; axis square; axis xy; hold on;
    % !!!!!! using the same image to have same BC thres
    lin_img = simu_sum_image{s}(:);
    set(gca,'CLim',[prctile(lin_img,00) prctile(lin_img,100).*1.1]);
    plot([ img_center_index img_center_index], [1 img_full_length],'Color',axis_line_color)
    plot([1 img_full_length],[ img_center_index img_center_index], 'Color',axis_line_color)
    plot([ 1 img_full_length], [1 img_full_length],'Color',axis_line_color)
    plot([ 1 img_full_length], [img_full_length 1],'Color',axis_line_color)
    
    title(sprintf([ ...
        'Simulated AVG image\n' ...
        sample_types{s} '\n' '(n = ' num2str(num_particles(s)) ')' ]));
    
    set(gca,'XLim',[0 400],'YLim',[0 400]);
    
    adj_image_coordinates(200)
end


%% Helper functions

function adj_image_coordinates(offset)

% Adjust image coordinates
xtick_num_new = cellfun(@str2double,get(gca,'XTickLabel'))-offset;
for k = 1:length(xtick_num_new)
    xtick_str_new{k,1} = num2str( xtick_num_new(k) );
end
set(gca,'XTickLabel',xtick_str_new)

ytick_num_new = cellfun(@str2double,get(gca,'YTickLabel'))-offset;
for k = 1:length(ytick_num_new)
    ytick_str_new{k,1} = num2str(ytick_num_new(k) );
end
set(gca,'YTickLabel',ytick_str_new)
end

function [all_accepted_radii_um_preEx, ...
    particle_avg_accepted_radii_um_preEx,...
    all_normalized_radii,...
    particle_avg_normalized_radii,...
    particle_std_normalized_radii ,...
    particle_std_radii ,...
    all_accepted_gaufit2_sigma_um_preEx,...
    particle_avg_accepted_gaufit2_sigma_um_preEx,...
    median_PAA_gaufit2_sigma_preEx_nm,...
    median_TG_gaufit2_sigma_preEx_nm,...
    median_pooled_gaufit2_sigma_preEx_nm,...
    sample_types] = ...
    load_grouped_data_and_compile_statistics(data_file_name)

load(data_file_name);

% Load Grouped_Data from virion particle analysis, as the model depends on
% statistics generated from such experimental data

sample_types = {'TG','PAA'};

empty_container = {[],[]};

all_accepted_radii_um_preEx = empty_container;
particle_avg_accepted_radii_um_preEx = empty_container;
all_normalized_radii = empty_container;
particle_avg_normalized_radii = empty_container;
particle_std_normalized_radii = empty_container;
particle_std_radii = empty_container;
all_accepted_gaufit2_sigma_um_preEx = empty_container;
particle_avg_accepted_gaufit2_sigma_um_preEx = empty_container;

for p = 1:length(particle_data)
    
    sample_type_index = find(strcmp(sample_types, particle_data(p).sample_type));
    if length(sample_type_index) ~= 1
        continue
        error('Incorrect sample type label.');
    end
    
    all_accepted_radii_um_preEx{sample_type_index} = [...
        all_accepted_radii_um_preEx{sample_type_index} ; ...
        particle_data(p).accepted_radii_um_preEx_vector] ;
    
    particle_avg_accepted_radii_um_preEx{sample_type_index} = [...
        particle_avg_accepted_radii_um_preEx{sample_type_index} ; ...
        particle_data(p).avg_of_accepted_radii_um_preEx];
    
    all_normalized_radii{sample_type_index} = [...
        all_normalized_radii{sample_type_index} ; ...
        particle_data(p).normalized_radii_vector];
    
    particle_avg_normalized_radii{sample_type_index} = [...
        particle_avg_normalized_radii{sample_type_index} ; ...
        particle_data(p).avg_of_normalized_radii];
    
    particle_std_normalized_radii{sample_type_index} = [...
        particle_std_normalized_radii{sample_type_index} ; ...
        particle_data(p).std_of_normalized_radii];
    
    particle_std_radii{sample_type_index} = [...
        particle_std_radii{sample_type_index} ; ...
        particle_data(p).std_of_radii];
    
    all_accepted_gaufit2_sigma_um_preEx{sample_type_index} = [...
        all_accepted_gaufit2_sigma_um_preEx{sample_type_index} ; ...
        particle_data(p).accepted_single_profile_sigma_um_preEx_vector];
    
    particle_avg_accepted_gaufit2_sigma_um_preEx{sample_type_index} = [...
        particle_avg_accepted_gaufit2_sigma_um_preEx{sample_type_index} ; ...
        particle_data(p).avg_of_accepted_single_profile_sigma_um_preEx_vector];
    
end

all_accepted_gaufit2_sigma_um_preEx_pooled_both_gels = [all_accepted_gaufit2_sigma_um_preEx{1};all_accepted_gaufit2_sigma_um_preEx{2}];
median_PAA_gaufit2_sigma_preEx_nm = median(all_accepted_gaufit2_sigma_um_preEx{1}) .* 1000;
median_TG_gaufit2_sigma_preEx_nm = median(all_accepted_gaufit2_sigma_um_preEx{2}) .* 1000;
median_pooled_gaufit2_sigma_preEx_nm = median(all_accepted_gaufit2_sigma_um_preEx_pooled_both_gels) .* 1000;

end
