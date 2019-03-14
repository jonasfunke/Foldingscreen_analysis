%%
close all, clear all, clc
s = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/AM_Monomer_v1/AM_Monomer_v1_data.mat';
load(s)

profileData = analyze_gel_fractions(profileData, gelData);
profileData = compute_exp_parameter(profileData, 1.5);
get_best_folding(profileData, gelInfo)



%%
close all, clear all, clc
s = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/CS_T1v1/CS_T1v1_data.mat';
load(s)

profileData = analyze_gel_fractions(profileData, gelData);
profileData = compute_exp_parameter(profileData, 1.5);
get_best_folding(profileData, gelInfo)


%%

close all, clear all, clc
[fname pname] = uigetfile('Select mat file', '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/*.mat')

load([pname fname])

profileData = analyze_gel_fractions(profileData, gelData);

profileData = compute_exp_parameter(profileData, 1.5);
get_best_folding(profileData, gelInfo)










%%
%s = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/CS_T1v1/CS_T1v1_data.mat';

%t = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/CS_T1v1/170117_T3DI_folding_screen_300V-[EtBr].tif';
%%


%%
close all
i = 10;

plot( fliplr(profileData.profiles{i}')), hold on

plot(cumsum(fliplr(profileData.profiles{i}')))


%%


close all

for i=1:length(profileData.profiles)

    %min_peak_height = median(profileData.profiles{i});
    n_peaks = 0;
    p_start = 0.000001;
    while n_peaks~=1
     
        min_peak_height = 0.005*sum(profileData.profiles{i}(100:end));
        findpeaks(profileData.profiles{i}(100:end), 'MinPeakHeight', min_peak_height, 'MinPeakWidth', 1)
      
        
    end
    pause
end
%%

%%

plot_image_ui(gelData.images{1});
h = imrect;

wait(h);
selectedArea = int32(getPosition(h));


close all

%%
fits = cell(length(profileData.profiles),1);
tmp = zeros(2, length(profileData.profiles));
for i=1:length(profileData.profiles)
    %subplot(length(profileData.profiles), 1, i)
    figure
    
    tmp(i,:) = max(profileData.fullProfiles{i}(selectedArea(2):selectedArea(2)+selectedArea(4)));
    
    y = profileData.fullProfiles{i}(selectedArea(2):selectedArea(2)+selectedArea(4));
    x = double(selectedArea(2):selectedArea(2)+selectedArea(4));
    
    fits{i} = fit(x', y, 'gauss1');
    plot(x, profileData.fullProfiles{i}(selectedArea(2):selectedArea(2)+selectedArea(4))), hold on

    plot(x, fit_val(x))
    
    
end
%%

subplot(4, 1, 1)
tmp = zeros(length(profileData.profiles),1);
for i=1:length(profileData.profiles)
    tmp(i) = fits{i}.a1;
end
plot(tmp)

subplot(4, 1, 2)
tmp = zeros(length(profileData.profiles),1);
for i=1:length(profileData.profiles)
    tmp(i) = fits{i}.c1;
end
plot(tmp)

subplot(4, 1, 3)
tmp = zeros(length(profileData.profiles),1);
for i=1:length(profileData.profiles)
    tmp(i) = fits{i}.b1;
end
plot(tmp)
%%
subplot(4, 1, 4)
tmp = zeros(length(profileData.profiles),1);
for i=1:length(profileData.profiles)
    tmp(i) = fits{i}.a1/fits{i}.c1;
end
plot(tmp)
%%





