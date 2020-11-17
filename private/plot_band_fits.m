
function plot_band_fits(gelData, profileData)
%% plotting lanes and pocket fits
    imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on
    n = length(profileData.profiles);

    % plot pocket fits
    mu = profileData.aggregateFit(2);
    sig = profileData.aggregateFit(3) * profileData.sigma_integrate;
    for i=1:n
        x = profileData.lanePositions(i,1:2);
        plot(mean(x) , mu, 'r.');
        plot([mean(x) mean(x)], [mu-sig mu+sig], 'r');
    end

    % plot leading band fits
    for i=1:n
        mu = profileData.monomerFits(i,2);
        sig = profileData.monomerFits(i,3) * profileData.sigma_integrate;
        x = profileData.lanePositions(i,1:2);
        plot([x(1) x(2)] , [mu mu], 'r');
        plot([mean(x) mean(x)], [mu-sig mu+sig], 'r');
    end
end


