
function plot_band_fits(gelData, profileData)
%% plotting lanes and pocket fits
    img_L = profileData.aggregateSelectedArea(1);
    img_R = img_L + profileData.aggregateSelectedArea(3);
    img_plot = gelData.images{1}(:,img_L : img_R);
    imagesc(img_plot, [0 3.*std(img_plot(:))]), axis image, colormap gray, hold on
    n = length(profileData.profiles);

    % plot pocket fits
    mu_p = profileData.aggregateFit(2);
    sig_p = profileData.aggregateFit(3) * profileData.sigma_integrate;
    limit = profileData.has_ladder +1;
    for i=1:n
        x = profileData.lanePositions(i,1:2) - double(img_L);
        mu = profileData.monomerFits(i,2);
        sig = profileData.monomerFits(i,3) * profileData.sigma_integrate;
        if i>limit && i <= n-limit
            c = [1.0, 0.0, 0.0];
            plot(mean(x), profileData.stapleLine(i-limit,2), '.', 'color', [1.0 0.77 0.13]);
        elseif mu < mu_p
            c = [1.0 1.0 1.0];
        else        
            c = [1.0 0.77 0.13];
        end
        plot(mean(x), mu_p, '.', 'color',  c);
        plot([mean(x) mean(x)], [mu_p-sig_p mu_p+sig_p], 'color',  c);
        % plot leading band fits
        plot([x(1) x(2)] , [mu mu], 'color',  c);
        plot([mean(x) mean(x)], [mu-sig mu+sig],  'color', c);
        
    end
end


