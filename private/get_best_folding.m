function [data_out, cur_fig] = get_best_folding(profileData, gelInfo, gelData, show_summary_figure)
% @ step3
% compute various metrics from gel data
% determine best folding conditions

    %process data & normalize
    stapleNorm = double(profileData.stapleLine(1,2)) ./ double(profileData.stapleLine(:,2));
    %TODO find propper spreadNormfactor for plotting and quality metric balance
    spreadNormfactor = 15.0;
    
    % extrapolating staple norm
    offL =  stapleNorm(2) - stapleNorm(1);
    nL = stapleNorm(1) - offL;
    offR =  stapleNorm(end-1) - stapleNorm(end);
    nR = stapleNorm(end) - offR;
    if profileData.has_ladder
        nL2 = stapleNorm(1) - 2* offL;
        nR2 = stapleNorm(end) - 2* offR;
        stapleNorm = [nL2 nL stapleNorm' nR nR2];
        migrateNorm = profileData.monomerFits(1,2);
    else
        stapleNorm = [nL stapleNorm' nR];
        %TODO uses scaffold but no proportionality factor to correct for it
        %TODO proportionality factor form database average of each scaffold 
        migrateNorm = profileData.monomerFits(1,2);
    end
    %NOTE: sqrt of stapleNorm to reduce its effect on the normalisation.
    %Full effect seems to overcorrect migration distances (theoretical
    %explanation missing!)
    normfactor =  sqrt(stapleNorm') ./ migrateNorm';
    

    %% get lane indices
    function indices = get_lane_index(lanes, identifier)
        indices = [];
        for l=1:length(lanes)
            name = upper(strtrim(lanes{l}));
            comp = strcmpi(name(1:length(identifier)), identifier);
            if comp
                indices = [indices l];
            end
        end
    end
    index_Tscrn = get_lane_index(gelInfo.lanes, 'T');
    index_Mgscrn = get_lane_index(gelInfo.lanes, 'M');
    index_RM = get_lane_index(gelInfo.lanes, 'RM');
    index_foldings = sort([index_Tscrn index_Mgscrn index_RM]);
    % TODO better solution with find?
    % index M20
    for i=1:length(gelInfo.lanes)
        cur_name = strtrim(gelInfo.lanes{i});
        if strcmpi(cur_name, 'M20')
            index_M20 = i;
        end
    end 


    %% metrics

    mono_spread = profileData.monomerFits(:,3) .* normfactor .* spreadNormfactor;
    mono_migrate = profileData.monomerFits(:,2) .* normfactor;
     
    % calculate amount of monomer, smear and aggreagtes for best folding
    total_band = (profileData.monomerTotal+profileData.pocketTotal+profileData.smearTotal);
    fraction_monomer = profileData.monomerTotal ./ total_band;
    fraction_smear = profileData.smearTotal./ total_band;
    fraction_pocket = profileData.pocketTotal./ total_band;
    
    % compute quality metric based on monomer fraction and band width
    %NOTE currently unbalanced as components are of different size (spread<<fraction)
    %NOTE 1-mono_spread is only ok if monospread in [0,1], current normalisation factor does not guarantee that
    folding_quality_metric = fraction_monomer .* (1.0 - mono_spread);
    if profileData.has_ladder
        folding_quality_metric(1:2) = 0.0;
        folding_quality_metric(end-1:end) = 0.0;
    else
        folding_quality_metric(1) = 0.0;
        folding_quality_metric(end) = 0.0;
    end
    
    % compute relative migration distance
    if profileData.has_ladder
        mono_migrate_best = max(mono_migrate(3:end-2));
        mono_spread_best = max(mono_spread(3:end-2));
    else
        mono_migrate_best = max(mono_migrate(2:end-1));
        mono_spread_best = max(mono_spread(3:end-2));
    end
    rel_mono_migrate = mono_migrate./mono_migrate_best;
    rel_mono_spread = mono_spread./mono_spread_best;
    
    % compute migration_distance corrected quality metric
    ladder_migrate_error = abs(rel_mono_migrate(1) - rel_mono_migrate(end));
    ladder_spread_error = abs(rel_mono_spread(1) - rel_mono_spread(end));

    n = length(profileData.profiles);
    folding_quality_metric_migrate = zeros(n, 1);
    cutoff = 0.75; 
    tolerance = 1.0 * ladder_migrate_error;
    
    for i = 1:n
        if rel_mono_migrate(i) < cutoff
            folding_quality_metric_migrate(i) = 0.0;
        elseif (1-rel_mono_migrate(i)) < tolerance
            folding_quality_metric_migrate(i) = folding_quality_metric(i);
        else
            %NOTE: ^2 increasses the harsheness of the migration penalty to
            %balance the migration quality towards quality (no explanation
            %theoretical explanation given)
            migrate_penalty = (rel_mono_migrate(i) - cutoff) / (1.0 - cutoff);
            folding_quality_metric_migrate(i) = migrate_penalty^2 .* folding_quality_metric(i); 
        end
    end
            

    %% find best folding using migrate quality metric
    function index = best_lane(metric, indices)
        if ~isempty(indices)
            [~, i_sort] = sort(metric(indices), 'descend');
            index = indices(i_sort(1));
        end
    end
    index_best = best_lane(folding_quality_metric_migrate, index_foldings);
    index_best_Tscrn = best_lane(folding_quality_metric_migrate, index_Tscrn);
    index_Mgscrn = index_Mgscrn(index_Mgscrn~=0); % remove zeros if people did not include all Mg samples
    index_best_Mgscrn = best_lane(folding_quality_metric_migrate, index_Mgscrn);
    index_RM = index_RM(index_RM~=0); % remove zeros if people did not include all RM samples
    index_best_RM = best_lane(folding_quality_metric_migrate, index_RM);

   
    % check if M20 is better than best 4h T
    if folding_quality_metric_migrate(index_M20) > folding_quality_metric_migrate(index_best_Tscrn)
     	M20_better_than_bestT = true;
    else
        M20_better_than_bestT = false;
    end


   
    %% save 
    data_out.fractionMonomer = fraction_monomer;
    data_out.fractionSmear = fraction_smear;
    data_out.fractionPocket = fraction_pocket;
    data_out.bestFolding = gelInfo.lanes{index_best};
    data_out.bestFoldingIndex = index_best;
    data_out.bestTscrn = gelInfo.lanes{index_best_Tscrn}; 
    data_out.bestTscrnIndex = index_best_Tscrn; 
    data_out.bestMgscrn = gelInfo.lanes{index_best_Mgscrn};
    data_out.bestMgscrnIndex = index_best_Mgscrn;
    data_out.bestRM = gelInfo.lanes{index_best_RM};
    data_out.bestRMIndex = index_best_RM;

    data_out.M20BetterThanTscrn = M20_better_than_bestT;
    data_out.bandWidthNormalized = mono_spread;
    data_out.migrationDistanceNormalized = mono_migrate;
    data_out.qualityMetric = folding_quality_metric;

    
    %% print results
    %TODO: move to seperate function   
   
    % calculate amount of monomer, smear and aggreagtes for best folding condition and M20
    function string_percentage = p2str(value)
        string_percentage = num2str(round(100*value));
    end
    disp(['Best folding: Lane ' num2str(index_best) ' (' gelInfo.lanes{index_best} ')'])
    disp(['Best Mg-Screen: Lane ' num2str(index_best_Mgscrn) ' (' gelInfo.lanes{index_best_Mgscrn} ')'])
    disp(['Best T-Screen: Lane ' num2str(index_best_Tscrn) ' (' gelInfo.lanes{index_best_Tscrn} ')'])
    disp(['Best RM: Lane ' num2str(index_best_RM) ' (' gelInfo.lanes{index_best_RM} ')'])
    disp(['Is M20 better than best 4h folding (' gelInfo.lanes{index_best_Tscrn} '): ' num2str(M20_better_than_bestT)])

    disp(['Fraction of monomers/smear/pocket for the best foldings: ' p2str(fraction_monomer(index_best)) '%, ' ...
       p2str(fraction_smear(index_best)) '%, ', p2str(fraction_pocket(index_best)) '%'])
   
    %% figure
    if show_summary_figure
        cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 20 30], 'PaperSize', [20 30]);
        subplot(5,1,1:2)
        % AGE gel with band fits (best highlighted)
        plot_band_fits(gelData, profileData)
        shift = profileData.aggregateSelectedArea(1);
        plot(mean([profileData.lanePositions(index_best,1) profileData.lanePositions(index_best,2) ]) - shift, ...
                profileData.monomerFits(index_best,2), 'go');
        plot(mean([profileData.lanePositions(index_best_Tscrn,1) profileData.lanePositions(index_best_Tscrn,2) ]) - shift , ...
                profileData.monomerFits(index_best_Tscrn,2), 'g+');
        plot(mean([profileData.lanePositions(index_best_Mgscrn,1) profileData.lanePositions(index_best_Mgscrn,2) ]) - shift , ...
                profileData.monomerFits(index_best_Mgscrn, 2), 'gx');
        title(['Band positions with sigma=' num2str(profileData.sigma_integrate)]);
        
        n = length(profileData.profiles);
        % fraction per lane plot
        subplot(5,1,3)
        plot(fraction_monomer, '.-'), hold on
        plot(fraction_smear, '.-'), hold on
        plot(fraction_pocket, '.-'), hold on
        ylabel('Fraction')
        set(gca, 'XTick', (1:n), 'XTickLabels', gelInfo.lanes, 'XLim', [1 n])
        legend({'monomer', 'smear', 'pocket'}, 'location', 'best')
        grid on
        
        % migration distance & spread plot
        err_migrate = ladder_migrate_error * ones(length(rel_mono_migrate),1);
        err_spread = ladder_spread_error * ones(length(rel_mono_spread),1);
        subplot(5,1,4)
        errorbar(rel_mono_spread, err_migrate, '.-'), hold on
        errorbar(rel_mono_migrate, err_spread, '.--'), hold on
        ylabel({'Normalized ', 'rel_spread or rela. migr. distance'})
        set(gca, 'XTick', (1:n), 'XTickLabels', gelInfo.lanes, 'XLim', [1 n], 'YLim', [0 rel_mono_migrate(1)])
        legend({'Band spread', 'Migr. distance'}, 'location', 'best')
        grid on 
        
        % quality metric plot
        subplot(5,1,5)
        plot(folding_quality_metric, '.-'), hold on
        plot(folding_quality_metric_migrate, '.--'), hold on
        plot(index_best, folding_quality_metric_migrate(index_best), 'o'), hold on
        plot(index_best_Tscrn, folding_quality_metric_migrate(index_best_Tscrn), '+'), hold on
        plot(index_best_Mgscrn, folding_quality_metric_migrate(index_best_Mgscrn), 'x'), hold on
        ylabel({'Quality metric ', 'Monomer fraction/norm_width'})
        set(gca, 'XTick', (1:n), 'XTickLabels', gelInfo.lanes, 'XLim', [1 n])
        legend({'Quality metric', 'Quality metric migrate'}, 'location', 'best')
        grid on
    end
end
