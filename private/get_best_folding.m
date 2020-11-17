function [data_out, cur_fig] = get_best_folding(profileData, gelInfo, gelData, show_summary_figure)
% @ step3
% compute various metrics from gel data
% determine best folding conditions

    %process data & normalize
    stapleNorm = double(profileData.stapleLine(1,2)) ./ double(profileData.stapleLine(:,2));
    %TODO find propper spreadNormfactor for plotting and quality metric balance
    spreadNormfactor = 15.0;
    %TODO extrapolating staple_line to ladder might be even better
    if profileData.has_ladder
        stapleNorm = [stapleNorm(1) stapleNorm(1) stapleNorm' stapleNorm(end) stapleNorm(end)];
        migrateNorm = profileData.monomerFits(1,2);
    else
        stapleNorm = [stapleNorm(1) stapleNorm' stapleNorm(end)];
        %TODO uses scaffold but no proportionality factor to correct for it
        %TDOD proportionality factor form database average of each scaffold 
        migrateNorm = profileData.monomerFits(1,2);
    end
    normfactor = stapleNorm' ./ migrateNorm';


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
    % TODO add ladder_migrate_error helps qunatify migrate error dispite normalisation
    % ladder_migrate_error = mono_migrate(1) - mono_migrate(end);  %mono_migrate(1) = 1 due to normalisation
     
    % calculate amount of monomer, smear and aggreagtes for best folding
    total_band = (profileData.monomerTotal+profileData.pocketTotal+profileData.smearTotal);
    fraction_monomer = profileData.monomerTotal ./ total_band;
    fraction_smear = profileData.smearTotal./ total_band;
    fraction_pocket = profileData.pocketTotal./ total_band;
    
    %TODO currently unbalanced as components are of different size (spread<<fraction)
    %TODO 1-mono_spread is only ok if monospread in [0,1], current normalisation factor does not guarantee that
    % compute quality metric based on monomer fraction and band width
    folding_quality_metric = fraction_monomer .* (1.0 - mono_spread);
    % compute quality metric based on monomer fraction and band width and migration distance
    if profileData.has_ladder
        mono_migrate_best = max(mono_migrate(3:end-2));
    else
        mono_migrate_best = max(mono_migrate(2:end-1));
    end
    folding_quality_metric_migrate = (mono_migrate./mono_migrate_best).^2 .* folding_quality_metric;


    %% find best folding
    function index = best_lane(metric, indices)
        if ~isempty(indices)
            [~, i_sort] = sort(metric(indices), 'descend');
            index = indices(i_sort(1));
        end
    end
    index_best = best_lane(folding_quality_metric, index_foldings);
    index_best_Tscrn = best_lane(folding_quality_metric, index_Tscrn);
    index_Mgscrn = index_Mgscrn(index_Mgscrn~=0); % remove zeros if people did not include all Mg samples
    index_best_Mgscrn = best_lane(folding_quality_metric, index_Mgscrn);
    index_RM = index_RM(index_RM~=0); % remove zeros if people did not include all RM samples
    index_best_RM = best_lane(folding_quality_metric, index_RM);

   
    % check if M20 is better than best 4h T
    if folding_quality_metric(index_M20) > folding_quality_metric(index_best_Tscrn)
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
        
        % AGE gel with band fits (best highlighted 
        plot_band_fits(gelData, profileData)
        plot(mean([profileData.lanePositions(index_best,1) profileData.lanePositions(index_best,2) ]) , ...
                profileData.monomerFits(index_best,2), 'go');
        plot(mean([profileData.lanePositions(index_best_Tscrn,1) profileData.lanePositions(index_best_Tscrn,2) ]) , ...
                profileData.monomerFits(index_best_Tscrn,2), 'g+');
        plot(mean([profileData.lanePositions(index_best_Mgscrn,1) profileData.lanePositions(index_best_Mgscrn,2) ]) , ...
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
        % TODO add ladder_migrate_error to migrate as error bar
        subplot(5,1,4)
        plot(mono_spread, '.-'), hold on
        rel_mono_migrate = mono_migrate./mono_migrate_best;
        plot(rel_mono_migrate, '.--'), hold on
        ylabel({'Normalized ', 'spread or rela. migr. distance'})
        set(gca, 'XTick', (1:n), 'XTickLabels', gelInfo.lanes, 'XLim', [1 n], 'YLim', [0 rel_mono_migrate(1)])
        legend({'Band spread', 'Migr. distance'}, 'location', 'best')
        grid on 
        
        % quality metric plot
        subplot(5,1,5)
        plot(folding_quality_metric, '.-'), hold on
        plot(folding_quality_metric_migrate, '.--'), hold on
        plot(index_best, folding_quality_metric(index_best), 'o'), hold on
        plot(index_best_Tscrn, folding_quality_metric(index_best_Tscrn), '+'), hold on
        plot(index_best_Mgscrn, folding_quality_metric(index_best_Mgscrn), 'x'), hold on
        ylabel({'Quality metric ', 'Monomer fraction/norm_width'})
        set(gca, 'XTick', (1:n), 'XTickLabels', gelInfo.lanes, 'XLim', [1 n])
        legend({'Quality metric', 'Quality metric migrate'}, 'location', 'best')
        grid on
    end
end
