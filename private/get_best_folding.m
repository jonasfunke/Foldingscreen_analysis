function [data_out, cur_fig] = get_best_folding(profileData, gelInfo, gelData, show_summary_figure)
% @ step3
% compute various metrics from gel data
% determine best folding conditions

    %process data & normalize
    stapleNorm = double(profileData.stapleLine(:,2)) ./ double(profileData.stapleLine(1,2));
    spreadNormfactor = 3.0;
    if profileData.has_ladder
        stapleNorm = [stapleNorm(1) stapleNorm(1) stapleNorm' stapleNorm(end) stapleNorm(end)];
        migrateNorm = profileData.monomerFits(1,2);
    else
        stapleNorm = [stapleNorm(1) stapleNorm' stapleNorm(end)];
        %TODO uses scaffold but no proportionality factor to correct for it
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
    %index_scaffold = get_lane_index(gelInfo.lanes, "SC", true);

    index_foldings = sort([index_Tscrn index_Mgscrn index_RM]);
    % TODO find?
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
    fraction_monomer = profileData.monomerTotal./ total_band;
    fraction_smear = profileData.smearTotal./ total_band;
    fraction_pocket = profileData.pocketTotal./ total_band;
    
    % compute quality metric based on monomer fraction and band width
    folding_quality_metric = profileData.monomerTotal ./ (profileData.monomerTotal+profileData.pocketTotal+profileData.smearTotal) ./ mono_spread;
    % compute quality metric based on monomer fraction and band width and migration distance
    folding_quality_metric_2 = mono_migrate .* folding_quality_metric;


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

    
    %% figure
       
   
    % calculate amount of monomer, smear and aggreagtes for best folding
    % condition and M20
    disp(['Best folding: Lane ' num2str(index_best) ' (' gelInfo.lanes{index_best} ')'])
    disp(['Best Mg-Screen: Lane ' num2str(index_best_Mgscrn) ' (' gelInfo.lanes{index_best_Mgscrn} ')'])
    disp(['Best T-Screen: Lane ' num2str(index_best_Tscrn) ' (' gelInfo.lanes{index_best_Tscrn} ')'])

    disp(['Best RM: Lane ' num2str(index_best_RM) ' (' gelInfo.lanes{index_best_RM} ')'])
    disp(['Is M20 better than best 4h folding (' gelInfo.lanes{index_best_Tscrn} '): ' num2str(M20_better_than_bestT)])

    disp(['Fraction of monomers/smear/pocket for the best foldings: ' num2str(round(100*fraction_monomer(index_best))) '%, ' ...
       num2str(round(100*fraction_smear(index_best))) '%, ', ...
       num2str(round(100*fraction_pocket(index_best))) '%'])

    disp(['Fraction of monomers/smear/pocket for M20: ' num2str(round(100*fraction_monomer(index_M20))) '%, ' ...
       num2str(round(100*fraction_smear(index_M20))) '%, ', ...
       num2str(round(100*fraction_pocket(index_M20))) '%'])
 
   
    if show_summary_figure
  
       cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 20 30], 'PaperSize', [20 30]);

        %subplot(4, 1, 1)
        %plot(height_width, '.-')
        %xlabel('Lane')
        %ylabel('Height/width')
        %set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)] )

        %subplot(4,1,2)
        %plot(profileData.monomerTotal, '.-')
        %xlabel('Lane')
        %ylabel('monomer absolute')
        %set(gca, 'XTick', [1:length(profileData.profiles) ] )

    %     subplot(4,1,1)
    %     plot(width, '.-')
    %     xlabel('Lane')
    %     ylabel('Width')
    %     set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)])
    % 
    %     subplot(4,1,2)
    %     plot(widths, '.-')
    %     xlabel('Lane')
    %     ylabel('Width')
    %     set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)])

       subplot(5,1,1:2)
       imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on

       % plot pocket fits
        for i=1:length(profileData.profiles)
            plot(mean(profileData.lanePositions(i,1:2)), profileData.aggregateFit(2), 'r.')
            plot(mean(profileData.lanePositions(i,1:2))*[1 1], ...
                [profileData.aggregateFit(2)-profileData.sigma_integrate*profileData.aggregateFit(3) ...
                profileData.aggregateFit(2)+profileData.sigma_integrate*profileData.aggregateFit(3)], 'r')
        end

        % plot leading band fits
        for i=1:length(profileData.profiles)
            plot([profileData.lanePositions(i,1) profileData.lanePositions(i,2) ] , ...
                [profileData.monomerFits(i,2) profileData.monomerFits(i,2)], 'r')
            plot(mean(profileData.lanePositions(i,1:2))*[1 1], ...
                [profileData.monomerFits(i,2)-profileData.sigma_integrate*profileData.monomerFits(i,3) ...
                profileData.monomerFits(i,2)+profileData.sigma_integrate*profileData.monomerFits(i,3)], 'r')
        end


        plot(mean([profileData.lanePositions(index_best,1) profileData.lanePositions(index_best,2) ]) , ...
                profileData.monomerFits(index_best,2), 'go')
        plot(mean([profileData.lanePositions(index_best_Tscrn,1) profileData.lanePositions(index_best_Tscrn,2) ]) , ...
                profileData.monomerFits(index_best_Tscrn,2), 'g+')
        plot(mean([profileData.lanePositions(index_best_Mgscrn,1) profileData.lanePositions(index_best_Mgscrn,2) ]) , ...
                profileData.monomerFits(index_best_Mgscrn, 2), 'gx')
        title(['Band positions with sigma=' num2str(profileData.sigma_integrate)])
        
%         subplot(4,1,3)
%         plot(fraction_monomer, '.-'), hold on
%         plot(fraction_smear, '.-'), hold on
%         plot(fraction_pocket, '.-'), hold on
%         xlabel('Lane')
%         ylabel('monomer/(monomer+pocket+smear)')
%         set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)])
%         legend({'monomer', 'smear', 'pocket'})
        %subplot(4,1,4)
        %mig_distance = zeros(length(profileData.profiles),1);
        %for i=1:length(profileData.profiles)
        %    mig_distance(i) = profileData.monomerFits{i}.b1-profileData.aggregateFit.b1;
        %end
        %plot(mig_distance, '.-')
        %xlabel('Lane')
        %ylabel('Migration distance (pocket to monomer) [px]')
        %set(gca, 'XTick', [1:length(profileData.profiles) ] )

        subplot(5,1,3)
        plot(fraction_monomer, '.-'), hold on
        plot(fraction_smear, '.-'), hold on
        plot(fraction_pocket, '.-'), hold on
        ylabel('Fraction')
        set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)])
        legend({'monomer', 'smear', 'pocket'}, 'location', 'best')
        grid on
        
        
        subplot(5,1,4)
        plot(mono_spread, '.-'), hold on
        plot(mono_migrate, '.--')
        ylabel({'Normalized ', 'band width or migr. distance'})
        set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)], ...
            'YLim', [0 1.])
        grid on 
        legend({'Band width', 'Migr. distance'}, 'location', 'best')
        
        subplot(5,1,5)
        plot(folding_quality_metric, '.-'), hold on
        plot(folding_quality_metric_2, '.--')
        plot(index_best, folding_quality_metric(index_best), 'o'), hold on
        plot(index_best_Tscrn, folding_quality_metric(index_best_Tscrn), '+'), hold on
        plot(index_best_Mgscrn, folding_quality_metric(index_best_Mgscrn), 'x'), hold on

        ylabel({'Quality metric ', 'Monomer fraction/norm_width'})
        set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)])
        grid on
        legend({'Quality metric', 'Quality metric 2'}, 'location', 'best')
    end
end



