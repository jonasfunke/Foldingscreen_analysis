function [data_out] = get_best_folding(profileData, gelInfo, gelData, show_summary_figure)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    %%
    %keyboard
    %%
    height_width = zeros(length(profileData.profiles),1);
    width = zeros(length(profileData.profiles),1);
    for i=1:length(profileData.profiles)
        height_width(i) = profileData.monomerFits{i}.a1/profileData.monomerFits{i}.c1;
        width(i) = profileData.monomerFits{i}.c1;
    end
        
    % T-screen indices
    index_Tscrn = [];
    for i=1:length(gelInfo.lanes)
        cur_name = strtrim(gelInfo.lanes{i});
        if strcmpi(cur_name(1), 'T')
            index_Tscrn(str2num(cur_name(2))) = i;
        end
    end
    % Mg-screen indices
    index_Mgscrn =[];
    for i=1:length(gelInfo.lanes)
        cur_name = strtrim(gelInfo.lanes{i});
        if strcmpi(cur_name, 'M5')
            index_Mgscrn(1) = i;
        end
        if strcmpi(cur_name, 'M10')
            index_Mgscrn(2) = i;
        end
        if strcmpi(cur_name, 'M15')
            index_Mgscrn(3) = i;
        end
        if strcmpi(cur_name, 'M20')
            index_Mgscrn(4) = i;
        end
        if strcmpi(cur_name, 'M25')
            index_Mgscrn(5) = i;
        end
        if strcmpi(cur_name, 'M30')
            index_Mgscrn(6) = i;
        end
    end
            
    % RM indices
    index_RM = [];
    for i=1:length(gelInfo.lanes)
        cur_name = strtrim(gelInfo.lanes{i});
        if strcmpi(cur_name, 'RM1')
            index_RM(1) = i;
        end
        if strcmpi(cur_name, 'RM1_diluted')
            index_RM(1) = i;
        end
        if strcmpi(cur_name, 'RM2')
            index_RM(2) = i;
        end
    end 
    % all indices except for scaffold and ladder
    index_foldings = [];
    for i=1:length(gelInfo.lanes)
        cur_name = strtrim(gelInfo.lanes{i});
        if contains(cur_name,'scaffold','IgnoreCase',true) || strcmpi(cur_name, '1kb_ladder') || contains(cur_name,'ladder','IgnoreCase',true)
            
        else
            index_foldings = [index_foldings, i];
        end
    end 
    
    % index M20
    index_M20 = [];
    for i=1:length(gelInfo.lanes)
        cur_name = strtrim(gelInfo.lanes{i});
        if strcmpi(cur_name, 'M20')
            index_M20 = i;
        end
    end 
    
    %index scaffold
    index_scaffold = [];
    for i=1:length(gelInfo.lanes)
        cur_name = strtrim(gelInfo.lanes{i});
        if contains(cur_name, 'sc','IgnoreCase',true)
            index_scaffold = [index_scaffold, i];
        end
    end 
    
   %keyboard
   %%
   
    migration_distance = zeros(length(gelInfo.lanes),1);
    for i = 1:length(gelInfo.lanes)
        migration_distance(i) = profileData.monomerFits{i}.b1;
    end
    
    normalized_migration_distance = (migration_distance-min(migration_distance(index_foldings)))./(max(migration_distance(index_foldings))-min(migration_distance(index_foldings)));
    %%
   folding_quality_metric_old = profileData.monomerTotal./(profileData.monomerTotal+profileData.pocketTotal+profileData.smearTotal)./width;
   folding_quality_metric = normalized_migration_distance.*profileData.monomerTotal./(profileData.monomerTotal+profileData.pocketTotal+profileData.smearTotal)./width;
   %%
   
   
   
   % find best folding
   [~, i_sort] = sort(folding_quality_metric(index_foldings), 'descend');
   index_best = index_foldings(i_sort(1));
   disp(['Best folding: Lane ' num2str(index_best) ' (' gelInfo.lanes{index_best} ')'])
    %get_ranking(height_width(index_foldings))
    
   % find best Temperature interval
   if ~isempty(index_Tscrn)
       [~, i_sort] = sort(folding_quality_metric(index_Tscrn), 'descend');
       index_best_Tscrn = index_Tscrn(i_sort(1));
       disp(['Best T-Screen: Lane ' num2str(index_best_Tscrn) ' (' gelInfo.lanes{index_best_Tscrn} ')'])
   end
   % find best Mg concentratio
   if ~isempty(index_Mgscrn)
       index_Mgscrn = index_Mgscrn(index_Mgscrn~=0); % remove zeros if people did not include all Mg samples
       [~, i_sort] = sort(folding_quality_metric(index_Mgscrn), 'descend');
       index_best_Mgscrn = index_Mgscrn(i_sort(1));
       disp(['Best Mg-Screen: Lane ' num2str(index_best_Mgscrn) ' (' gelInfo.lanes{index_best_Mgscrn} ')'])
   end
   
   % find best staple-scaffold ratio
   if ~isempty(index_RM)
       index_RM = index_RM(index_RM~=0); % remove zeros if people did not include all RM samples
       [~, i_sort] = sort(folding_quality_metric(index_RM), 'descend');
       index_best_RM = index_RM(i_sort(1));
       disp(['Best RM: Lane ' num2str(index_best_RM) ' (' gelInfo.lanes{index_best_RM} ')'])
       data_out.bestRM = gelInfo.lanes{index_best_RM};
       data_out.bestRMIndex = index_best_RM;
   end
   
   % compare M20 to best Temp folding
   % RM indices
    
    % check if M20 is better than best 4h
    if ~isempty(index_M20) && ~isempty(index_best_Tscrn)
       if folding_quality_metric(index_M20) > folding_quality_metric(index_best_Tscrn)
           M20_better_than_bestT = true;
       else
           M20_better_than_bestT = false;
       end
       disp(['Is M20 better than best 4h folding (' gelInfo.lanes{index_best_Tscrn} '): ' num2str(M20_better_than_bestT)])

    end
   
   
    %%
   %% 
   
   % calculate amount of monomer, smear and aggreagtes for best folding
   % condition and M20
   
   fraction_monomer = profileData.monomerTotal./(profileData.monomerTotal+profileData.pocketTotal+profileData.smearTotal);
   fraction_smear = profileData.smearTotal./(profileData.monomerTotal+profileData.pocketTotal+profileData.smearTotal);
   fraction_pocket = profileData.pocketTotal./(profileData.monomerTotal+profileData.pocketTotal+profileData.smearTotal);
   disp(['Fraction of monomers/smear/pocket for the best foldings: ' num2str(round(100*fraction_monomer(index_best))) '%, ' ...
       num2str(round(100*fraction_smear(index_best))) '%, ', ...
       num2str(round(100*fraction_pocket(index_best))) '%'])
   
   disp(['Fraction of monomers/smear/pocket for M20: ' num2str(round(100*fraction_monomer(index_M20))) '%, ' ...
       num2str(round(100*fraction_smear(index_M20))) '%, ', ...
       num2str(round(100*fraction_pocket(index_M20))) '%'])
 
   
    % normalize width to scaffold width
    
    width_normalized = width/mean(width(index_scaffold));
    
    data_out.fractionMonomer = fraction_monomer;
    data_out.fractionSmear = fraction_smear;
    data_out.fractionPocket = fraction_pocket;
    data_out.bestFolding = gelInfo.lanes{index_best};
    data_out.bestFoldingIndex = index_best;
    data_out.bestTscrn = gelInfo.lanes{index_best_Tscrn}; 
    data_out.bestTscrnIndex = index_best_Tscrn; 
    data_out.bestMgscrn = gelInfo.lanes{index_best_Mgscrn};
    data_out.bestMgscrnIndex = index_best_Mgscrn;
    
    data_out.M20BetterThanTscrn = M20_better_than_bestT;
    data_out.bandWidthNormalized = width_normalized;
    
    

    %show_summary_figure = false;
    
    if show_summary_figure
    

        cur_fig = figure('units','normalized','outerposition',[0 0 1 1]);
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
    %     plot(width_normalized, '.-')
    %     xlabel('Lane')
    %     ylabel('Normalized Width')
    %     set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)])

       subplot(4,1,1:3)
       imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, colorbar, hold on

       % plot pocket fits
        for i=1:length(profileData.profiles)
            plot(mean(profileData.lanePositions(i,1:2)), profileData.aggregateFit.b1, 'r.')
            plot(mean(profileData.lanePositions(i,1:2))*[1 1], ...
                [profileData.aggregateFit.b1-profileData.sigma_integrate*profileData.aggregateFit.c1 ...
                profileData.aggregateFit.b1+profileData.sigma_integrate*profileData.aggregateFit.c1], 'r')
        end

        % plot leading band fits
        for i=1:length(profileData.profiles)
            plot([profileData.lanePositions(i,1) profileData.lanePositions(i,2) ] , ...
                [profileData.monomerFits{i}.b1 profileData.monomerFits{i}.b1], 'r')
            plot(mean(profileData.lanePositions(i,1:2))*[1 1], ...
                [profileData.monomerFits{i}.b1-profileData.sigma_integrate*profileData.monomerFits{i}.c1 ...
                profileData.monomerFits{i}.b1+profileData.sigma_integrate*profileData.monomerFits{i}.c1], 'r')
        end


        plot(mean([profileData.lanePositions(index_best,1) profileData.lanePositions(index_best,2) ]) , ...
                profileData.monomerFits{index_best}.b1, 'go')
        plot(mean([profileData.lanePositions(index_best_Tscrn,1) profileData.lanePositions(index_best_Tscrn,2) ]) , ...
                profileData.monomerFits{index_best_Tscrn}.b1, 'g+')
        plot(mean([profileData.lanePositions(index_best_Mgscrn,1) profileData.lanePositions(index_best_Mgscrn,2) ]) , ...
                profileData.monomerFits{index_best_Mgscrn}.b1, 'gx')
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

        subplot(4,1,4)
        plot(folding_quality_metric, '.-'), hold on
        plot(folding_quality_metric_old, '.--')
        plot(index_best, folding_quality_metric(index_best), 'o'), hold on
        plot(index_best_Tscrn, folding_quality_metric(index_best_Tscrn), '+'), hold on
        plot(index_best_Mgscrn, folding_quality_metric(index_best_Mgscrn), 'x'), hold on

        xlabel('Lane')
        ylabel('monomer/(monomer+pocket+smear)/width')
        set(gca, 'XTick', [1:length(profileData.profiles) ], 'XTickLabels', gelInfo.lanes, 'XLim', [1 length(profileData.profiles)])
        
        pause
        close all
    end

 
end

