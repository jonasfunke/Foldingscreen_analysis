function [gelData, gelInfo, profileData] = compute_profiles(root_path, dir_name, gel_info_name, image_name)
% @ step1
% calculate profiles from gels

    path_out = [root_path filesep dir_name ];
    file_out = [path_out filesep dir_name];

    % load gel_info and check for errors
    [gelInfo, ~] = parse_gel_info([path_out filesep gel_info_name], [path_out filesep gel_info_name(1:end-4) '_log.log']);
    check_parsed_gel_info(gelInfo);

    %% via MATLAB_TOOLBOX
    % load gelData
    gelData_raw = load_gel_image('data_dir', path_out, 'n_images', 1,...
        'paths_to_images', {[path_out filesep image_name]});
    % check and correct raw data
    gelData_raw = check_gel_saturation(gelData_raw);
    gelData = background_correct_gel_image(gelData_raw, 'histogram_background', 'on');
    % get profileData
    profileData = get_gel_lanes(gelData, 'display', 'on', 'cutoff', 0.05, ...
        'number_of_lanes', length(gelInfo.lanes), 'display', 'off', 'move_min_zero', 'Yes');

    
    %% check and link data
    if length(profileData.fullProfiles)==length(gelInfo.lanes)

        % plot areas 
        cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','points','PaperPosition', [0 0 1000 500], 'Position', [0 1000 1000 500]);
        imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on

        areas = [profileData.lanePositions(:,1)   profileData.lanePositions(:,3) profileData.lanePositions(:,2)-profileData.lanePositions(:,1)  profileData.lanePositions(:,4)-profileData.lanePositions(:,3)];
        for i=1:length(profileData.profiles)
            rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
            text(areas(i,1)+areas(i,3)/2, areas(i,2) , gelInfo.lanes{i}, 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
        end
        set(gca, 'XTickLabel', [], 'YTickLabel', [])
        print(cur_fig, '-dpng', '-r 300' , [file_out '_lanes.png']); %save figure

        % plot profiles
    
        cur_fig = figure;
        set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
            'PaperPosition', [0 0 20 3*length(profileData.profiles)], ...
            'PaperSize', [ 20 3*length(profileData.profiles)]);
        
        for i=1:length(profileData.profiles)
            subplot(length(profileData.profiles), 1, i)
            plot(profileData.lanePositions(i,3):profileData.lanePositions(i,4), profileData.profiles{i}), hold on
            legend(gelInfo.lanes{i})
            %myleg = [myleg, {gelInfo.lanes{i}}];
        end
        %legend(myleg)
        ylabel('Raw Intensity')
        xlabel('Migration distance [px]')
        print(cur_fig, '-dpdf', [file_out '_profiles.pdf']); %save figure

    else
        disp('Warning: number of lanes in gel_info.txt does noth match number of detected lanes.' )
    end

    %% save data
    disp('Saving data... please wait')
    save([file_out  '_data.mat'], 'gelData', 'gelInfo', 'profileData' )
    disp(['Data written to: ' file_out  '_data.mat'])

end

