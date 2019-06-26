%% execute this to select a folder and analyze a folding screen

close all, clear all, clc
root_path = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS';

[path_selected] = uigetdir(root_path, 'Select an Directory with intial foling screen');
[~,dir_name,ext] = fileparts(path_selected);

if ~isempty(ext)
    dir_name = [dir_name ext];
end

txt_files = dir([path_selected filesep '*.txt']);
tif_files = dir([path_selected filesep '*.tif']);

if length(tif_files)==1 && length(txt_files)==1
    compute_profiles(root_path, dir_name, txt_files(1).name, tif_files(1).name)
else
    disp(['Error. Found ' num2str(length(tif_files)) ' tif files and ' num2str(length(txt_files)) ' txt files in ' dir_name])
end

%% analyze the profiles

close all, clear all, clc
[fname pname] = uigetfile('/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/*.mat', 'Select mat file');

data = load([pname fname]); % load data


if ~isfield(data.profileData, 'aggregateSelectedArea')
    % select aggregates and monomer bands
    data.profileData = select_species(data.profileData, data.gelData, 1.5); 

    % integrate aggregates, smear, and monomer bands
    data.profileData = integrate_species(data.profileData, 1.5);
    %get_best_folding(data.profileData, data.gelInfo, data.gelData, true);
    
    save([pname fname], '-struct','data')
    disp(['Data saved to ' pname fname])

else
    display('Pocket and monomers have already been selected.')
    


end

[~, cur_fig] = get_best_folding(data.profileData, data.gelInfo, data.gelData, true);
tmp = strsplit(pname, filesep);
print(cur_fig, '-dpdf', [pname tmp{end-1} '_analysis.pdf']); %save figure

%% -------------------- Analyze all screens -----------------------------------------

%log_out = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/data.out';
root_path = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/';

%logfile_ID = fopen(log_out,'a');
%fprintf(logfile_ID,'%s\n', '------');
%fclose(logfile_ID);

%logfile_ID = fopen(log_out,'a');


folders = dir(root_path);

discard = {'.', '..', 'AAA_TEMPLATE', 'ZZZfolder_of_shame_aka_missing_data', 'ZZZnon_standard_folding_screens'};

i_discard = [];
for i=1:length(folders)
    if any(strcmp(discard, folders(i).name))
        disp(['Discarding folder ' folders(i).name])
        i_discard = [i_discard, i];
    end
    if ~folders(i).isdir
        i_discard = [i_discard, i];
    end
end
folders(i_discard) = [];

metrics = cell(length(folders),1);
for i=1:length(folders)
    if folders(i).isdir
        disp('---')
        %fprintf(logfile_ID,'%s\n', '------');
        txt_files = dir([folders(i).folder filesep folders(i).name filesep '*.txt']);
        tif_files = dir([folders(i).folder filesep folders(i).name filesep '*.tif']);
        mat_files = dir([folders(i).folder filesep folders(i).name filesep '*.mat']);
        
        if length(txt_files)==1
           % [tmp, warnings] = parse_gel_info([txt_files(1).folder filesep txt_files(1).name], log_out);
           %     check_parsed_gel_info(tmp);
        else
            disp(['More or less than one txt file found in ' folders(i).name ' (' num2str(length(txt_files)) ' found).' ])
            %fprintf(logfile_ID,'%s\n', ['More or less than one txt file found in ' folders(i).name ' (' num2str(length(txt_files)) ' found).' ]);
        end
        
        if length(tif_files)==1
            %[tmp, warnings] = parse_gel_info_simple([txt_files(1).folder filesep txt_files(1).name]);
            %    check_parsed_gel_info(tmp);
        else
            disp(['More or less than one tif file found in ' folders(i).name ' (' num2str(length(tif_files)) ' found).' ])
            %fprintf(logfile_ID,'%s\n', ['More or less than one tif file found in ' folders(i).name ' (' num2str(length(tif_files)) ' found).' ]);
        end
        
        if length(mat_files)==1
            %[tmp, warnings] = parse_gel_info_simple([txt_files(1).folder filesep txt_files(1).name]);
            %    check_parsed_gel_info(tmp);
            cur_data = load([mat_files(1).folder filesep mat_files(1).name]);
            if isfield(cur_data.profileData, 'monomerFits')
                metrics{i} = get_best_folding(cur_data.profileData, cur_data.gelInfo, cur_data.gelData, false);
            end
            
        else
            disp(['More or less than one .mat file found in ' folders(i).name ' (' num2str(length(mat_files)) ' found).' ])
            %fprintf(logfile_ID,'%s\n', ['More or less than one tif file found in ' folders(i).name ' (' num2str(length(tif_files)) ' found).' ]);
        end
       
       
       
    end
end
%fclose(logfile_ID);


d = [];
for i=1:length(metrics)
    if ~isempty(metrics{i})
        d = [d; metrics{i}.fractionMonomer(metrics{i}.bestFoldingIndex)];
    end
end
%%

