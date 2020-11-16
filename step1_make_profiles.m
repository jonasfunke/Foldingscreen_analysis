%% execute this to select a folder and analyze a folding screen
[path_selected] = uigetdir(pwd, 'Select an Directory with intial foling screen');

step1(path_selected)

function step1(path_selected) 

    [~, dir_name, ext] = fileparts(path_selected);

    if ~isempty(ext)
        dir_name = [dir_name ext];
    end

    txt_files = dir([path_selected filesep '*.txt']);
    tif_files = dir([path_selected filesep '*.tif']);

    if length(tif_files)==1 && length(txt_files)==1
        compute_profiles(pwd, dir_name, txt_files(1).name, tif_files(1).name);
    else
        disp(['Error. Found ' num2str(length(tif_files)) ' tif files and ' num2str(length(txt_files)) ' txt files in ' dir_name])
    end

    disp('Done.')
end