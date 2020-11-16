%% execute this to select a folder and analyze a folding screen
[pname] = uigetdir(pwd, 'Select an Directory with intial foling screen');

[fname] = step1_make_profiles(pname);
step2_anotate_profiles(fname, pname);
step3_analyse_profiles(fname, pname);
