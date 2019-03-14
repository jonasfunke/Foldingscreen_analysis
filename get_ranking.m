function [r] = get_ranking(values)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [values_sort, i_sort] = sort(values, 'descend'); 
    r = 1:length(values);
    r(i_sort) = r;
        
end

