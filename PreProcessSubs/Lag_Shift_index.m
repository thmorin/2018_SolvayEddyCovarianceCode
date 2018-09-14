function [V_lag_1, V_lag_end , V_1 , V_end] = ...
    Lag_Shift_index (lag , concatenating_files , Current_Window , Number_o_Timestamps_window , Windows_per_File)

% function [V_lag_1, V_lag_end , V_1 , V_end] = ...
%   Lag_Shift_index (lag , concatenating_files , Current_Window , Number_o_Timestamps_window , Windows_per_File)
%
% Calculates the 4 indices of the shift after calling subrouting lag_test.
% Used when NOT concatenating files. In this case we might have to insert
% nan's if this is the first window of the file and the lag is positive, OR
% if this is the last window of the file and the lag is negative.
%
% lag = the output of the subrouting lag_test
% concatenating_files = Binary variable, '0' not concatenating files , '1' concatenating files.
% Current_Window = Number of current window in file.
% Number_o_Timestamps_window = Number of timestamps per window.
% Windows_per_File = Number of windows per file.
%
% V_lag_1 = first index of the current window's data vector.
% V_lag_end = last index of the current window's data vector.
% V_1 = first index of the data source vector.
% V_end = last index of the data source vector.

    if concatenating_files == 0
        if Current_Window == 1 && lag < 0
            V_lag_1 = -lag +1;
            V_lag_end = Number_o_Timestamps_window ;
            V_1 = 1;
            V_end = Number_o_Timestamps_window + lag ;
        elseif Current_Window == Windows_per_File && lag > 0
            V_lag_1 = 1;
            V_lag_end = Number_o_Timestamps_window - lag ;
            V_1 = (Current_Window-1)*Number_o_Timestamps_window + 1 + lag ;
            V_end = Current_Window * Number_o_Timestamps_window ;
        else 
            V_lag_1 = 1;
            V_lag_end = Number_o_Timestamps_window ;
            V_1 = (Current_Window-1)*Number_o_Timestamps_window + 1 + lag ;
            V_end = Current_Window * Number_o_Timestamps_window + lag ;
        end
    else
        disp( 'Do not use Lag_Shift_index when concatenating files' )
    end
    
end