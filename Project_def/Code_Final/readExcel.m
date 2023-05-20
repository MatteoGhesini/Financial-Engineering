function [year, DR_SG, DR_All_Rated, RR] = readExcel( filename )
%Reads data from excel
% 
%INPUTS:
% filename: excel file name where data are stored
% 
%OUTPUTS:
% year: 
% DR_SG: 
% DR_All_Rated:
% RR:
%
%USING: 
% readmatrix
%
%

    %% Dates from Excel
    
    %ESG 2021
    year = zeros(37,1);
    DR_SG = zeros(37,1);
    DR_All_Rated = zeros(37,1);
    RR = zeros(37,1);
    
    year         = readmatrix(filename,'Sheet',1,'Range','A2:A38');
    DR_SG        = readmatrix(filename,'Sheet',1,'Range','B2:B38');
    DR_All_Rated = readmatrix(filename,'Sheet',1,'Range','C2:C38');
    RR           = readmatrix(filename,'Sheet',1,'Range','D2:D38');


end % readExcel