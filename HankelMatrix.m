function [X1, X2] = HankelMatrix(Data,n,m, type)
%HankelBlock - Generates Hankel blocks from input data matrix.
%
%   [X1, X2] = HankelMatrix(Data, n, m, type)
%
%   This function takes input data matrix 'Data' and generates Hankel blocks
%   based on the specified parameters 'n', 'm', and 'type'.
%
%   Parameters:
%   - Data (matrix): Input data matrix.
%   - n (integer): Number of cols for the Hankel matrix.
%   - m (integer): Number of rows for the Hankel matrix.
%   - type (string): Type of block generation ('row' or 'column').
%
%   Returns:
%   - X1 (matrix): First set of generated Hankel blocks.
%   - X2 (matrix): Second set of generated Hankel blocks.
% Adapted from "Data-driven spectral analysis for coordinative structures 
% in periodic human locomotion" Fujii (2019)

index1 = 1:n;
index2 = n:n+m-1;

X1 = []; X2=[];

for ir = 1:size(Data,1)
   
    % Hankel blocks ()
    c = Data(ir,index1).'; r = Data(ir,index2);
    H = hankel(c,r).';
    c = Data(ir,index1+1).'; r = Data(ir,index2+1);
    UH= hankel(c,r).';
    
    % input data matrices for exact DMD
    if strcmp(type,'row')
        X1=[X1,H]; X2=[X2,UH];
    elseif strcmp(type,'column')
        X1=[X1;H]; X2=[X2;UH];
    end
    
end
end