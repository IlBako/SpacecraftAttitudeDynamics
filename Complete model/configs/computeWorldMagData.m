clc; clearvars; close all;
%{
    Compute the Guassian normalized form of IGRF coefficients up to 
    order 13 and pre-compute the K matrix needed for the on-line calculation of
    the Gaussian normalized form of the lagrange polynomials P and their
    partial derivatives.
    
    The IGRF coefficients are taken from [1], while the model definition and the
    recursive formula for S, P and K are taken from [2].    
    
        
    Notes: 
    -All the matrices are [n x m];
    -g and h are in nT.

    References:
    [1] International Geomagnetic Reference Field (IGRF). 
        https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field.
    
    [2] M. Navabi, M. Barati, Mathematical modeling and simulation of the earth's magnetic field: 
        A comparative study of the models on the spacecraft attitude control application,
        Applied Mathematical Modelling, Volume 46, 2017, Pages 365-381, ISSN 0307-904X,
        https://doi.org/10.1016/j.apm.2017.01.040.
%}

%% Compute S matrix

N = 13;

S_nm = zeros(N, N+1);

% First value: S_(0,0)
S_nm(1,1) = 1;

% First column: S_(n,0)
for n = 2:N
    S_nm(n,1) = S_nm(n-1, 1) * ((2*n-1)/n);
end

% Remaining components: S_(n,m)
for n = 1:N
    for j = 2:n+1
        m = j-1;    % m starts from 0 and goes to n, the index j goes from 1 to n+1
        S_nm(n, j) = S_nm(n, j-1) * sqrt((n-m+1)*((m==1)+1)/(n+m));
    end
end

%% Load g and h coefficients from igrf13

coeffs = readtable("data\igrf13coeffs.txt", 'VariableNamingRule','preserve');
% As the file contains all coefficients from various years, select only the 2020 ones
coeffs = coeffs(2:end, [1:3 end-1]);

g = zeros(n, n+1);
h = zeros(n, n+1);

for i=1:size(coeffs, 1)
    if strcmp(cell2mat(table2array(coeffs(i, 1))), 'g')
        g(table2array(coeffs(i, 2)), table2array(coeffs(i, 3))+1) = table2array(coeffs(i, 4));
    end
    if strcmp(cell2mat(table2array(coeffs(i, 1))), 'h')
        h(table2array(coeffs(i, 2)), table2array(coeffs(i, 3))+1) = table2array(coeffs(i, 4));
    end
end

%% Gaussian normalization of g and h coefficients

for i=1:N
    for j=1:N+1
        g(i,j) = S_nm(i,j) * g(i,j);
        h(i,j) = S_nm(i,j) * g(i,j);
    end
end

%% K matrix pre-calculation

K = zeros(N, N+1);

for n = 2:N
    for j = 1:n+1
        m = j-1;
        K(n, j) = ((n-1)^2-m^2)/((2*n-1)*(2*n-3));
    end
end

%% Output file generation
% As this data is independent from the orbit, it can be computed once and
% stored in a .mat file so it can be fastly loaded for the simulations

save("data\WMD.mat", 'N', 'K', 'g', 'h', '-mat');