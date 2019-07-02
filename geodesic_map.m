function [ D ] = geodesic_map( name, ph, rows, cols )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
warning off
addpath(genpath('/Users/mike/Desktop/Maths/Backup 09102017/FastMarching/'));
addpath('FastMarching\');

if nargin == 1
n = max(size(name,1),size(name,2));
[M,W] = load_potential_map(name, n);
[start_points2] = pick_start_end_point(M);
start_points =[start_points2(2);start_points2(1)];
options.nb_iter_max = Inf;
[D,~] = perform_fast_marching(W, start_points, options);
figure;imagesc(D)
else
   n = max(size(name,1),size(name,2));
    [M,W] = load_potential_map(name, n,cols,rows);
    options.nb_iter_max = Inf;

rows_ref = rows;
cols_ref = cols;

start_points = [rows_ref,cols_ref]';
[D,~] = perform_fast_marching(W, start_points, options);
% figure;imagesc(D)
end

end

