clc ,
% clear all;
close all; format long

addpath(genpath('FastMarching'))

global use_cropping sigma1 xi theta_tai use_antimarkers CV_type eps2 lambda3 edit_image lambda_TV

sigma1 = 0;
use_antimarkers = 0;
edit_image = 1;
use_cropping =1;
xi = 0.1;

u0 = imread('Tumour3.jpg');

n = 512/2;

a = size(u0); if length(a)>=3, u0 = rgb2gray(u0); end
u0 = imresize(u0,[n n]); z = double(u0);

z = (z-min(z(:)))/(max(z(:))-min(z(:)));
if sigma1~=0
    z = imnoise(z,'gaussian',0,sigma1);
end
z_star = z;

Iters = 100000;

tau = 1e-2; % time step
utol = 1e-4; % Stopping criteria
c1 = 0; c2 = 1; % Intensity Constants
output = 517; % for file save

theta_store = [];
TC_store = [];
Iter_store = [];
Res_store =  struct();
u_stored = struct();
kk3 = 1;

lambda3 = 1;
lambda_TV = 0;
theta_tai = 0;

for CV_type = [1] % 0 =old way, 1=new way, 2 = LBF, 3= Hybrid, 4 = LCV, 5-Genave+new way
    for lambda = [3]
        for theta = [3]
            for eps2 = [1e-6]
                
                fprintf('-----TEST lambda = %1.2e --- theta = %1.2e \n', lambda, theta)
                output = output + 1;
                
                lambda3_store(kk3) = lambda3;
                lambda_store(kk3) = lambda;
                theta_store(kk3) = theta;
                eps2_store(kk3) = eps2;
                CV_store(kk3) = CV_type;
                
                [TC_store(kk3),Iter_store(kk3),Res_store(kk3).h,u_stored(kk3).h,cols,rows,u_data,R_min,R_max,C_min,C_max] =...
                    ConvexSeg_Run(output,z_star,c1,c2,lambda,tau,Iters,utol,theta);
                
                kk3 = kk3 + 1;
                
                close all
            end
        end
        
    end
end

if use_cropping
    new_u = zeros(size(z));
    new_u(R_min:R_max,C_min:C_max) = (u_stored(1).h>0.5);
else
    new_u = (u_stored(1).h>0.5);
end

figure;imagesc(z);colormap(gray);axis off;hold on;contour(new_u,[0.5,0.5],'r','LineWidth',2)

