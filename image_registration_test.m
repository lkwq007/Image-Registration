% Filename    : image_registration_test.m
% Author      : Liu Yuan
% Email       : i@llonely.com
% =============================================================================
% Description :
% Testing the image_registration_func.m
clc;
clear all;
close all;

image_ref=imread('img.jpg');
image_ref=rgb2gray(image_ref);
tl_step=50;
theta_step=60;
scaled_step=0.5;
noise_type={'none','gaussian','poisson'};
dx=[-size(image_ref,1)+tl_step:tl_step:size(image_ref,1)-tl_step];
dy=dx;
theta=[0:theta_step:360-theta_step];
scaled=[scaled_step:scaled_step:5];
result_count=1;
result=zeros(length(dx)*length(dy)*length(scaled)*length(theta)*length(noise_type),15);
tic
for n=1:length(noise_type)
for i=1:length(dx)
	for j=1:length(dy)
		for k=1:length(scaled)
			for l=1:length(theta)
				
					image_ch=image_operation(image_ref,dx(i),dy(j),scaled(k),theta(l),noise_type{n});
					[r_dx,r_dy,r_theta,r_scaled,phase_max,max_tr] = image_registration_func(image_ref,image_ch);
					result(result_count,:)=[dx(i),r_dx,abs((dx(i)-r_dx)/dx(i)),dy(j),r_dy,abs((dy(j)-r_dy)/dy(j)),...
					theta(l),r_theta,abs((theta(l)-r_theta)/theta(l)),scaled(k),r_scaled,abs((scaled(k)-r_scaled)/scaled(k)),...
					n,phase_max,max_tr];
					if abs((dx-r_dx)/dx)>0.2&abs((dy-r_dy)/dy)>0.2&abs((theta-r_theta)/theta)>0.2&abs((scaled-r_scaled)/scaled)>0.2
                        break;
                    end
                    if phase_max==-1&max_tr==-1
                        break;
                    end
                    result_count=result_count+1;
				end
			end
		end
	end
end
toc
save('result.mat',result);