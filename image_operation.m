% Filename    : image_operation.m
% Author      : Liu Yuan
% Email       : i@llonely.com
% =============================================================================
% Description :
% RGB to Gray
% Tranlate, resize, rotate the image, and add noise
% Argument -> image path or image, dx, dy, scaled, theta, noise_type
% return image & changed image
% function [image_ori,image_ch] = image_operation(image_path,dx,dy,scaled,theta,noise_type,varargin)
% 	if ischar(image_path)
% 		image_ori=imread(image_path);
% 	else
% 		image_ori=image_path;%for optimization
% 	end
% 	if ndims(image_ori)==3
% 		image_ori=rgb2gray(image_ori);
% 	end
% Rewrite for optimization
% Argument -> image, dx, dy, scaled, theta, noise_type
% return changed image
function [image_ch] = image_operation(image_ori,dx,dy,scaled,theta,noise_type,varargin)
	image_ch=image_ori;
	if(nargin>6)
		x_min=varargin{1};
		y_min=varargin{2};
		crop_width=varargin{3};
		crop_height=varargin{4};
		image_ch=imcrop(image_ch,[x_min,y_min,crop_width,crop_height]);
	end
	image_ch=imtranslate(image_ch,[dx,dy]);
	image_ch=imresize(image_ch,scaled);
	image_ch=imrotate(image_ch,theta,'bilinear','crop');
	if ~strcmp(noise_type,'none')
		image_ch=imnoise(image_ch,noise_type);
	end