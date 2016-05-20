% Filename    : zeropad.m
% Author      : Liu Yuan
% Email       : i@llonely.com
% =============================================================================
% Description :
% zero-pad image_1 and image_2
% in order to make that size(result_1)==size(result_2)
% Arguments:
% image_1,image_2 -> images to be zero-padded
% Return padded images
function [result_1,result_2]=zeropad(image_1,image_2)
	[height_1,width_1]=size(image_1);
	[height_2,width_2]=size(image_2);
	if height_1>height_2
		result_2=zeros(height_1,width_2);
		height_start=ceil((height_1-height_2)/2);
		result_2(height_start:height_start+height_2-1,:)=image_2;
		result_1=image_1;
	elseif height_1<height_2
		result_1=zeros(height_2,width_1);
		height_start=ceil((height_2-height_1)/2);
		result_1(height_start:height_start+height_1-1,:)=image_1;
		result_2=image_2;
	else
		result_1=image_1;
		result_2=image_2;
	end
	if width_1>width_2
		temp=zeros(size(result_2,1),width_1);
		width_start=ceil((width_1-width_2)/2);
		temp(:,width_start:width_start+width_2-1)=result_2;
		result_2=temp;
	elseif width_1<width_2
		temp=zeros(size(result_1,1),width_2);
		width_start=ceil((width_2-width_1)/2);
		temp(:,width_start:width_start+width_1-1)=result_1;
		result_1=temp;
	end