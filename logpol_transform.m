% Filename    : logpol_transform.m
% Author      : Liu Yuan
% Email       : i@llonely.com
% =============================================================================
% Description :
% Cartesian to Logarithmic-Polar Conversion
% Arguments: 
% original -> original image or matrix
% center   -> origin of the input matrix
% n_rho    -> the desired number of rows of transformed matrix
% n_theta  -> the desired number of columns of transformed matrix
% method   -> interpolation method
% shape    -> output size
% Return the transformed matrix
function [result,rho,base]=logpol_transform(original,center,n_rho,n_theta,method,shape)
	[m,n]=size(original);
	theta=linspace(0,2*pi,n_theta+1);
	theta(end)=[];
	switch(shape)
	case 'full'
		mat_corner=[1,1;m,1;m,n;1,n];
		base=max(sqrt(sum((repmat(center(:)',4,1)-mat_corners).^2,2)));
	case 'valid'
		base=min([n-center(1),center(1)-1,m-center(2),center(2)-1]);
	otherwise
		error(['Arguments not match: shape should be valid or full, not ',shape]);
	end
	minimum=1;
	rho=logspace(log10(minimum),log10(base),n_rho)';
	x=rho*cos(theta)+center(1);
	y=rho*sin(theta)+center(2);
	switch(method)
	case 'nearest'
		result=interp2(original,x,y,'nearest');
	case 'bilinear'
		result=interp2(original,x,y,'linear');
	case 'bicubic'
		result=interp2(original,x,y,'cubic');
	otherwise
		%error(['Arguments not match: method should be nearest or bilinear or bicubic, not ',method]);
		result=interp2(original,x,y,'linear');
	end
	outside=(x>n)|(x<1)|(y>m)|(y<1);
	result(outside)=0;