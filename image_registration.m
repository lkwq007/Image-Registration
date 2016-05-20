% Filename    : image_registration.m
% Author      : Liu Yuan
% Email       : i@llonely.com
% =============================================================================
% Description :
% An FFT-Based Technique for Translation, Rotation, and Scale-Invariant Image Registration
% clc;
clear all;
close all;
dbstop if error;
% get ref image
image_ref=imread('img.jpg');
if ndims(image_ref)==3
    image_ref=rgb2gray(image_ref);
end
image_ref_backup=image_ref;
% get the changed image
theta_ch=139;
dx_ch=21;
dy_ch=67;
scaled_ch=1.3;
crop_x=30;
crop_y=16;
crop_width=177;
crop_height=159;
noise_type='poisson';
image_ch=image_ref;
image_ch=imcrop(image_ch,[crop_x,crop_y,crop_width,crop_height]);
image_ch=imtranslate(image_ch,[dx_ch,dy_ch]);
image_ch=imresize(image_ch,scaled_ch);
image_ch=imrotate(image_ch,theta_ch,'bilinear','crop');
if ~strcmp('none',noise_type)
    image_ch=imnoise(image_ch,noise_type);
end
% image_ch=imread('anotherfile');
% image_ch=rgb2gray(image_ch);
image_ch_backup=image_ch;
tic;
[image_ref,image_ch]=zeropad(image_ref,image_ch);

% logpol size
logpol_m=size(image_ref,1);
logpol_n=size(image_ref,2);

% FFT and shift them to the center
fft_ref=fft2(image_ref);
fft_ref=fftshift(fft_ref);
fft_ch=fft2(image_ch);
fft_ch=fftshift(fft_ch);

% Highpass Filtering
mag_ref=highpass_filter(size(image_ref,1),size(image_ref,2)).*abs(fft_ref);
mag_ch=highpass_filter(size(image_ch,1),size(image_ch,2)).*abs(fft_ch);

% Cartesian to Logarithmic Polar Conversion
center=floor((size(image_ref)+1)/2);
[logpol_ref,rho,base]=logpol_transform(mag_ref,center,logpol_m,logpol_n,'bilinear','valid');
center=floor((size(image_ch)+1)/2);
logpol_ch=logpol_transform(mag_ch,center,logpol_m,logpol_n,'bilinear','valid');

% Phase Correlation
fft_logpol_ref=fft2(logpol_ref);
fft_logpol_ch=fft2(logpol_ch);
% 2 ways to get the cross power spectrum
fft_phase=(fft_logpol_ref.*conj(fft_logpol_ch)./(abs(fft_logpol_ref.*fft_logpol_ch)));
% or like this
% fft_phase=exp(1i*(angle(fft_logpol_ref)-angle(fft_logpol_ch)));
phase=real(ifft2(fft_phase));% reserve the real component
% Get the angle and scaled
phase_max=max(max(phase));
% disp(['Phase peak: ',num2str(phase_max)]);
[peak_x,peak_y]=find(phase==phase_max);
degrees_per_pixel=360/size(phase,2);
theta=(peak_y-1)*degrees_per_pixel;
% disp(['Theta:',num2str(theta)]);
% disp(abs(theta-theta_ch)/theta);
scaled_1=rho(peak_x);
% disp(scaled_1);
scaled_2=rho(logpol_m+1-peak_x);
% disp(1/scaled_2);
if(scaled_1>scaled_2)
	scaled=1/scaled_2;
else
	scaled=scaled_1;
end
% disp(['Scaled: ',num2str(scaled)]);

% try those theta 
image_ref_temp=image_ref_backup;
image_ch_temp=image_ch_backup;
if scaled>1
    image_ref_temp=imresize(image_ref_temp,scaled);
else
    image_ch_temp=imresize(image_ch_backup,1/scaled);
end
[image_ref_temp,image_ch_temp]=zeropad(image_ref_temp,image_ch_temp);
image_rech_1=imrotate(image_ch_temp,-theta,'bilinear','crop');
image_rech_2=imrotate(image_ch_temp,-theta-180,'bilinear','crop');
fft_ref_temp=fft2(image_ref_temp);
fft_ref_temp=fftshift(fft_ref_temp);
fft_rech_1=fft2(image_rech_1);
fft_rech_1=fftshift(fft_rech_1);
fft_rech_2=fft2(image_rech_2);
fft_rech_2=fftshift(fft_rech_2);
% 2 ways to get the cross power spectrum
fft_phase_rech_1=(fft_ref_temp.*conj(fft_rech_1)./(abs(fft_ref_temp.*fft_rech_1)));
fft_phase_rech_2=(fft_ref_temp.*conj(fft_rech_2)./(abs(fft_ref_temp.*fft_rech_2)));
% or like this
% fft_phase_rech_1=exp(1i*(angle(fft_ref_temp)-angle(fft_rech_1)));
% fft_phase_rech_2=exp(1i*(angle(fft_ref_temp)-angle(fft_rech_2)));

% Get the suitest theta
phase_rech_1=real(ifft2(fft_phase_rech_1));
phase_rech_2=real(ifft2(fft_phase_rech_2));
max_rech_1=max(max(phase_rech_1));
max_rech_2=max(max(phase_rech_2));
if max_rech_1>max_rech_2
	[dy,dx]=find(phase_rech_1==max_rech_1);
	phase_rech=phase_rech_1;
	image_rech=image_rech_1;
else
	[dy,dx]=find(phase_rech_2==max_rech_2);
	if(theta<180)
		theta=theta+180;
	else
		theta=theta-180;
	end
	phase_rech=phase_rech_2;
	image_rech=image_rech_2;
end
dx=dx-1;
% disp(dx)
% disp(size(phase_rech,1)-dx)
% disp(dx);
dy=dy-1;
% disp(dy)
% disp(size(phase_rech,2)-dy)
% disp(dy);

% Get the suitest tranlation
possible_tr=[-dx,-dy;-dx,size(phase_rech,2)-dy;size(phase_rech,1)-dx,-dy;size(phase_rech,1)-dx,size(phase_rech,2)-dy];
suitest=1;
max_tr=0;
for i=1:4
	dx_tr=-possible_tr(i,1);
	dy_tr=-possible_tr(i,2);
	temp=imtranslate(image_rech,[dx_tr,dy_tr]);
	fft_retr_temp=fft2(temp);
	fft_retr_temp=fftshift(fft_retr_temp);
	fft_phase_retr=(fft_ref_temp.*conj(fft_retr_temp)./(abs(fft_ref_temp.*fft_retr_temp)));
	phase_retr=real(ifft2(fft_phase_retr));
	max_tr_temp=max(max(phase_retr));
	if(max_tr_temp>max_tr)
		max_tr=max_tr_temp;
		image_retr=temp;
		suitest=i;
	end
end
toc;
dx=possible_tr(suitest,1);
dy=possible_tr(suitest,2);
% Show the results
disp(['Noise: ',noise_type]);
disp(['Phase peak: ',num2str(phase_max)]);
disp(['Scaled: ',num2str(scaled),', Real: ',num2str(scaled_ch),', Ratio: ',num2str(abs((scaled-scaled_ch)/scaled_ch))]);
disp(['Theta:',num2str(theta),', Real: ',num2str(theta_ch),', Ratio: ',num2str(abs((theta-theta_ch)/theta_ch))]);
disp(['Peak tr: ',num2str(max_tr)]);
disp(['dx: ',num2str(possible_tr(suitest,1)),',dy: ',num2str(possible_tr(suitest,2))]);
disp(['dx_scaled: ',num2str(possible_tr(suitest,1)/scaled),',dy_scaled: ',num2str(possible_tr(suitest,2)/scaled)]);
disp(['Real dx: ',num2str(dx_ch),', Ratio: ',num2str(abs((dx/scaled-dx_ch)/dx_ch)),', Real dy: ',num2str(dy_ch),', Ratio: ',num2str(abs((dy/scaled-dy_ch)/dy_ch))]);
% Show the results
figure;
subplot(2,2,1);
imshow(image_ref_backup,[]);
title('reference image');
subplot(2,2,2);
imshow(image_ch,[]);
title('image needs processing');
subplot(2,2,3);
imshow(image_retr,[]);
title('image retranslated & etc.');
subplot(2,2,4);
imshowpair(image_ref_temp,image_retr);
title('registrated image');
