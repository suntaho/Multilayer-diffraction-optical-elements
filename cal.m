% this program calculates masks for ultra-high-resolution square diffracction images
% version: Matlab R2017a; intel fortran 2017; windows 10
clear;clc;close all;

% user-defined parameters
ncpu = 4;                                                                  % number of threads
nr = 2;                                                                    % 1/2/3; multiple of resolution
mwl = 0;                                                                   % 0:Off, 1:On. switch multi-wavelength optimization 
phi_var = 10.0;                                                            % percentage; if mwl=1, set two wavelengths phi+phi_var and phi-phi_var for optimization
fstr = 2.0;                                                                % fluctuation strength.
fprd = 2;                                                                  % fluctuation period in unit of iteration.
nitr = 7;                                                                  % number of solving iterations
loadmask=[''];
% define number of grating layer
nlay = nr*nr;                                                              % for nr=1/2/3 respectively
%define offset vector
ushift=zeros(nlay,nr); vshift=zeros(nlay,nr);
layi=0;
for w = 0:nr-1
	for s = 0:nr-1
		layi=layi+1;
		ushift(layi,nr)=w;
		vshift(layi,nr)=s;
	end
end
% define phase difference for each binary phase layer
if(nr==1)
    dphi=pi;
else
    dphi=2.0*pi/(nlay-1.0);
end
% read diffraction
imgfp=[num2str(nr) '.bmp']; img=imread(imgfp); img=mean(img,3);            % dimensions must be multiple of nr



imgx = size(img,1); imgy = size(img,2); 
maskx = imgx/nr; masky = imgy/nr;
% define arrays
imgi = img;                                                                % original images
xc = round(0.5*imgx);  yc = round(0.5*imgy);
img = circshift(img,-xc,1); img = circshift(img,-yc,2);                    % shift image
% guess initial masks
if (size(loadmask,2)==0)
    mask=round(rand(imgx,imgy,nlay));
	for i=1:nr
		for j=1:nr
			mask(i:nr:end,j:nr:end,:)=mask(1:nr:end,1:nr:end,:);
		end
	end
	% shift the mask
	for i=1:nlay
		mask(:,:,i)=circshift(mask(:,:,i),ushift(i,nr),1);
		mask(:,:,i)=circshift(mask(:,:,i),vshift(i,nr),2);
	end
else
    load(loadmask);
end 
I0 = sum(img(:));
theta = sum(mask,3)*dphi;
g = ifft2(exp(sqrt(-1.0)*theta));
I = real(g).^2+imag(g).^2;
alpha = I0/sum(I(:));
bufm = (img-alpha*I).^2;
Etold = sum(bufm(:));

% export variables
vars = [ncpu nr mwl phi_var fstr fprd nitr nlay imgx imgy maskx masky dphi I0 Etold];
fid = fopen('tmp1','w');   % export vars
for i = 1:size(vars,2)
    fprintf(fid,'%f\n',vars(1,i));
end
for i = 1:size(ushift,1)
    fprintf(fid,'%g\t',ushift(i,:));
    fprintf(fid,'\n');
end
for i = 1:size(vshift,1)
    fprintf(fid,'%g\t',vshift(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
fid2 = fopen('tmp2','w');   % export img
for i = 1:size(img,1)
    fprintf(fid2,'%g\t',img(i,:));
    fprintf(fid2,'\n');
end
fclose(fid2);
fid3 = fopen('tmp3','w');   % export g: real part
for i = 1:size(g,1)
    fprintf(fid3,'%g\t',real(g(i,:)));
    fprintf(fid3,'\n');
end
fclose(fid3);
fid4 = fopen('tmp4','w');   % export g: imaginary part
for i = 1:size(g,1)
    fprintf(fid4,'%g\t',imag(g(i,:)));
    fprintf(fid4,'\n');
end
fclose(fid4);
fid5 = fopen('tmp5','w');   % export mask
for m = 1:size(mask,3)
    maski = mask(:,:,m);
    for i = 1:size(g,1)
        fprintf(fid5,'%g\t',maski(i,:));
        fprintf(fid5,'\n');
    end
end
fclose(fid5);

% call fortran program
mex COMPFLAGS=" /O2 /QxHost /Qparallel /free /real_size:64 /Qopenmp /MD /fp:source /assume:bscc -I"D:\Program\\Files\MATLAB\R2017a\extern\include"" -compatibleArrayDims my3d.f90;
my3d();


 
% read results
errhist = textread('errhist.txt');
maskoAll = textread('masko.txt');
g0real = textread('g0real.txt');
g0imag = textread('g0imag.txt');
g0 = g0real+sqrt(-1.0)*g0imag;
masko=round(rand(imgx,imgy,nlay));
for i=1:nlay
    masko(:,:,i)=maskoAll((i-1)*imgx+1:i*imgx,:);
end
%save mask
mask=masko;
save mask.mat mask;
% show image
theta = sum(masko,3)*dphi;
g = ifft2(exp(sqrt(-1.0)*theta));
I = real(g).^2+imag(g).^2;
I = circshift(I,xc,1); I = circshift(I,yc,2);  
subplot(1,2,1); imagesc(100*I); colormap gray; axis off; box off; axis equal; title('by FFT');
subplot(1,2,2); imagesc(100*circshift(circshift(abs(g0).^2,xc,1),yc,2)); colormap gray; axis off; box off; axis equal; title('by Fortran');

