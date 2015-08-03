clear
clc
close all
beep off
%orig=double(imread('sofi_result.analysis_orig.tif'));
orig=double(imread('lena.tif'));
orig2=double(imread('noise_compare.tif'));
%orig=orig(:,:,2);
den=double(imread('noise.tif'));
den2=double(imread('result_exact.tif'));
%den=den(:,:,2);
forig=fft2(orig);
forig2=fft2(orig2);
fden=fft2(den);
fden2=fft2(den2);

[xx yy]=meshgrid(1:size(forig,1),1:size(forig,2));
xx=min(xx,size(forig,1)-xx);
yy=min(yy,size(forig,2)-yy);
r=sqrt(xx.^2+yy.^2);
r=floor(r)+1;
rmax=max(max(r));
psdo=zeros(rmax,1);
psdo2=zeros(rmax,1);
psdd=zeros(rmax,1);
psdd2=zeros(rmax,1);
forig=abs(forig(:));
forig2=abs(forig2(:));
fden=abs(fden(:));
fden2=abs(fden2(:));
r=r(:);
for rr=1:rmax
    psdo(rr)=mean(forig(r==rr).^2);
    psdo2(rr)=mean(forig2(r==rr).^2);
    psdd(rr)=mean(fden(r==rr).^2);
    psdd2(rr)=mean(fden2(r==rr).^2);
end
psdo=psdo/sum(psdo);
psdo2=psdo2/sum(psdo2);
psdd=psdd/sum(psdd);
psdd2=psdd2/sum(psdd2);
loglog(1:rmax,psdo,1:rmax,psdo2,1:rmax,psdd,1:rmax,psdd2)
legend('Original Image','Noisy Image','Approx','Exact')
grid on
xlabel('|k|_2')
ylabel('power spectral density')

%%
x = ones(length(psdd), 1);
for n = 1:length(psdd)
    x(n) = n;
end

A = [x, psdo, psdo2, psdd, psdd2];
dlmwrite('psd_lena_data', A, '\t')
