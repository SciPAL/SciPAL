a=double(imread('result_approx.tif'));
b=double(imread('result_exact.tif'));
c=double(imread('lena.tif'));
tva=(circshift(a,[1 0])-a-(circshift(c,[1 0])-c)).^2+(circshift(a,[0 1])-a-(circshift(c,[0 1])-c)).^2;
tve=(circshift(b,[1 0])-b-(circshift(c,[1 0])-c)).^2+(circshift(b,[0 1])-b-(circshift(c,[0 1])-c)).^2;
tva=mean(mean(tva));
tve=mean(mean(tve));
l2e=mean(mean((b-c).^2));
l2a=mean(mean((a-c).^2));
disp(['Exact error ' num2str(l2e+tve) ' Approx Error ' num2str(l2a+tva)])