a=imread('sofi_result.analysis.tif');
a=20*double(a(:,:,2));
imwrite(uint8(a),'sofi_result_analysis.png');
a=imread('sofi_result.analysis_orig.tif');
a=20*double(a(:,:,2));
imwrite(uint8(a),'sofi_result_analysis_orig.png');