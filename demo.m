clear,close all,clc;

load data.mat;
% Write the covariance matrix in the vector form
img(:,:,1)=C11;
img(:,:,2)=C22;
img(:,:,3)=C33;
img(:,:,4)=real(C12);
img(:,:,5)=real(C13);
img(:,:,6)=real(C23);
img(:,:,7)=imag(C12);
img(:,:,8)=imag(C13);
img(:,:,9)=imag(C23);
% Filtering
ENL=2.97;
tic;
img_filtered = POLSAR_POSSC(img,'ENL',ENL);
toc;
% Get the filtered covariance matrix
C11_f=img_filtered(:,:,1);
C22_f=img_filtered(:,:,2);
C33_f=img_filtered(:,:,3);
C12_f=img_filtered(:,:,4)+1i*img_filtered(:,:,7);
C13_f=img_filtered(:,:,5)+1i*img_filtered(:,:,8);
C23_f=img_filtered(:,:,6)+1i*img_filtered(:,:,9);

% The filtering result is displayed with Pauli decomposition coefficients. 
img_pauli=fPauliImShow(C2T(img));
img_filtered_pauli=fPauliImShow(C2T(img_filtered));
figure,imshow([img_pauli,img_filtered_pauli]);