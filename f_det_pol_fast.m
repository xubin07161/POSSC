function [ dataC_det ] = f_det_pol_fast( dataC )
%
[m,n,pag]=size(dataC);
if pag==9
    dataC_det=dataC(:,:,1).*dataC(:,:,2).*dataC(:,:,3)-...
        dataC(:,:,1).*(dataC(:,:,6).^2+dataC(:,:,9).^2)-...
        dataC(:,:,3).*(dataC(:,:,4).^2+dataC(:,:,7).^2)-...
        dataC(:,:,2).*(dataC(:,:,5).^2+dataC(:,:,8).^2)+...
        2*(dataC(:,:,4).*dataC(:,:,5).*dataC(:,:,6)+dataC(:,:,6).*dataC(:,:,7).*dataC(:,:,8)+...
        dataC(:,:,4).*dataC(:,:,8).*dataC(:,:,9)-dataC(:,:,5).*dataC(:,:,7).*dataC(:,:,9));
end
if pag==4
    dataC_det=dataC(:,:,1).*dataC(:,:,2)-(dataC(:,:,3).^2+dataC(:,:,4).^2);
end
if pag==1
    dataC_det=dataC;
end
% figure,imshow(dataT_det./mean2(dataT_det));

end

