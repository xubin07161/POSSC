function [IOut] = POLSAR_POSSC(dataC,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  POSSC is an algorithm for POLSAR image filtering. 
%  This algorithm reproduces the results from the article:
%  [1] B. Xu et al. 'Polarimetric SAR Image Filtering based on Patch
%      Ordering and Simultaneous Sparse Coding'
%  Please refer to this paper for a more detailed description of the algorithm.
%
%  BASIC USAGE EXAMPLES:
%
%     1) Using the default parameters
% 
%      img_filtered = POLSAR_POSSC(img,'ENL',ENL)
% 
%  INPUT ARGUMENTS (OPTIONAL):
%
%     1) img : The input POLSAR image. Each pixel should be the vector form
%     of covariance matrix.
%
%     2) ENL : The equivalent number of looks. The ENL can be obtained by 
%              by supervised or unsupervised estimation. For a homogeneous 
%              region, the ENL can be calculated by 
%              ENL=Tr(<C>)^2/(<Tr(CC)>-Tr(<C><C>)).
%
%  OUTPUTS:
%     1) img_filtered  : The filtered image. Each pixel is the vector form
%     of covariance matrix.                                              
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2014 Bin Xu.
% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% AUTHORS:
%     Bin Xu, email: xubin07161@gmail.com
%     Bin Zuo, email: zuob2009@126.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.C1 = 1.0;
par.C2 = 1.0;
par.slidingDis = 2;
par.n = 8;
par.n_2 = par.n^2;
par.N = par.n*3;
par.K=256;
par.errT = 1.0;
par.groupSz=8;
par.maxOfGroups=2000;
par.learnDic=0;
ENL=1;

for argI = 1:2:length(varargin)
    if (strcmp(varargin{argI}, 'ENL'))
        ENL = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'SD'))
        par.slidingDis = varargin{argI+1};
    end
end

if ENL<1.2
    ENL=0.8*ENL;
end

% consider the odd case
[IRow,ICol,~]=size(dataC);
if mod(IRow,2)==1
    dataC(IRow+1,:,:)=dataC(IRow,:,:);
end
if mod(ICol,2)==1
    dataC(:,ICol+1,:)=dataC(:,ICol,:);
end
[NN1,NN2,NN3] = size(dataC);

slidingDis=par.slidingDis;
n=par.n;
n_2=par.n_2;
errT=par.errT;
% Create an initial dictionary from the DCT frame
Pn=ceil(sqrt(par.K));
DCT=zeros(n,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:n-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);

for pag=1 : NN3
    [blocks(pag*n_2-n_2+1:pag*n_2,:),idx] = f_im2col(dataC(:,:,pag),[n,n],slidingDis);
end

nn=3;
dataC_M=imfilter(dataC,ones(nn,nn)/nn^2,'symmetric');
dataC_det=f_det_pol_fast(dataC_M);
[blocks_det,~] = f_im2col(dataC_det,[n,n],slidingDis);
imgSz=[NN1 NN2];           %原始图像的大小
subImgSz=floor((imgSz-[n n])/slidingDis+[1 1]);   %每个子图的大小
blkOrder=c_PatchSort( blocks_det,9,subImgSz );
clear blocks_det dataC_det;

% ENL_thr=ENL/4;
nn=5;
dataC_M=imfilter(dataC,ones(nn,nn)/nn^2,'symmetric');
for pag=1 : 3    
    dataC_Medf(:,:,pag)=medfilt2(dataC(:,:,pag),[nn nn],'symmetric');
    temp1=dataC(:,:,pag);
    temp2=dataC_Medf(:,:,pag);
    temp1(temp1./temp2>5)=0;
%     dataC_M(:,:,pag)=imfilter(dataC(:,:,pag),ones(nn,nn)/nn^2,'symmetric');
    dataC_MNew(:,:,pag)=imfilter(temp1,ones(nn,nn)/nn^2,'symmetric');
end
clear dataC_Medf

% for pag=4 : NN3
%     dataC_M(:,:,pag)=imfilter(dataC(:,:,pag),ones(nn,nn)/nn^2,'symmetric');
% end

for pag=1 : 3
    dataC_std=par.C1*dataC_MNew(:,:,pag)/sqrt(ENL);
    [blocks_std(pag*n_2-n_2+1:pag*n_2,:),~] = f_im2col(dataC_std,[n,n],slidingDis);
end
clear dataC_MNew;

C=par.C2;
pag=4;
dataC_std=C*sqrt((dataC_M(:,:,4).^2-dataC_M(:,:,7).^2+dataC_M(:,:,1).*dataC_M(:,:,2))/(2*ENL));
[blocks_std(pag*n_2-n_2+1:pag*n_2,:),~] = f_im2col(dataC_std,[n,n],slidingDis);
pag=5;
dataC_std=C*sqrt((dataC_M(:,:,5).^2-dataC_M(:,:,8).^2+dataC_M(:,:,1).*dataC_M(:,:,3))/(2*ENL));
[blocks_std(pag*n_2-n_2+1:pag*n_2,:),~] = f_im2col(dataC_std,[n,n],slidingDis);
pag=6;
dataC_std=C*sqrt((dataC_M(:,:,6).^2-dataC_M(:,:,9).^2+dataC_M(:,:,2).*dataC_M(:,:,3))/(2*ENL));
[blocks_std(pag*n_2-n_2+1:pag*n_2,:),~] = f_im2col(dataC_std,[n,n],slidingDis);
pag=7;
dataC_std=C*sqrt((-dataC_M(:,:,4).^2+dataC_M(:,:,7).^2+dataC_M(:,:,1).*dataC_M(:,:,2))/(2*ENL));
[blocks_std(pag*n_2-n_2+1:pag*n_2,:),~] = f_im2col(dataC_std,[n,n],slidingDis);
pag=8;
dataC_std=C*sqrt((-dataC_M(:,:,5).^2+dataC_M(:,:,8).^2+dataC_M(:,:,1).*dataC_M(:,:,3))/(2*ENL));
[blocks_std(pag*n_2-n_2+1:pag*n_2,:),~] = f_im2col(dataC_std,[n,n],slidingDis);
pag=9;
dataC_std=C*sqrt((-dataC_M(:,:,6).^2+dataC_M(:,:,9).^2+dataC_M(:,:,2).*dataC_M(:,:,3))/(2*ENL));
[blocks_std(pag*n_2-n_2+1:pag*n_2,:),~] = f_im2col(dataC_std,[n,n],slidingDis);
clear dataC_M;

% patch ordering
blocks = blocks(:,blkOrder);
blocks_std = blocks_std(:,blkOrder);

numOfBlks=size(blocks,2);
vecOfMeans = mean(blocks);
blocks = blocks-ones(size(blocks,1),1)*vecOfMeans;

X=zeros(size(blocks,1)/NN3,size(blocks,2)*NN3);
X_std=zeros(size(blocks,1)/NN3,size(blocks,2)*NN3);
for pag=1 : NN3
    X(:,pag:NN3:end)=blocks(pag*n_2-n_2+1:pag*n_2,:);
    X_std(:,pag:NN3:end)=blocks_std(pag*n_2-n_2+1:pag*n_2,:);
end
clear blocks_std;

blkSD=par.groupSz*NN3;
blkWidth=par.groupSz*NN3;
blkW2SD=blkWidth/blkSD;
blkTimes=ceil(numOfBlks*9/blkSD);
%==========================================================================
k=blkTimes-blkW2SD;
ind_groups=double([0:blkWidth:blkWidth*k]);
maxOfGroups=par.maxOfGroups;
count=floor(length(ind_groups)/maxOfGroups);

if count>0
    ind_groups_tmp=double([0:blkWidth:blkWidth*(maxOfGroups-1)]);
    len=maxOfGroups*blkWidth;
    for k=1 : count
        Coefs = wsomp2(DCT,X(:,len*k-len+1:len*k),...
            1./X_std(:,len*k-len+1:len*k),...
            DCT'*(X(:,len*k-len+1:len*k)./X_std(:,len*k-len+1:len*k).^2),...
            (DCT.*DCT)'*(1./X_std(:,len*k-len+1:len*k).^2),...
            ind_groups_tmp, errT);
        X(:,len*k-len+1:len*k)=DCT*Coefs;
    end
    ind_groups_tmp=double([0:blkWidth:blkWidth*(length(ind_groups)-2000*count-1)]);
    Coefs = wsomp2(DCT,X(:,len*k+1:end),...
        1./X_std(:,len*k+1:end),...
        DCT'*(X(:,len*k+1:end)./X_std(:,len*k+1:end).^2),...
        (DCT.*DCT)'*(1./X_std(:,len*k+1:end).^2),...
        ind_groups_tmp, errT);
    X(:,len*k+1:end)=DCT*Coefs;
end

if count==0
    Coefs = wsomp2(DCT,X,1./X_std,DCT'*(X./X_std.^2),(DCT.*DCT)'*(1./X_std.^2),...
        ind_groups, errT);
    X=DCT*Coefs;
end

clear X_std Coefs;

for pag=1 : NN3
    blocks(pag*n_2-n_2+1:pag*n_2,:)=X(:,pag:NN3:end);    
end
clear X;

blocks = blocks + ones(size(blocks,1),1)*vecOfMeans;
blocks(:,blkOrder(1 : numOfBlks))=blocks(:,1 : numOfBlks);
blocks(1:3*n_2,:)=abs(blocks(1:3*n_2,:));
%% subimage averaging

count = 1;
Weight= zeros(NN1,NN2,NN3);
IMout = zeros(NN1,NN2,NN3);
[rows,cols] = ind2sub(size(dataC(:,:,1))-n+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);
    for pag=1 : NN3
        block =reshape(blocks(pag*n_2-n_2+1:pag*n_2,count),[n,n]);
        IMout(row:row+n-1,col:col+n-1,pag)=IMout(row:row+n-1,col:col+n-1,pag)+block;
        Weight(row:row+n-1,col:col+n-1,pag)=Weight(row:row+n-1,col:col+n-1,pag)+ones(n);
    end    
    count = count+1;
end;
IOut = IMout./Weight;
IOut=IOut(1:IRow,1:ICol,:);
end

