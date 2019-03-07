
function MAR(threshold,method,dicom_directory,opt)

% threshod: used to segment metallic implants

% dicom_directory:  directory of dicom images

% method: [1] 'NMAR',  [2] 'iMAR'
% [1] Meyer, E., et al. (2010). "Normalized metal artifact reduction (NMAR) in computed tomography." Med Phys 37(10): 5482-5493.
% [2] Mehranian, A., et al. (2013). "3D Prior Image Constrained Projection Completion for X-ray CT Metal Artifact Reduction." IEEE Transactions on Nuclear Science 60(5): 3318-3332.

% opt: options for iMAR method
% opt.nIter: number of iterations
% opt.alpha: controls the impact of prior sinogram

if nargin==0
    method = 'iMAR';
    threshold = 3500; %
    dicom_directory = uigetdir();
end

if nargin==2
    dicom_directory = uigetdir();
end

if nargin<4
    opt.nIter = 200;
    opt.alpha = 1;
end


% read Dicom files
[imgPatientOriginal,I] = dicom23D(dicom_directory);

% find affected slices
k = [];
for i = 1:size(imgPatientOriginal,3)
    isAffected = imgPatientOriginal(:,:,i)> threshold;
    if any(isAffected(:))
        k = [k,i];
    end
end

if isempty(k)
    error('no slice was found affected, change the threshod\n');
end

% segment and dilate the metal implants
imgPatient = imgPatientOriginal(:,:,k);
imgMetals = imgPatient> threshold;

if 0
    [x,y,z] = meshgrid(-1:1:1);
    r = (x.^2+y.^2+z.^2)<2;
    
    imgMetals = imdilate(single(imgMetals),r)>0;
end

% generate metal sinograms, do linear interpolation MAR
A = @(x) fanbeam(x,600,'FanSensorSpacing',0.1,'FanSensorGeometry','arc','FanRotationIncrement',0.1);
At = @(x) ifanbeam(x,600,'FanSensorSpacing',0.1,'FanSensorGeometry','arc','OutputSize',512);

sinogramSize = size(A(ones(size(imgPatient(:,:,1)))));
sinoMetals = zeros([sinogramSize, length(k)],'logical');
sinoPatient = zeros([sinogramSize, length(k)],'single');
sinoPatientLiMar = sinoPatient;
imgPatientLiMAR = 0*imgPatient;

h = waitbar(0,'Linear Interpolation...','WindowStyle','modal');
N = length(k);
for i = 1: N
    sinoMetals(:,:,i) = A(imgMetals(:,:,i)) > 0;
    sinoPatient(:,:,i) = single(A(imgPatient(:,:,i)));
    sinoPatientLiMar(:,:,i) = linear_iterpolation (sinoPatient(:,:,i),sinoMetals(:,:,i));
    imgPatientLiMAR(:,:,i) = At(sinoPatientLiMar(:,:,i));
    waitbar(i/N,h);
end

close(h);
imgPatientLiMAR = max(0,imgPatientLiMAR);

% generate the prior image and sinograms
soft =  (imgPatientLiMAR -1024) >-501;
se = strel('disk',1);
soft = imerode(imdilate(soft,se),se);

metal_di = imdilate(imgMetals,strel('disk',3));

bone=(imgPatientLiMAR.*~metal_di -1024)>350;

imgPrior = 1100*(soft).*~bone + min(2500,imgPatientLiMAR.*bone);

imgPrior = imfilter(imgPrior,fspecial('gaussian',[5,5],0.8));
% %
if strcmpi(method,'NMAR')
    imgPatientNMAR = 0*imgPatientLiMAR;
    
    for i = 1: length(k)
        sinoPrior = A(imgPrior(:,:,i));
        sinoPatientNorm = sinoPatient(:,:,i)./sinoPrior;
        sinoPatientNorm(isinf(sinoPatientNorm))=1;
        sinoPatientNorm(isnan(sinoPatientNorm))=0;
        
        sinoPatientNorm = linear_iterpolation (sinoPatientNorm,sinoMetals(:,:,i));
        
        sinoPatientNMAR = sinoPatientNorm .* sinoPrior;
        
        sinoPatientNMAR(isnan(sinoPatientNMAR))= 0 ;
        sinoPatientNMAR = sinoPatient(:,:,i).*~sinoMetals(:,:,i) + sinoPatientNMAR.*sinoMetals(:,:,i);
        
        imgPatientNMAR(:,:,i) = At(sinoPatientNMAR);
    end
    imgPatientMAR = imgPatientNMAR;
    
elseif strcmpi(method,'iMAR')
    h = waitbar(0,'Iterative MAR...','WindowStyle','modal');
    imgPatientIMAR = 0*imgPatientLiMAR;
    N = length(k);
    for i = 1: N
        sinoPrior = A(imgPrior(:,:,i));
        
        sinoPatientIMAR =  iMAR(sinoPatient(:,:,i), sinoMetals(:,:,i),sinoPatient(:,:,i),sinoPrior,opt);
        
        sinoPatientIMAR(isnan(sinoPatientIMAR))= 0 ;
        sinoPatientIMAR = sinoPatient(:,:,i).*~sinoMetals(:,:,i) + sinoPatientIMAR.*sinoMetals(:,:,i);
        
        imgPatientIMAR(:,:,i) = At(sinoPatientIMAR);
        waitbar(i/N,h);
    end
    close(h);
    imgPatientMAR = imgPatientIMAR;
    
else
    error('Unknown MAR method')
end

% %
imgPatientCorrected = imgPatientOriginal;
imgPatientCorrected(:,:,k) =imgPatientMAR.*~imgMetals+ imgMetals.*imgPatientOriginal(:,:,k);% min(,3071+1024);

% re-order the slices
n = size(imgPatientCorrected,3);
temp = 0*imgPatientCorrected;
for i = 1:n
    temp(:,:,I(i))= (imgPatientCorrected(:,:,i));
    
end

imgPatientCorrected= temp;

% write it back into dicom
files=dir([dicom_directory '\*.ima']);
if strcmpi(method,'iMAR')
    dirMar = [dicom_directory '_iMAR'];
elseif strcmpi(method,'NMAR')
    dirMar = [dicom_directory '_NMAR'];
end
mkdir(dirMar);
h = waitbar(0,'Writing DICOM MAR images...','WindowStyle','modal');
for ii=1:n
    info=dicominfo([dicom_directory,'\',files(ii).name]);
    dicomwrite(int16(imgPatientCorrected(:,:,ii)),[dirMar '\mar_' files(ii).name],info);
    waitbar(ii/n,h);
end
close(h);


fprintf('Done\n');


function sii = linear_iterpolation (sii,omega)
sii(omega)=1;
[~,c]=size(sii);
for i=1:c
    d=[];
    a=find(omega(:,i)==1);
    if isempty(a)
        continue
    end
    for ii =1:length(a)-1
        if a(ii+1)~=a(ii)+1
            d=[d,a(ii)+1,a(ii+1)-1];
        end
    end
    s=[a(1)-1,d,a(end)+1];
    for j=1:length(s)/2
        p=s(2*j-1):s(2*j);
        sii(p,i)=linspace(sii(s(2*j-1),i), sii(s(2*j),i),length(p));
    end
end

function sinoCorrected = iMAR(sinoPatient, sinoMetals,sinoPatientLiMar,sinoPrior,opt)

% trim the sinograms for faster computation
[UpperLimit, LowerLimit] = findUpperLowerLimits(sinoMetals);

sinoMetals = sinoMetals(UpperLimit:LowerLimit,:);
sinoPrior = sinoPrior(UpperLimit:LowerLimit,:);
sinoCorrected = sinoPatient(UpperLimit:LowerLimit,:);
tmp = sinoPatientLiMar(UpperLimit:LowerLimit,:);

tau = 0.8/(4);
t=1;
for i = 1:opt.nIter
    
    t0 = t;
    sinoCorrectedOld = sinoCorrected;
    sinoCorrected = tmp +  tau* sinoMetals.*div2(grad2((sinoCorrected-opt.alpha*sinoPrior)));
    t = (1+sqrt(1+4*t0^2))/2;
    tmp = sinoCorrected+(t0-1)/t*(sinoCorrected - sinoCorrectedOld);
    
end
tmp = sinoCorrected;
sinoCorrected = sinoPatient;
sinoCorrected(UpperLimit:LowerLimit,:) = tmp;


function d = div2(G)

% divergence operator, backward differences with symmetric boundry
Gx = G(:,:,1);
Gy = G(:,:,2);

fx = Gx - Gx([1 1:end-1],:);
fx(1,:) = Gx(1,:);        % boundary
fx(end,:) = -Gx(end-1,:);
clear Gx
fy = Gy - Gy(:,[1 1:end-1]);
fy(:,1) = Gy(:,1);    % boundary
fy(:,end) = -Gy(:,end-1);
clear Gy
d = fx + fy;

function G = grad2(M)
% gradient, backward differences with symmetric boundry


Gx = M([2:end end],:)-M;
Gy = M(:,[2:end end])-M;
clear M
[x,y]=size(Gx);
G = single(zeros(x,y,2));
G(:,:,1) = Gx;
G(:,:,2) = Gy;




function [UpperLimit, LowerLimit] = findUpperLowerLimits (y)
% to find the upper and lower radial positions of metal traces for sinogram
% trimming

Limits = zeros(size(y,2),2);
for i = 1:size(y,2)
    
    r = find(y(:,i)==1);
    if ~isempty(r)
        Limits(i,:) = [r(1),r(end)];
    end
end

UpperLimit = max(1,floor(min(Limits(:,1)) - min(Limits(:,1))*0.15));
LowerLimit = min(size(y,1), floor(max(Limits(:,2)) + max(Limits(:,2))*0.15));

