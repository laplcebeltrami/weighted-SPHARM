%% Weighted Spherical Harmonic Representation (SPHARM)
%
% The representation and the online estimation procedure is published in
% the following three papers. If you are using codes below, please
% reference one of the following papers:
%
% [1] Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. 
% Weighted Fourier series representation and its application to quantifying the amount 
% of gray matter. Special Issue of  IEEE Transactions on Medical Imaging, on 
% Computational Neuroanatomy. 26:566-581 
% https://pages.stat.wisc.edu/~mchung/papers/TMI.SPHARM.2007.pdf
%
% [2] Chung, M.K., Dalton, K.M., Davidson, R.J. 2008. . Tensor-based cortical surface 
% morphometry via weighed spherical harmonic representation. IEEE Transactions on 
% Medical Imaging 27:1143-1151 
% https://pages.stat.wisc.edu/~mchung/papers/TMI.2008.pdf
%
% [3] Chung, M.K. Hartley, R., Dalton, K.M., Davidson, R.J. 2008. Encoding cortical 
% surface by spherical harmonics.  Satistica Sinica 18:1269-1291
%https://pages.stat.wisc.edu/~mchung/papers/sinica.2008.pdf
%
% (C) 2022 Moo K. Chung 
% University of Wisconsin-Madison
% mkchung@wisc.edu
%
% The codes are downloaded from https://github.com/laplcebeltrami/weighted-SPHARM
% The code tested with MATLAB 2019b. SPHARM basis may not work with
% older versions of MATLAB. 
%
%% Autism brain surface dataset
% The dataset is published in Chung, M.K., Robbins,S., Dalton, K.M., Davidson, 
% Alexander, A.L., R.J., Evans, A.C. 2005. Cortical thickness analysis in autism 
% via heat kernel smoothing. NeuroImage 25:1256-1265.
%
% Load full data and then display the inner and outer surface of the 1st autistic subject

load autism.surface.mat

surf.vertices=squeeze(autisminner(1,:,:))';
surf.faces=tri;
figure; subplot(1,2,1); figure_wire(surf, 'yellow', 'white')

surf.vertices=squeeze(autismouter(1,:,:))';
surf.faces=tri;
subplot(1,2,2); figure_wire(surf, 'yellow', 'white')
%% Spherical mesh with 40962 triangles
% Project the spherical angles $\theta$ and $\varphi$ onto the brain surface

% Load spherical mesh
load sphere40962.mat

% Compute two spherical angles
[theta varphi]=SPHARMangles(sphere40962);

figure; subplot(1,2,1); figure_trimesh(sphere40962,theta,'rwb');
subplot(1,2,2); figure_trimesh(surf,theta,'rwb');

figure; subplot(1,2,1); figure_trimesh(sphere40962, varphi,'rwb');
subplot(1,2,2); figure_trimesh(surf,varphi,'rwb');

%% Coordinate functions on sphere
% Map the x-,y-,z-coordinates of the brain surface onto the spheree

x=surf.vertices(:,1);
y=surf.vertices(:,2);
z=surf.vertices(:,3);

figure; subplot(1,3,1); figure_trimesh(sphere40962, x,'rwb');
subplot(1,3,2); figure_trimesh(sphere40962, y,'rwb');
subplot(1,3,3); figure_trimesh(sphere40962, z,'rwb');

%% Spherical harmonics (SPHARM)
% Compute spherical harmonic basis $Y_{31}(\theta, \varphi)$ and $Y_{30,-20}(\theta, \varphi)$ on the sphere

Ylm=Y_lm(3,1,theta',varphi'); 
figure; subplot(1,2,1); figure_trimesh(sphere40962,Ylm,'rwb');
Ylm=Y_lm(30,-20,theta',varphi'); 
subplot(1,2,2); figure_trimesh(sphere40962,Ylm,'rwb');


%% Spherical harmonic (SPHARM) expansion 
% SPHARM expands the x-,y-,z-coordinates using spherical harmonic basis

%Degree 5 expansion
[surfsmooth, fourier]=SPHARMsmooth3(surf,sphere40962,20,0);
figure; subplot(1,3,1); figure_wire(surfsmooth,'yellow','white');

%Degree 50 expansion
[surfsmooth, fourier]=SPHARMsmooth3(surf,sphere40962,50,0);
subplot(1,3,2); figure_wire(surfsmooth,'yellow','white');

%Degree 100 expansion (computation is slow)
[surfsmooth, fourier]=SPHARMsmooth3(surf,sphere40962,100,0);
subplot(1,3,3); figure_wire(surfsmooth,'yellow','white');

