function [theta,varphi]=EULERangles(sphere)
%function [theta,varphi]=EULERangles(sphere)
%
%Given a spherical mesh, it generates Euler angles \theta and \varphi used
%in Chung et al. (2007), which follows the most widely used mathematical 
%convention and it is different from MATLAB convenstion.
%
% sphere: spherical mesh
% theta, varphi: spherical angles
%
%
% References 
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
% (C) 2006- Moo K. Chung 
%  mkchung@wisc.edu
%  Department of Biostatisics and Medical Informatics
%  University of Wisconsin, Madison
%
% The codes are downloaded from https://github.com/laplcebeltrami/weighted-SPHARM
% The code tested with MATLAB 2019b. SPHARM basis may not work with
% older versions of MATLAB. 

n_vertex=size(sphere.vertices,1);
c=mean(sphere.vertices);  %mass center
sphere.vertices=sphere.vertices-kron(ones(n_vertex,1),c);  % translation

[theta,varphi,r] = cart2sph(sphere.vertices(:,1),sphere.vertices(:,2),sphere.vertices(:,3));

% MATLAB coordinate systems are different from the convention used in the
% TMI paper.
temp = theta;
theta = pi/2 - varphi;
varphi = pi + temp;

%figure_wire(sphere,'yellow')

