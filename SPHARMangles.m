function [theta,varphi]=SPHARMangles(surf);
% function [theta,varphi]=SPHARMangles(surf);
%
% Given a mesh, it generates spherical angles \theta and \varphi with
% respect to the mass center. The spehrical angles are introduced in
% in Chung et al. (2007), which follows the most widely used mathematical 
% convention and it is different from MATLAB convention.
%
%  INPUT
%  surf: surface mesh
%  theta, varphi: spherical angles
%
% OUTPUT
% angles \theta, \varphi
%
% (C) 2006- Moo K. Chung 
%  mkchung@wisc.edu
%  Department of Biostatisics and Medical Informatics
%  University of Wisconsin, Madison
%
% The SPAHRM representation and the online estimation procedure is published in
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
% The codes are downloaded from
% https://github.com/laplcebeltrami/weighted-SPHARM or
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
% The code tested with MATLAB 2019b. SPHARM basis may not work with
% older versions of MATLAB. 


n_vertex=size(surf.vertices,1);
%c=mean(surf.vertices);  %mass center
%surf.vertices=surf.vertices-kron(ones(n_vertex,1),c);  % translation

[theta,varphi,r] = cart2sph(surf.vertices(:,1),surf.vertices(:,2),surf.vertices(:,3));

% MATLAB coordinate systems are different from the convention used in the
% TMI paper.
temp = theta;
theta = pi/2 - varphi;
varphi = pi + temp;

%figure_wire(surf,'yellow')

