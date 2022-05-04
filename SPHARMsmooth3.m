function [surf_smooth, fourier]=SPHARMsmooth3(surf,sphere, L,sigma)
%[surf_smooth, fourier]=SPHARMsmooth2(surf,sphere, L,sigma)
%
% INPUT: You need to read references [1],[2] or [3] to understand the input prameters.
%
% surf            : A surface mesh format identical to the MATLAB format
%                
%  sphere         : Spherical mesh where the representation is constructed
%
%    L            : The maximal degree of the weighted-SPHARM representation.
%                      Read the paper below to find it optimally.
%
%  sigma          : bandwith of weighted-SPHARM representation
%                    It is the bandwidth of heat kernel on unit sphere.
%                    When sigma=0, it is the traditional SPHARM representation.
%                    range beween 0.0001 and 0.01 will be sufficient for cortical surfaces.
%
% surf_smooth     : The weighted-SPHARM result.
%
% fourier         : The estimated SPHARM coefficients (Fourier coefficients) given as a structured array
%                   containg coeff.x, coeff.y, coeff.z
%                   coeff.x is the SPHARM coefficients of x-cooridinate given as (L+1) by (2*L+1) matrix
%                   coeff.x(3,:) is the 2nd degree SPHARM coefficients of all order.
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
%
% (C) 2006- Moo K. Chung
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%
% email://mkchung@wisc.edu
%
% Update history Sept 19 2006; July 5, 2007; January 22, 2008; July 18, 2010
% WARNING: This version is different from SPHARMsmooth.m in
% http://www.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/SPHARMsmooth.m
% SPHARMsmooth.m is based on a fixed connectivity across all surface. If
% every subjects have different mesh topology, use SPHARM.smooth2.m.

coord=surf.vertices;
n_vertex = size(coord,1);   % the number of vertices in a mesh.

[theta varphi]=EULERangles(sphere);

x=coord(:,1);
y=coord(:,2);
z=coord(:,3);

%INITIALIZATION

% xestiamte is the weighted-SPHARM of x-coordinate.
xestimate=zeros(n_vertex,1);
yestimate=zeros(n_vertex,1);
zestimate=zeros(n_vertex,1);

% betax is the Fourier coefficients of x-coordinate.
betax=zeros(L+1,2*L+1);
betay=zeros(L+1,2*L+1);
betaz=zeros(L+1,2*L+1);


%0-TH DEGREE. 
%Step 2 in the iterative resiual fitting (IRF) algorithm. See reference [1].

Y=Y_l(0,theta',varphi')';
Ycommon=inv(Y'*Y)*Y';

betal=Ycommon*x;
betax(1,1)=betal';
xsmooth=Y*betal;
xestimate=xsmooth;

betal=Ycommon*y;
betay(1,1)=betal';
ysmooth=Y*betal;
yestimate=ysmooth;

betal=Ycommon*z;
betaz(1,1)=betal';
zsmooth=Y*betal;
zestimate=zsmooth;


%l-TH DEGREE ITERATION


for l=1:L
 
    % Step 4: residual. See the paper for detail
    x_j = x-xestimate;
    y_j = y-yestimate;
    z_j = z-zestimate;
 
    Y=Y_l(l,theta',varphi')';
    % real(Y_l) gives negative harmonics imag(Y_l) gives positive harmonics
   
    Y=[real(Y) imag(Y(:,2:(l+1)))];
    Y=real(Y);
    Ycommon=inv(Y'*Y)*Y';  %changed from Ycommon=inv(Y'*Y)*Y';

   %Y(1,:)
   
    % Step 5: refitting the residual. See the paper for detail
    betal=Ycommon*x_j;
    betax(l+1,1:2*l+1)=betal';
    xsmooth=Y*betal;
    xestimate=xestimate + exp(-l*(l+1)*sigma)*xsmooth;

    betal=Ycommon*y_j;
    betay(l+1,1:2*l+1)=betal';
    ysmooth=Y*betal;
    yestimate=yestimate + exp(-l*(l+1)*sigma)*ysmooth;

    betal=Ycommon*z_j;
    betaz(l+1,1:2*l+1)=betal';
    zsmooth=Y*betal;
    zestimate=zestimate + exp(-l*(l+1)*sigma)*zsmooth;
    
% animation. remove next 4 lines if you don't want animation
%     temp=[xestimate; yestimate; zestimate];
%     surf_smooth.vertices=squeeze(reshape(temp,n_vertex,3));
%     surf_smooth.faces=surf.faces;
%     M(l) = getframe; 
%     figure_wire(surf_smooth,'yellow','white');

    
end;


%output the results in a proper shape
temp=[xestimate; yestimate; zestimate];
surf_smooth.vertices=squeeze(reshape(temp,n_vertex,3));
surf_smooth.faces=surf.faces;

fourier.x=betax;
fourier.y=betay;
fourier.z=betaz;

%---------------------------------------------------------------
function [theta,varphi]=EULERangles(surf);

n_vertex=size(surf.vertices,1);
c=mean(surf.vertices);  %mass center
surf.vertices=surf.vertices-kron(ones(n_vertex,1),c);  % translation

[theta,varphi,r] = cart2sph(surf.vertices(:,1),surf.vertices(:,2),surf.vertices(:,3));

% MATLAB coordinate systems are different from the convention used in the
% TMI paper.
temp = theta;
theta = pi/2 - varphi;
varphi = pi + temp;

%figure_wire(surf,'yellow')

%-----------------------------------------------------------------
function Y_l=Y_l(l,theta,varphi)
% computes spherical harmonics of degree l.
sz=length(theta);

m=0:l;
CLM=[];
exp_i_m=[];
sign_m=[];
SIGNM=[];
Pn=[];

for k = 0:(2*l)
    fact(k+1) = factorial(k);
end
%clm = sqrt(((2*l+1)/(2*pi))*(fact(l-abs(m)+1)./fact(l+abs(m)+1))); % this clm will be close to 0 when l is large.
clm=(-1).^m*sqrt(1/pi); % clm needs to be modified because Pn is normalized
CLM=kron(ones(1,sz),clm');

for k = 0:l
    exp_i_m(k+1,:)= exp(i*k*varphi);
    sign_m(k+1) = (-1)^k;
end

c=sqrt(2);
exp_i_m(1,:)=exp_i_m(1,:)/sqrt(2);

SIGNM=kron(ones(1,sz),sign_m');
%Pn=legendre(l,cos(theta));     %this Pn goes to Inf when l is large
Pn=legendre(l,cos(theta),'norm'); % Pn is normalized to avoid Inf, and thus clm needs to be modified
Y_l=CLM.*SIGNM.*Pn.*exp_i_m;
