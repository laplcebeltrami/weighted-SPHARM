function Y_lm=Y_lm(l,m,theta,varphi)
%Y_lm=Y_lm(l,m,theta,varphi)
%
% Computes the spherical harmonic function based on the method published in
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
% Update history 2005; July 2007; Mar 2008; 

Y=(-1)^m*legendre(l,cos(theta));	
Y=Y(abs(m)+1,:);
Ylm=Y;

if m>0 
	Ylm=Ylm.*cos(abs(m)*varphi);
end;
if m<0
	Ylm=Ylm.*sin(abs(m)*varphi);
end;
if m==0
    Ylm=Ylm.*cos(abs(m)*varphi)/sqrt(2);
end;

% old line 
%clm=sqrt((2*l+1)*factorial(l-abs(m))/(2*pi*factorial(l+abs(m))));
% factorial routine is very time consuming

for k = 0:(2*l)
    fact(k+1) = factorial(k);
end
clm = sqrt(((2*l+1)/(2*pi))*(fact(l-abs(m)+1)./fact(l+abs(m)+1)));
Y_lm=clm*Ylm;
	