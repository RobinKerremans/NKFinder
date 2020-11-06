%	Copyright 2019 Robin Kerremans, Paul Meredith, Ardalan Armin, Swansea University
%   
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%__________________________________________________________________________ 
%   This Transfer Matrix script is capable of Transmission and Reflection
%   calculation of any stack under any angle of incidence. 
%   It also has the capability to simulate incoherent layers (arbitrarily taken as any layer above 10.000nm thick) 
%   in any position within the mulitlayer stack.
% 
%   Inputs are the wavelength, optical constant vector at said wavelength, 
%   thickness vector of the stack and the incident angle of light.
% 
%	The mathematical procedure for calculating the optical field via the transfer matrix method 
%	was taken from:
%   Harbecke B., "Coherent and incoherent reflection and transmission of multilayer structures." Applied Physics B 1986, 39:165-170.

%% Main function code
% Get full Reflection and Transmission from the polarized parts.
function [R,T] = Tmat(lambda,n,L,theta_in)
[Rs,Ts] = Tmats(lambda,n,L,theta_in);
[Rp,Tp] = Tmatp(lambda,n,L,theta_in);
R=Rs/2+Rp/2;
T=Ts/2+Tp/2;
end

%% Tmat function code (s-polarized)

% Combines Tstack with optional incoherent layers
% For proof of this full matrix, see Appl. Phys. B 39,165-170 (1986) 
% "Coherent and Incoherent Reflection and Transmission of Multilayer Structures"

function [R,T,rcd,tcd] = Tmats(lambda,n,L,theta_in)
% layers are assumed (somewhat arbitrarily) to become incoherent at +0.1mm
Incohtest1=logical(L(1)>=100000);
Incohtest2=logical(L(length(L))>=100000);

% theta calculation (Snell's law)
theta(1)=theta_in;
for i=2:length(n)
    theta(i)=asin(n(i-1)/n(i)*sin(theta(i-1)));
end

if(Incohtest1) %incoherent layer at the front
    phi1=exp(-2*pi*imag(n(2))/lambda*L(1)*cos(theta(2)));
    rAG1=(n(1)*cos(theta(1))-n(2)*cos(theta(2)))/(n(1)*cos(theta(1))+n(2)*cos(theta(2)));
    rGA1=-rAG1;
    tAG1=2*n(1)*cos(theta(1))/(n(1)*cos(theta(1))+n(2)*cos(theta(2)));
    tGA1=2*n(2)*cos(theta(2))/(n(1)*cos(theta(1))+n(2)*cos(theta(2)));
    mat1=[1 -abs(rAG1)^2;abs(rAG1)^2 abs(tAG1*tGA1)^2-abs(rAG1*rGA1)^2]*[1 0;0 abs(phi1^2)^2];
    %trim out the variables for incoherent layer 
    n(1)=[];
    L(1)=[];
    theta(1)=[];
end

if(Incohtest2) %incoherent layer at the back 
    phi2=exp(-2*pi*imag(n(length(n)-1))/lambda*L(length(L))*cos(theta(length(n)-1)));
    rAG2=(n(length(n)-1)*cos(theta(length(n)-1))-n(length(n))*cos(theta(length(n))))/(n(length(n)-1)*cos(theta(length(n)-1))+n(length(n))*cos(theta(length(n))));
    rGA2=-rAG2;
    tGA2=2*n(length(n)-1)*cos(theta(length(n)-1))/(n(length(n)-1)*cos(theta(length(n)-1))+n(length(n))*cos(theta(length(n))));
    tAG2=2*n(length(n))*cos(theta(length(n)))/(n(length(n)-1)*cos(theta(length(n)-1))+n(length(n))*cos(theta(length(n))));
    mat2=[1 0;0 abs(phi2^2)^2]*[1 -abs(rAG2)^2;abs(rAG2)^2 abs(tAG2*tGA2)^2-abs(rAG2*rGA2)^2];
    theta(length(n))=[];
    n(length(n))=[];
    L(length(L))=[];  
end

[rcd,tcd]=Tstacks(lambda,n,L,theta);

%reverse layers
nn=n;
tt=theta;
for i=1:length(n)                                   
    n(i)=nn(length(n)+1-i);
    theta(i)=tt(length(n)+1-i);
end
LL=L;
for i=1:length(L)
    L(length(L)+1-i)=LL(i);
end

[rdc,tdc]=Tstacks(lambda,n,L,theta);

% Full Matrix construction
matmain=[1 -abs(rdc)^2;abs(rcd)^2 abs(tcd*tdc)^2-abs(rcd*rdc)^2];
if(Incohtest1&&Incohtest2)
    matfull=mat1*matmain*mat2;
    z=matfull(1,1);
    R=matfull(2,1)/z;
    T=abs(tAG1*phi1*phi2*tcd*tGA2)^2/z;
elseif(Incohtest2)
    matfull=matmain*mat2;
    z=matfull(1,1);
    R=matfull(2,1)/z;
    T=abs(phi2*tcd*tGA2)^2/z;
elseif(Incohtest1)
    matfull=mat1*matmain;
    z=matfull(1,1);
    R=matfull(2,1)/z;
    T=abs(phi1*tcd*tAG1)^2/z;
else
    matfull=matmain;
    z=matfull(1,1);
    R=matfull(2,1)/z;
    T=abs(tcd)^2/z;
end
end

%% Tmat function code = p-polarized

function [R,T] = Tmatp(lambda,n,L,theta_in)
% layers are assumed (somewhat arbitrarily) to become incoherent at +0.1mm
Incohtest1=logical(L(1)>=100000);
Incohtest2=logical(L(length(L))>=100000);

% theta calculation (Snell's law)
theta(1)=theta_in;
for i=2:length(n)
    theta(i)=asin(n(i-1)/n(i)*sin(theta(i-1)));
end

if(Incohtest1) %incoherent layer at the front
    phi1=exp(-2*pi*imag(n(2))/lambda*L(1)*cos(theta(2)));
    rAG1=(n(1)*cos(theta(2))-n(2)*cos(theta(1)))/(n(1)*cos(theta(2))+n(2)*cos(theta(1)));
    rGA1=-rAG1;
    tAG1=2*n(1)*cos(theta(1))/(n(1)*cos(theta(2))+n(2)*cos(theta(1)));
    tGA1=2*n(2)*cos(theta(2))/(n(1)*cos(theta(2))+n(2)*cos(theta(1)));
    mat1=[1 -abs(rAG1)^2;abs(rAG1)^2 abs(tAG1*tGA1)^2-abs(rAG1*rGA1)^2]*[1 0;0 abs(phi1^2)^2];
    %trim out the variables for incoherent layer 
    n(1)=[];
    L(1)=[];
    theta(1)=[];
end

if(Incohtest2) %incoherent layer at the back 
    phi2=exp(-2*pi*imag(n(length(n)-1))/lambda*L(length(L))*cos(theta(length(n)-1)));
    rAG2=(n(length(n)-1)*cos(theta(length(n)))-n(length(n))*cos(theta(length(n)-1)))/(n(length(n)-1)*cos(theta(length(n)))+n(length(n))*cos(theta(length(n)-1)));
    rGA2=-rAG2;
    tGA2=2*n(length(n)-1)*cos(theta(length(n)-1))/(n(length(n)-1)*cos(theta(length(n)))+n(length(n))*cos(theta(length(n)-1)));
    tAG2=2*n(length(n))*cos(theta(length(n)))/(n(length(n)-1)*cos(theta(length(n)))+n(length(n))*cos(theta(length(n)-1)));
    mat2=[1 0;0 abs(phi2^2)^2]*[1 -abs(rAG2)^2;abs(rAG2)^2 abs(tAG2*tGA2)^2-abs(rAG2*rGA2)^2];
    theta(length(n))=[];
    n(length(n))=[];
    L(length(L))=[];  
end

[rcd,tcd]=Tstackp(lambda,n,L,theta);

%reverse layers
nn=n;
tt=theta;
for i=1:length(n)                                   
    n(i)=nn(length(n)+1-i);
    theta(i)=tt(length(n)+1-i);
end
LL=L;
for i=1:length(L)
    L(length(L)+1-i)=LL(i);
end

[rdc,tdc]=Tstackp(lambda,n,L,theta);

% Full Matrix construction
matmain=[1 -abs(rdc)^2;abs(rcd)^2 abs(tcd*tdc)^2-abs(rcd*rdc)^2];
if(Incohtest1&&Incohtest2)
    matfull=mat1*matmain*mat2;
    z=matfull(1,1);
    R=matfull(2,1)/z;
    T=abs(tAG1*phi1*phi2*tcd*tGA2)^2/z;
elseif(Incohtest2)
    matfull=matmain*mat2;
    z=matfull(1,1);
    R=matfull(2,1)/z;
    T=abs(phi2*tcd*tGA2)^2/z;
elseif(Incohtest1)
    matfull=mat1*matmain;
    z=matfull(1,1);
    R=matfull(2,1)/z;
    T=abs(phi1*tcd*tAG1)^2/z;
else
    matfull=matmain;
    z=matfull(1,1);
    R=matfull(2,1)/z;
    T=abs(tcd)^2/z;
end
end

%% Tstack function code (s-polarized)

% Basic == coherent == Transformation Matrix function where the input is a wavelength and two vectors
% (Refractive indices and lengths of the respective layers)
% Output is reflection and transmission of ELECTRIC FIELD (not Intensity)
% APPLIED OPTICS / Vol. 29, No. 13 / 1 May 1990
% "Matrix formalism for calculation of electric field intensity of light in stratified multilayered films"

function [rE,tE] = Tstacks(lambda,n,L,theta) 
lay=length(L);                              % Lay is the number of layers
for i=1:lay+2                               % Convert n to k vectors, note that surrounding medium indices are included as first and last
    k(i)=2*pi*n(i)/lambda*cos(theta(i));
end

Mtot=eye(2);                                % Mtot will be the transfer matrix of the entire stack

for i=1:lay+1                               % Calculate the total transfer matrix to be returned, Abeles formalism, theta_incident=0
    if(i==1)
        d(1)=0;
    else
        d(i)=k(i)*L(i-1);                   % Pay attention to indices, k(i) and n(i) have extra entries for surrounding layers   
    end
    r(i)=(n(i)*cos(theta(i))-n(i+1)*cos(theta(i+1)))/(n(i)*cos(theta(i))+n(i+1)*cos(theta(i+1)));
    t(i)=2*n(i)*cos(theta(i))/(n(i)*cos(theta(i))+n(i+1)*cos(theta(i+1)));
    M{i}=Mcalc(r(i),d(i));            
    Mtot= (Mtot*M{i});                      % note that Mtot does not include the t(i) yet.
    if(i==1)
        tTOTa=t(i);
    else
        tTOTa=tTOTa*t(i);
    end
end

%Calculate tE and rE
tE=tTOTa/Mtot(1,1);
rE=Mtot(2,1)/Mtot(1,1);                            
end

%% Tstack function code = p-polarized

function [rE,tE] = Tstackp(lambda,n,L,theta) 
lay=length(L);                              % Lay is the number of layers
for i=1:lay+2                               % Convert n to k vectors, note that surrounding medium indices are included as first and last
    k(i)=2*pi*n(i)/lambda*cos(theta(i));
end

Mtot=eye(2);                                % Mtot will be the transfer matrix of the entire stack

for i=1:lay+1                               % Calculate the total transfer matrix to be returned, Abeles formalism, theta_incident=0
    if(i==1)
        d(1)=0;
    else
        d(i)=k(i)*L(i-1);                   % Pay attention to indices, k(i) and n(i) have extra entries for surrounding layers   
    end
    r(i)=(n(i)*cos(theta(i+1))-n(i+1)*cos(theta(i)))/(n(i)*cos(theta(i+1))+n(i+1)*cos(theta(i)));
    t(i)=2*n(i)*cos(theta(i))/(n(i)*cos(theta(i+1))+n(i+1)*cos(theta(i)));
    M{i}=Mcalc(r(i),d(i));            
    Mtot= (Mtot*M{i});                      % note that Mtot does not include the t(i) yet.
    if(i==1)
        tTOTa=t(i);
    else
        tTOTa=tTOTa*t(i);
    end
end

%Calculate tE and rE
tE=tTOTa/Mtot(1,1);
rE=Mtot(2,1)/Mtot(1,1);                            
end

%% Transfer matrix function
%calculate one transfer matrix (1 layer)
%Abeles formalism is used
function M = Mcalc(r,d)                                                   
M=ones(2);
M(1,1)=exp(-j*d);
M(1,2)=r*exp(-j*d);
M(2,1)=r*exp(j*d);
M(2,2)=exp(j*d);
end