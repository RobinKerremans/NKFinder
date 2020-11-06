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
%   This Kramers Kronig script is capable of calculating the real optical
%   constant n of a material, starting from the imaginary
%   component k as a spectral function of wavelength, in nm (not frequency!).
% 
%   Note that k is expanded spectrally past the fed data by assuming 3
%   deeper bands underneat the valence band, using the Tauc model to derive the
%   absorption coefficient (and hence k) from the bandgap at lower wavelengths towards 0.
% 
%   Inputs are n_cauchy (at an arbitrary wavelength in the cauchy regime, 
%   the spectral k data, and the stepsize in wavelength of the k data
%  
%	For the mathematical derivation of the Kramer Kronig equations for optical
%	constants, and for discussion of the Tauc equations, see for example:
%   "Bohren CF. What did Kramers and Kronig do and how did they do it?
%   European Journal of Physics 2010, 31:573."
%   "Viezbicke BD, Patel S, Davis BE, Birnie III DP, physica status solidi
%   (b) 2015, 252:1700-1710."



function n_out = KK(n_cauchy,k_data,spacing_KK,B)

global step
global spectrum
global spectrum_c

% Constants
c = 2.99792458e8;       % m/s speed of light
h = 6.62606957e-34; 	% Js Planck's constant

lambda_c=spectrum_c(round(length(spectrum_c)/2)); %puts the sample wavelength for n_cauchy in middle of the defined cauchy region

%expand k data for KK with gap estimate (Tauc equation)
alpha_crit=4*pi*k_data(1)/spectrum(1)/1e-7; %Use absorption coefficient at a point in the bandgap to determine A constant for full absorption coef. estimation alpha_ll
Eg=h*c/spectrum_c(1);
Eg2=h*c/(spectrum(1)-50);
A=alpha_crit/sqrt(h*c/spectrum(1)-Eg);

ls=length(spectrum);
spc=step:step:spectrum(ls);
lex=length(step:step:spectrum(1)-step);

k_extra=k_data(1)*ones(1,lex);
i=0;
for ll=step:step:spectrum(1)-step
    i=i+1;
    alpha_ll=(A*sqrt(h*c/ll-Eg))+real(A*sqrt(h*c/ll-Eg2)); % this assumes a simplified deeper valence band structure, with an extra band 50nm below lower bound of the spectral range                                                                                             
    k_extra(i)=ll*1e-7*alpha_ll/4/pi;   %expands k spectrally past the available Transmission data, using the absorption coefficient estimation
end
k_init=zeros(1,length(k_data)+length(k_extra)); %full spectral k up to 1nm wavelength
for i=1:length(k_extra)
    k_init(i)=k_extra(i);
end
for i=length(k_extra)+1:length(k_data)+length(k_extra)
    k_init(i)=k_data(i-length(k_extra));
end

spec_original=spectrum;  %expand spectrum
spectrum=spc;
ls=length(spectrum);

%Calculate off-set with the known cauchy value to determine integral
%constant of KK
j=0;
sum_c=0;
for lamb=spectrum(1):spacing_KK:lambda_c-spacing_KK
    j=j+1;
    sum_c=sum_c+4*c*(2*pi*c)/lamb^3*k_init(j)/((2*pi*c/lamb)^2-(2*pi*c/lambda_c)^2)*spacing_KK; %note lambda_c
end
j=j+1; %skips the lambda_c to prevent infinite value (special type of integral for KK)
for lamb=lambda_c+spacing_KK:spacing_KK:spectrum(ls)
    j=j+1;
    sum_c=sum_c+4*c*(2*pi*c)/lamb^3*k_init(j)/((2*pi*c/lamb)^2-(2*pi*c/lambda_c)^2)*spacing_KK;
end
cte=n_cauchy+B/(lambda_c^2)-1-sum_c;  %where B in the non-constant cauchy factor

%Calculate n from Kramer's Kronig
i=0;
for lamb_i=spectrum(1):spacing_KK:spectrum(ls)
    i=i+1;
    sum_i=0;
    j=0;
    for lamb=spectrum(1):spacing_KK:lamb_i-spacing_KK
        j=j+1;
        sum_i=sum_i+4*c*(2*pi*c)/lamb^3*k_init(j)/((2*pi*c/lamb)^2-(2*pi*c/lamb_i)^2)*spacing_KK;
    end
    j=j+1; %skips lamb_i
    for lamb=lamb_i+spacing_KK:spacing_KK:spectrum(ls)
        j=j+1;
        sum_i=sum_i+4*c*(2*pi*c)/lamb^3*k_init(j)/((2*pi*c/lamb)^2-(2*pi*c/lamb_i)^2)*spacing_KK;
    end
    n_ex(i)=1+sum_i+cte;
end

spectrum=spec_original; %go back to original spectrum
ls=length(spectrum);
lnex=length(n_ex);
n_out=n_ex(lex+1:lnex); %final rough n