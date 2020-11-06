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
%   This script implements a global thickness fitting in the Cauchy Regime,
%   using Transfer Matrix calculations.
%	The mathematical procedure for calculating the optical field via the transfer matrix method 
%	was taken from:
%   Harbecke B., "Coherent and incoherent reflection and transmission of multilayer structures." Applied Physics B 1986, 39:165-170.

%% new instance of GtFinder Script
function [t_sample1,t_sample2,nc_firstorder,B_f,offset] = GT_Finder()  

global step
global spectrum
global spectrum_c
global Rexp
global Rexp2
global R_angle
global n_c_range
global thickness_range
global thickness_range2
global layers
global thicknesses
global active_layer
global nk_filename
global B_cauchy_range

ls=length(spectrum);
lsc=length(spectrum_c);

[t,t2]=deal(thicknesses);

% Load in index of refraction for each material before active layer
ntotal = zeros(length(layers),ls);
for index = 1:active_layer
    ntotal(index,:) = LoadRefrIndex(layers{index},spectrum,nk_filename);
end

%remaining refractive indices after active layer
for index = active_layer+2:length(layers)
    ntotal(index,:) = LoadRefrIndex(layers{index},spectrum,nk_filename);
end

%% Compare the R,T data in the cauchy regime to varying Transfer Matrix t and n_cauchy, to derive the sample thickness and n_cauchy

Error=42; %meaningless high error value

for nc=n_c_range
   for B=1000000*B_cauchy_range
    [R1,R2]=deal(zeros(1,lsc));    %allocate Rs,Rp for speed
    for t_act=thickness_range 
    t(active_layer)=t_act;
    for t_act2=thickness_range2
        t2(active_layer)=t_act2;    
        i=0; 
        for lambda=spectrum_c
            ntotal(active_layer+1,:)=nc+B/(lambda^2);
            n = ntotal(:,lambda/step-spectrum(1)/step+1).';
            i=i+1;
            R1(i) = Tmat(lambda,n,t,R_angle(1));
            R2(i) = Tmat(lambda,n,t2,R_angle(2));
        end
        Dmean= mean(Rexp(spectrum_c/step-spectrum(1)/step+1))-mean(R1);
        Dmean2= mean(Rexp2(spectrum_c/step-spectrum(1)/step+1))-mean(R2);
          R1=R1+Dmean;
          R2=R2+Dmean2;
        err=norm((R1-Rexp(spectrum_c/step-spectrum(1)/step+1)))...
        +norm((R2-Rexp2(spectrum_c/step-spectrum(1)/step+1)));
        if(err<Error)
            Error=err;
            t_active=t_act;
            t_active2=t_act2;
            n_cauchy=nc;
            Dmean_f=Dmean;
            Dmean_f2=Dmean2;
            B_f=B;
        end
    end
    end
   end
end

%print out the findings
t_sample1=t_active %thickness of sample 1
t_sample2=t_active2 %thickness of sample 2
nc_firstorder=n_cauchy %First order Cauchy term

offset=Dmean_f %average deviation of T data
B_cauchy=B_f/1000000 %Second order Cauchy term


%Following can be used to see quality of fit
figure('name','Cauchy Regime Model to Data quality match')

t(active_layer)=t_active;
i=0;
    for lambda=spectrum_c
        i=i+1;
        ntotal(active_layer+1,:)=n_cauchy+B_f/(lambda^2);
        n = ntotal(:,lambda/step-spectrum(1)/step+1).';
        [~,T(i)] = Tmat(lambda,n,t,R_angle(1));
        [R(i),~] = Tmat(lambda,n,t,R_angle(2));
        T(i)=T(i)-Dmean_f;
        R(i)=R(i)+Dmean_f;
    end
t(active_layer)=t_active2;
i=0;
    for lambda=spectrum_c
        i=i+1;
        ntotal(active_layer+1,:)=n_cauchy+B_f/(lambda^2);
        n = ntotal(:,lambda/step-spectrum(1)/step+1).';
        [~,T2(i)] = Tmat(lambda,n,t,R_angle(1));
        [R2(i),~] = Tmat(lambda,n,t,R_angle(2));
        T2(i)=T2(i)-Dmean_f2;
        R2(i)=R2(i)+Dmean_f2;
    end

%Quality of Fit plot
orange=1/255*[255, 165, 0];
plot(spectrum_c, T,'b--','LineWidth',1)
hold on
plot(spectrum_c, 1-Rexp(spectrum_c/step-spectrum(1)/step+1),'b','LineWidth',2)
plot(spectrum_c,T2,'--','Color',orange,'LineWidth',1)
plot(spectrum_c,1-Rexp2(spectrum_c/step-spectrum(1)/step+1),'Color',orange,'LineWidth',2)
ylabel('Transmission')
xlabel('Wavelength [nm]')
legend('Model_T S1', 'Data_T S1', 'Model_T S2', 'Data_T S2')
hold off
drawnow

end

%% Function LoadRefrIndex 
% The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function ntotal = LoadRefrIndex(name,wavelengths,file)

%Data in IndRefr, Column names in IndRefr_names
[IndRefr,IndRefr_names]=xlsread(file);

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr(:,strmatch(strcat(name,'_lambda'),IndRefr_names));
n=IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
k=IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));   

Nan=find(isnan(file_wavelengths));
if(~isempty(Nan))
    file_wavelengths=file_wavelengths(1:Nan(1)-1);
    n=n(1:Nan(1)-1);
    k=k(1:Nan(1)-1);
end


% Interpolate/Extrapolate data linearly to desired wavelengths
n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
ntotal = n_interp+1i*k_interp; 
end