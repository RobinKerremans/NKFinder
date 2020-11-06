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
%   This script aims to approximate a material's optical constants, 
%   derived from transmission data of the sample on a transparent substrate. 
%   (2 samples of different thicknesses are required).
%   To do so it uses automated global fitting in the Cauchy Regime to determine
%   the layer thicknesses, as well as transfer matrix and Kramers Kronig
%   calculations to find n an k.
%
%   Inputs for the spectral range, data files and other variables are
%   entered manually at the top (see commentary)
%   
%   For a more detailed description and for citation please refer to: 
%   "The Optical Constants of Solution Processed Semiconductors – New Challenges with Perovskites and Non-Fullerene Acceptors"
%   Robin Kerremans, Christina Kaiser, Wei Li, Nasim Zarrabi, Paul Meredith* and Ardalan Armin, Advanced Optical Materials, DOI:10.1002/adom.202000319, 2020. 
%_________________________________________________________________________
%	The mathematical procedure for calculating the optical field via the transfer matrix method 
%	was taken from:
%   Harbecke B., "Coherent and incoherent reflection and transmission of multilayer structures." Applied Physics B 1986, 39:165-170.

function NKfinder()  
%% new instance of NKfinder Script
close all
clc
clear

global step
global spectrum
global spectrum_c
global Texp
global Texp2
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

%%====================================================================================================================
%% Manual Input Section (Initiating variables)
%%====================================================================================================================

step = 5;                            % Wavelength step-size of the R,T measurement
spectrum=370:step:1600;            % Examined spectrum, in nm. Recommended to start above glass absorption edge, 350nm
spectrum_c=800:step:1600;          % Cauchy regime spectrum, in nm. Any run will display a T data plot so cauchy regime can be estimated.

ns=0.05;                         % Stepsize when looking for n_cauchy. Recommended at 0.1-0.05.
n_c_range= 1.5:ns:2.5;          % Range and sensitivity when looking for n_cauchy.
k_range= 0:0.001:2.5;             % Range and sensitivity when looking for k in general

ts=5;                            % Stepsize when looking for sample thickness. Recommended at 5-10 for final run, but can be larger for estimates of the range.

%A data fit will be displayed after every completion of Global Thickness finder,
%and can be used to adjust range and spectrum_c accordingly

thickness_range =145:ts:165;     % Range when looking for thickness of sample 1. Recommended sample thickness = 120nm-400nm
thickness_range2 = 105:ts:135;    % Range when looking for thickness of sample 2. Can be set 0:ts:0 if only 1 sample.
B_cauchy_range = 0:0.1:0.3;       % This parameter allows a non-constant Cauchy Regime approximation (n_cauchy = nc + B/(lambda^2) in units of cm2

Texp=LoadTransmission(spectrum,'PCE12 IT4F 2.csv');         % Filenames of the T csv data files IN PERCENTAGE (in same map as script)
Texp2=LoadTransmission(spectrum,'PCE12 IT4F 1.csv');
   
Output_file='PCE12_IT4F.csv';                          % Name of the output csv file (will also show n_cauchy and thickness in this file)
loopnum=1;                                            % Number of iterations. 1-2 is usually enough.

%%====================================================================================================================
%% Only ADVANCED input below here (Using Reflection data or inputting more layers)
%%====================================================================================================================

layers = {'Air' 'x' 'Glass_AA' 'Air'};         % Names of layers of materials starting from side light is incident from during measurement
                                            % Names of layers must match the excel n,k source file
                                            % 'x' marks the position of the unknown NK material.
                                            
thicknesses = [0 7000000];        % Thicknesses of layers on substrate in nm. Substrate is treated as incoherent if >100,000nm
                                  % Note  this leaves out the medium (Air)
                                  % Unknown material thickness is initiated as '0' by default
                                                   
active_layer=1;                   % in reference to 'thicknesses' vector, the index of the layer to be examined ('0')                                                    

nk_filename='Optical_Constants.xls';        % filename of the n,k excel data for the known layers (usually only substrate and medium)

Rexp=ones(1,length(Texp))-Texp;                        % R can be derived from T in cauchy regime
Rexp2=ones(1,length(Texp))-Texp2; 
%Rexp2=LoadReflection(spectrum,'PCP 800 R.csv');        % R data file can also be referenced IN ADDITION to T data.  
                                                              % This enhances GT_Finder only (poor sensitivity to k)
                                                              
%disp('R DATA INPUT IS ON')                                    % Don't forget to comment out again, or to restore the angle values! 

R_angle=[0,0].*pi/180;                      % Incident Angle of the measurement, default 0 degrees for Transmission (R_angle(1)) and 0 degrees for Reflectance (R_angle(2))

%% Manual Input ENDS here =================================================================================================================================================

figure('name','Transmission Data (Cauchy Regime identification)')
plot(spectrum,Texp,spectrum,Texp2,'LineWidth',2)              %Used to display R,T data to pinpoint Cauchy region.
legend('T data 1','T data 2')
ylabel('Transmission')
xlabel('Wavelength [nm]')
drawnow

%open Cauchy Dialog box
answer = questdlg('Continue?','Cauchy regime pause','Yes','No','Yes');
if(strmatch(answer,'No'))
    return
end

ls=length(spectrum);
lsc=length(spectrum_c);

%make sure two spectral ranges match
if(spectrum(ls)~=spectrum_c(lsc))
    disp('Warning: Upper bound of cauchy regime range does not match upper bound of spectrum.')
end
t=thicknesses;

lambda_c=spectrum_c(round(length(spectrum_c)/2)); %puts the sample wavelength for n_cauchy in middle of the defined cauchy region

% Load in index of refraction for each material before active layer
ntotal = zeros(length(layers),ls);
for index = 1:active_layer
    ntotal(index,:) = LoadRefrIndex(layers{index},spectrum,nk_filename);
end

%remaining refractive indices after active layer
for index = active_layer+2:length(layers)
    ntotal(index,:) = LoadRefrIndex(layers{index},spectrum,nk_filename);
end


%% Global thickness finder. Finds n_cauchy and sample 1&2 thicknesses, with a spectral data fit in the Cauchy Regime.

[t_init,t2,nc,B,offset]= GT_Finder();  %Ouputs thickness, n cauchy and 'Dmean_f', which is the average offset due to scattering, baseline, ...
nstring=num2str(nc+B/lambda_c^2);
tstring=num2str(t_init);
tstring2=num2str(t2);
Dialog=strcat('ncauchy=',nstring,', t1=',tstring,',nm t2=',tstring2,'nm. Continue?');

%open Global Thickness Dialog box
answer = questdlg(Dialog,'Thickness Finder pause','Yes','No','Yes');
if(strmatch(answer,'No'))
    return
end

%% Refine the thickness for 1 sample only (arbitrarily the first). 
%Final thickness and n_cauchy initialization. 

Error=10; %arbitrary high error value
for ncauchy=nc-ns:0.01:nc+ns
for t_act=t_init-round(ts/2):1:t_init+round(ts/2) %nc + thickness finder
    t(active_layer)=t_act;
    R=zeros(1,lsc);    %allocate Rs,Rp for speed
        i=0;
        for lambda=spectrum_c
            ntotal(active_layer+1,:)=ncauchy + B/(lambda^2); %note B is in micrometer^2
            n = ntotal(:,lambda/step-spectrum(1)/step+1).';
            i=i+1;
            R(i) = Tmat(lambda,n,t,R_angle(2)); 
        end
        R=R+offset;
        err=norm((R-Rexp(spectrum_c/step-spectrum(1)/step+1)).^5);
        if(err<Error)
            Error=err;
            t_active_f=t_act;
            n_c_f=ncauchy;
        end   
end
end

t(active_layer)=t_active_f;
t_refined=t_active_f
n_cauchy=n_c_f+B/(lambda_c^2)  %this is just to signal the actual n_cauchy, including B
n_cauchy=n_c_f;

Error=ones(1,ls); %return Error to an arbitrary high value
k_active=zeros(1,ls);
% n_test=zeros(1,ls);

%% Rough n,k constant determination, using n_cauchy and a free k variable. Comparison through Transfer Matrix Transmission simulation.

%k iteration
n_act=n_cauchy;
for lambda=spectrum(1):step:spectrum_c(1)
    for k_act=k_range
        ntotal(active_layer+1,lambda/step-spectrum(1)/step+1)=n_act+B/(spectrum_c(1)^2)+1i*k_act;
        n = ntotal(:,lambda/step-spectrum(1)/step+1).';
        [R,T] = Tmat(lambda,n,t,R_angle(1));
        %T=T-offset; %baseline/scattering factor
        er=abs(T-Texp(lambda/step-spectrum(1)/step+1));
        if(er<Error(lambda/step-spectrum(1)/step+1))
            Error(lambda/step-spectrum(1)/step+1)=er;
            k_active(lambda/step-spectrum(1)/step+1)=k_act;
%             n_test(lambda/step-spectrum(1)/step+1)=n(2);
        end
    end
end

%% Kramer's Kronig for accurate n determination from k

k_data=k_active;
spacing_KK=step;

n_1=KK(n_cauchy,k_data,spacing_KK,B); %Kramers Kronig function, gives n from k

%plot n and k and output the results to a csv file
figure('name',['Optical Constants of ',Output_file,' - rough'])
plot(spectrum,k_active,spectrum,n_1)
legend('k Rough','n Rough')
ylabel('Optical Constants')
xlabel('Wavelength [nm]')
drawnow

%% Fine-tunes the found constants by repeating the previous steps, but now with the more accurate n
for loop=1:loopnum
%complete scattering factor with knowledge of k outside cauchy regime
j=0;
n_act=n_1;
Error=ones(1,ls);
for lambda=spectrum(1):step:spectrum_c(1)
    j=j+1;
    for k_act=k_range
        ntotal(active_layer+1,lambda/step-spectrum(1)/step+1)=n_act(j)+1i*k_act;
        n = ntotal(:,lambda/step-spectrum(1)/step+1).';
        [R,T] = Tmat(lambda,n,t,R_angle(1));
        %T=T-offset;
        err=abs(T-Texp(lambda/step-spectrum(1)/step+1));
        if(err<Error(lambda/step-spectrum(1)/step+1))
            Error(lambda/step-spectrum(1)/step+1)=err;
            k_active(lambda/step-spectrum(1)/step+1)=k_act;
        end
    end
end
%% Kramer's Kronig for accurate n determination from k

k_data=k_active;
spacing_KK=step;

n_1=KK(n_cauchy,k_data,spacing_KK,B); %Kramers Kronig function, gives n from k

%plot n and k and output the results to a csv file
str=num2str(loop);
figure('name',['Optical Constants of ',Output_file,' - Loop=',str])
plot(spectrum,k_active,spectrum,n_1,'LineWidth',2)
legend('k Final','n Final')
ylabel('Optical Constants')
xlabel('Wavelength [nm]')
drawnow
tncauchy=zeros(ls,1);
tncauchy(1)=t(active_layer);
tncauchy(2)=n_cauchy+B/(lambda_c^2);
String_fine=['Optical Constants of ',Output_file];
M=[spectrum',n_1',k_active',tncauchy];
csvwrite(String_fine,M)  %make sure old csv file of same name is closed
end

end


%% Function LoadRefrIndex 
% The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function ntotal = LoadRefrIndex(name,wavelengths,file)

%Data in IndRefr, Column names in IndRefr_names
[IndRefr,IndRefr_names]=xlsread(file);

if(~strmatch(name,IndRefr_names))
disp('Error occured: Material name does not match the name in Optical Constants file.')
end

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

function Refl = LoadReflection(wavelengths,file)

M=csvread(file,1,0);
file_wavelengths=M(:,1);
file_R=1/100*M(:,2);

mw=max(wavelengths);
mfw=max(file_wavelengths);
if(mw>mfw)
disp('Error occured: Specified Spectrum range exceeds Spectral range in Reflection data file.')
return
end

% Interpolate/Extrapolate data linearly to desired wavelengths
Refl=interp1(file_wavelengths, file_R, wavelengths, 'linear', 'extrap');
end

function Trans = LoadTransmission(wavelengths,file)

M=csvread(file,1,0);
file_wavelengths=M(:,1);
file_T=1/100*M(:,2);

mw=max(wavelengths);
mfw=max(file_wavelengths);
if(mw>mfw)
disp('Error occured: Specified Spectrum range exceeds Spectral range in Transmission data file.')
return
end

% Interpolate/Extrapolate data linearly to desired wavelengths
Trans=interp1(file_wavelengths, file_T, wavelengths, 'linear', 'extrap');
end