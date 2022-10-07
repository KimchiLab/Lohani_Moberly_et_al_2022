function[S1,S2]=computeBeerLamberCoeffs(lambda_EX,lambda_EM,lambda_R1,lambda_R2)
%% compute cofficients for hemodynamics subtraction based on 
%% simplified beer lambert model (Allen...Waters 2020, similar to Ma...Hillman 2016)

% load hillman pathlengths and hemoglobin abosrption coefficients
load('PrahlExtinctionCoeffs.mat');
load('HillmanPathlengths.mat');

%extinction coefficients for oxygenated and deoxygenated haemoglobin at specific wavelengths
EHbO_EM=.1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,2),lambda_EM,'pchip','extrap');
EHbO_EX=.1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,2),lambda_EX,'pchip','extrap');
EHbO_R1=.1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,2),lambda_R1,'pchip','extrap');
EHbO_R2=.1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,2),lambda_R2,'pchip','extrap');
EHbR_EM=.1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,3),lambda_EM,'pchip','extrap');
EHbR_EX=.1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,3),lambda_EX,'pchip','extrap');
EHbR_R1=.1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,3),lambda_R1,'pchip','extrap');
EHbR_R2=.1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,3),lambda_R2,'pchip','extrap');

%path lengths in um
x_EM= interp1(hillman_pathlengths.data(:,1),hillman_pathlengths.data(:,2),lambda_EM,'pchip','extrap')./2;%need to divide the excitation wavelenght path length by half according to Waters (2020)
x_EX= interp1(hillman_pathlengths.data(:,1),hillman_pathlengths.data(:,2),lambda_EX,'pchip','extrap')./2;%need to divide the emission wavelength path length by half according to Waters (2020) 
x_R1= interp1(hillman_pathlengths.data(:,1),hillman_pathlengths.data(:,2),lambda_R1,'pchip','extrap');
x_R2= interp1(hillman_pathlengths.data(:,1),hillman_pathlengths.data(:,2),lambda_R2,'pchip','extrap');

%solve S1 coefficient
Num1=(EHbO_EM*EHbR_R2)-(EHbR_EM*EHbO_R2);
Den=(EHbO_R1*EHbR_R2)-(EHbR_R1*EHbO_R2);
Num2=(EHbO_EX*EHbR_R2)-(EHbR_EX*EHbO_R2);
S1=((Num1/Den)*(x_EM/x_R1))+((Num2/Den)*(x_EX/x_R1));

%solve S2 coefficient
Num1=(EHbR_EM*EHbO_R1)-(EHbO_EM*EHbR_R1);
Den=(EHbO_R1*EHbR_R2)-(EHbR_R1*EHbO_R2);
Num2=(EHbR_EX*EHbO_R1)-(EHbO_EX*EHbR_R1);
S2=((Num1/Den)*(x_EM/x_R2))+((Num2/Den)*(x_EX/x_R2));


%alternative way to solve S1 and S2 coefficents based on Waters code (you get identical results to the above method but just keep here for sanity check)
%spectra.lambda=[lambda_EX,lambda_EM,lambda_R1,lambda_R2];
% % Interpolate extinction coeffs and Hillman pathlengths onto spectra wavelengths
% %   multiply by .1 to convert from /cm to /mm
% E(:,1) = .1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,2),spectra.lambda,'pchip','extrap');
% E(:,2) = .1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,3),spectra.lambda,'pchip','extrap');
% X = interp1(hillman_pathlengths.data(:,1),hillman_pathlengths.data(:,2),spectra.lambda,'pchip','extrap');
% % compute the "characteristic wavelength approximation"
% M=(-X.*E')';
% % we need to compress the excitation and emission (add .5*M_Fex + .5*M_Fem = M_F)
% numF=1; numR=2;
% AF = repelem(.5*eye(numF),1,2);
% AR = eye(numR);
% M = blkdiag(AF,AR)*M; % dimensions of M are now (S.numF+S.numR) x 2
% % separate fluor and refl blocks
% MF = M(1:numF,:);
% MR = M((numF+1):end,:);
% C = MF/MR; % [numF x numR]
% C = C'; % transpose to get [numR x numF]
%S1=C(1); S2=C(2); 
