clear;
clc;

%add some path 
addpath(genpath('../'));

AirCon = importdata('spectrum_output.csv');

wavelength = AirCon.data(:,1);

Energy_total = 6.626e-34 * 3.0e8 / 1.6e-19 ./ AirCon.data(:,1) / 1e-6;

plot(wavelength * 1000 , AirCon.data(:,2:13),':k',':k',':k',':k',':k',':k',':k',':k',':k',':k',':k',':k');
xlabel('Wavelength (nm)');
ylabel('Spectral irradiance (W/(m^2��nm)');
legend(AirCon.textdata(2:13));

Energy_gap = 0.31:0.01:3.4;

m = length(AirCon.textdata);

n = length(Energy_gap);

results = zeros(n,2); 

results_j = zeros(2,m); %space for results

for j = 2:m
    
    for i = 1:n
    
        x1 = Energy_gap(i);
    
        Energy_loss_left = Energy_total - x1;
    
        negative = Energy_loss_left < 0;
    
        Energy_loss_left(negative) = 0; %remove negative values
    
        Power_loss_left = AirCon.data(:,j);
    
        Power_loss_left(negative) = 0; %remove negative values
    
        temp = Power_loss_left .* Energy_loss_left ./ Energy_total;
    
        Integ_loss_left = trapz(wavelength, temp);
    
    
        Power_loss_right = AirCon.data(:,j);
    
        Power_loss_right(~negative) = 0;
    
        Integ_loss_right = trapz(wavelength, Power_loss_right);
    
        Integ_loss_sum = Integ_loss_left + Integ_loss_right;
    
        results(i,1) = x1;
        results(i,2) = Integ_loss_sum;
        
        [Integration_sum, Index_sum_1] = min( squeeze(results(:,2)));

        Band_gap_1 = results(Index_sum_1,1);
        
    end
    
    results_j(:,j) = [Band_gap_1,Integration_sum];
    
    
end

