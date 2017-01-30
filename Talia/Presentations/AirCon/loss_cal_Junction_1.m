clear;
clc;

%add some path and imput data
addpath(genpath('../'));

AirCon = importdata('spectrum_output.csv');

%get wavelength and convert energy from nm to eV
wavelength = AirCon.data(:,1);

Energy_total = 6.626e-34 * 3.0e8 / 1.6e-19 ./ AirCon.data(:,1) / 1e-6;

%InGaN band gap range
Energy_gap = 3.4:-0.01:0.7;

%set some variables
m = length(AirCon.textdata);

n = length(Energy_gap);

results = zeros(n,2); 

results_j = zeros(2,m); %space for results

for j = 2:m % j is the serial number for each air condition
    
    for i = 1:n
    
        %left part calculation in the equation
        x1 = Energy_gap(i);
    
        Energy_loss_left = Energy_total - x1;
    
        negative = Energy_loss_left < 0;
    
        Energy_loss_left(negative) = 0; %remove negative values
    
        Power_loss_left = AirCon.data(:,j);
    
        Power_loss_left(negative) = 0; %remove negative values
    
        temp = Power_loss_left .* Energy_loss_left ./ Energy_total;
    
        Integ_loss_left = trapz(wavelength, temp);
    
    
        %right part calculation in the equation
        Power_loss_right = AirCon.data(:,j);
    
        Power_loss_right(~negative) = 0;
    
        Integ_loss_right = trapz(wavelength, Power_loss_right);
    
        Integ_loss_sum = Integ_loss_left + Integ_loss_right;
    
        results(i,1) = x1;
        results(i,2) = Integ_loss_sum;
        
        [Integration_sum, Index_sum_1] = min( squeeze(results(:,2)));

        Band_gap_1 = results(Index_sum_1,1);
        
    end
    
    %save results
    results_j(:,j) = [Band_gap_1,Integration_sum];
    
    
end

save result_1.mat