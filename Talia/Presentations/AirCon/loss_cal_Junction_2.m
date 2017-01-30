clear;
clc;

%add some path 
addpath(genpath('../'));

AirCon = importdata('spectrum_output.csv');

wavelength = AirCon.data(:,1);

Energy_total = 6.626e-34 * 3.0e8 / 1.6e-19 ./ AirCon.data(:,1) / 1e-6;

Energy_gap = 3.4:-0.01:0.7;

m = length(AirCon.textdata);

n = length(Energy_gap);

results_l = zeros(3,m); %space for results

for l = 2:m
    
    min_lost = inf;
    engy1 = 0;
    engy2 = 1;

    for i = 1:n
    
        x1 = Energy_gap(i);
        
        for j = (i+1):n
            
            x2 = Energy_gap(j);
            
            %left part
    
            Energy_loss_left = Energy_total(:,1) - x1;
    
            negative1 = Energy_loss_left < 0;
    
            Energy_loss_left(negative1) = 0; %remove negative values
    
            Power_loss_left = AirCon.data(:,l);
    
            Power_loss_left(negative1) = 0; %remove negative values
    
            temp = Power_loss_left .* Energy_loss_left ./ Energy_total;
    
            Integ_loss_left = trapz(wavelength, temp);
    
    
            % middle part
            
            Energy_loss_mid = Energy_total(:,1) - x2;
            
            negative2 = Energy_loss_mid < 0;
            
            Energy_loss_mid(~negative1) = 0;
            Energy_loss_mid(negative2) = 0; % remove values
            
            Power_loss_mid = AirCon.data(:,l);
            Power_loss_mid(~negative1) = 0;
            Power_loss_mid(negative2) = 0; %remove values
            
            temp = Power_loss_mid .* Energy_loss_mid ./ Energy_total;
            
            Integ_loss_mid = trapz(wavelength, temp);
            
            % right part
            
            Power_loss_right = AirCon.data(:,l);
               
            Power_loss_right(~negative2) = 0;
    
            Integ_loss_right = trapz(wavelength, Power_loss_right);
    
            Integ_loss_sum = Integ_loss_left + Integ_loss_mid + Integ_loss_right;
            
            if min_lost > Integ_loss_sum
                min_lost = Integ_loss_sum;
                engy1 = x1;
                engy2 = x2;
            end
        end
    end
            
    results_l(:,l) = [engy1,engy2,min_lost];

end


