clear;
clc;

%add some path 
addpath(genpath('../'));

AirCon = importdata('spectrum_output.csv');

wavelength = AirCon.data(:,1);

Energy_total = 6.626e-34 * 3.0e8 / 1.6e-19 ./ AirCon.data(:,1) / 1e-6;

Energy_gap = 0.31:0.01:3.4;

m = length(AirCon.textdata);

n = length(Energy_gap);

count_i = 0;

results = zeros((n*(n+1)/2),3); 

results_l = zeros(3,m); %space for results

for l = 2:m
    
    for i = 1:n
    
        x1 = Energy_gap(i);
        
        count_j = 1;
        
        for j = i:n
            
            % this line is only used for debugging
            fprintf('i = %03d -- j = %03d\n',i,j)
            x2 = Energy_gap(j);
            
            %left part
    
            Energy_loss_left = Energy_total - x1;
    
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
            
            results((count_i * n + count_j),:) = [x1,x2,Integ_loss_sum];
            
            count_j = count_j + 1;
        end
        
        count_i = count_i + 1;
        
    end
    
    [Integration_sum, Index_sum] = min( squeeze(results(:,3)));
    
    Band_gap_1 = results(Index_sum,1);
    Band_gap_2 = results(Index_sum,2);
    results_l(:,l) = [Band_gap_1,Band_gap_2,Integration_sum];

end


