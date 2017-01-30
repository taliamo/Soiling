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

results_l = zeros(6,m); %space for results

for l = 2:m
    
min_lost = inf;
engy1 = 0;
engy2 = 1;
engy3 = 2;
engy4 = 3;
engy5 = 4;

    for i = 1:n
    
        x1 = Energy_gap(i);
        
        for j = i+1:n
            
            x2 = Energy_gap(j);
            
            for k = j+1:n
                
                x3 = Energy_gap(k);
                
                for h = k+1:n
                    
                    x4 = Energy_gap(h);
                    
                    for g = h+1:n
                        
                        x5 = Energy_gap(g);
            
                %left part
    
                Energy_loss_left = Energy_total(:,1) - x1;
    
                negative1 = Energy_loss_left < 0;
    
                Energy_loss_left(negative1) = 0; %remove negative values
    
                Power_loss_left = AirCon.data(:,l);
    
                Power_loss_left(negative1) = 0; %remove negative values
    
                temp = Power_loss_left .* Energy_loss_left ./ Energy_total;
    
                Integ_loss_left = trapz(wavelength, temp);
    
    
                % middle part 1
            
                Energy_loss_mid = Energy_total(:,1) - x2;
            
                negative2 = Energy_loss_mid < 0;
            
                Energy_loss_mid(~negative1) = 0;
                Energy_loss_mid(negative2) = 0; % remove values
            
                Power_loss_mid = AirCon.data(:,l);
                Power_loss_mid(~negative1) = 0;
                Power_loss_mid(negative2) = 0; %remove values
            
                temp = Power_loss_mid .* Energy_loss_mid ./ Energy_total;
            
                Integ_loss_mid = trapz(wavelength, temp);
                
                % middle part 2
                
                Energy_loss_mid_2 = Energy_total(:,1) - x3;
                
                negative3 = Energy_loss_mid_2 < 0;
                
                Energy_loss_mid_2(~negative1) = 0;
                Energy_loss_mid_2(~negative2) = 0;
                Energy_loss_mid_2(negative3) = 0; % remove negative values
                
                Power_loss_mid_2 = AirCon.data(:,l);
                Power_loss_mid_2(~negative1) = 0;
                Power_loss_mid_2(~negative2) = 0;
                Power_loss_mid_2(negative3) = 0;
                
                temp = Power_loss_mid_2 .* Energy_loss_mid_2 ./ Energy_total;
                
                Integ_loss_mid_2 = trapz(wavelength, temp);
                
                % middle part 3
                
                Energy_loss_mid_3 = Energy_total(:,1) - x4;
                
                negative4 = Energy_loss_mid_3 < 0;
                
                Energy_loss_mid_3(~negative1) = 0;
                Energy_loss_mid_3(~negative2) = 0;
                Energy_loss_mid_3(~negative3) = 0;
                Energy_loss_mid_4(negative4) = 0;
                
                Power_loss_mid_3 = AirCon.data(:,l);
                Power_loss_mid_3(~negative1) = 0;
                Power_loss_mid_3(~negative2) = 0;
                Power_loss_mid_3(~negative3) = 0;
                Power_loss_mid_3(negative4) = 0;
                
                temp = Power_loss_mid_3 .* Energy_loss_mid_3 ./ Energy_total;
                
                Integ_loss_mid_3 = trapz(wavelength, temp);
                
                % middle part 4
                
                Energy_loss_mid_4 = Energy_total(:,1) - x5;
                
                negative5 = Energy_loss_mid_4 < 0;
                
                Energy_loss_mid_4(~negative1) = 0;
                Energy_loss_mid_4(~negative2) = 0;
                Energy_loss_mid_4(~negative3) = 0;
                Energy_loss_mid_4(~negative4) = 0;
                Energy_loss_mid_5(negative5) = 0;
                
                Power_loss_mid_4 = AirCon.data(:,l);
                Power_loss_mid_4(~negative1) = 0;
                Power_loss_mid_4(~negative2) = 0;
                Power_loss_mid_4(~negative3) = 0;
                Power_loss_mid_4(~negative4) = 0;
                Power_loss_mid_4(negative5) = 0;
                
                temp = Power_loss_mid_4 .* Energy_loss_mid_4 ./ Energy_total;
                
                Integ_loss_mid_4 = trapz(wavelength, temp);
                
                
                 % right part
            
                Power_loss_right = AirCon.data(:,l);
               
                Power_loss_right(~negative5) = 0;
    
                Integ_loss_right = trapz(wavelength, Power_loss_right);
    
                Integ_loss_sum = Integ_loss_left + Integ_loss_mid + Integ_loss_mid_2 + Integ_loss_mid_3 + Integ_loss_mid_4 + Integ_loss_right;
            
                if min_lost > Integ_loss_sum
                    
                    min_lost = Integ_loss_sum;
                    
                    engy1 = x1;
                    engy2 = x2;
                    engy3 = x3;
                    engy4 = x4;
                    engy5 = x5;
                    
                end
                
                    end
                    
                end
                
            end
            
        end
        
    end
    
    results_l(:,l) = [engy1,engy2,engy3,engy4,engy5,min_lost];
end
        