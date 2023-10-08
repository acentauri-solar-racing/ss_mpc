    t = OptResNLP.xx1(:,3)/60/15;
    t_real = 0:1:length(weather.temp_data)-1;
    temp = par.temp_1*(t).^2 + par.temp_2*(t) + par.temp_3;
    plot(t,temp), hold on,
    plot(t_real, weather.temp_data(:,1))