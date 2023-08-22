function result = helper_convert_trilayer_diagonal(type, params)
    if type == "delta2onsite"
        Delta1 = params(1);
        Delta2 = params(2);
        delta = params(3);
        
        result = helper_convert_delta2onsite(Delta1, Delta2, delta);
    elseif type == "onsite2delta"
        onsite_A1 = params(1);
        onsite_A2 = params(2);
        onsite_A3 = params(3);
        onsite_B1 = params(4);
        onsite_B2 = params(5);
        onsite_B3 = params(6);
        result = helper_convert_onsite2delta(onsite_A1, onsite_A2, onsite_A3, onsite_B1, onsite_B2, onsite_B3);
    else
        result = params;
    end
end


function result = helper_convert_delta2onsite(Delta1, Delta2, delta)
    onsite_A1 = Delta1 + Delta2;
    onsite_A3 = - Delta1 + Delta2;
    onsite_A2 = - 2 * Delta2 + delta;

    onsite_B1 = Delta1 + Delta2 + delta;
    onsite_B3 = - Delta1 + Delta2 + delta;
    onsite_B2 = - 2 * Delta2;
    
    result = zeros(6,1);
    result = [onsite_A1, onsite_A2, onsite_A3, onsite_B1, onsite_B2, onsite_B3];
end


function result = helper_convert_onsite2delta(onsite_A1, onsite_A2, onsite_A3, onsite_B1, onsite_B2, onsite_B3)
    Delta2 = (onsite_A1 + onsite_A3) / 2;
    Delta1 = (onsite_A1 - onsite_A3) / 2;
    delta = onsite_A1 + onsite_A2 + onsite_A3;

    Delta2_ = -onsite_B2 / 2;
    delta_ = (onsite_B1 + onsite_B2 + onsite_B3) / 2;
    Delta1_ = (onsite_B1 - onsite_B3) / 2; 
    
    result = zeros(6,1);
    result = [Delta1, Delta2, delta, Delta1_, Delta2_, delta_];
end


% % example
% delta = 35.8;
% Delta1 = 10;
% Delta2 = 18;
% 
% params = [Delta1, Delta2, delta];
% type = "delta2onsite";
% result = helper_convert_trilayer_diagonal(type, params);
% 
% params = result;
% type = "onsite2delta";
% result2 = helper_convert_trilayer_diagonal(type, params);