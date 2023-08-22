delta = 35.8;
Delta1 = 10;
Delta2 = 18;

% onsite_A1 = Delta1 + Delta2;
% onsite_A3 = - Delta1 + Delta2;
% onsite_A2 = - 2 * Delta2 + delta;
% 
% onsite_B1 = Delta1 + Delta2 + delta;
% onsite_B3 = - Delta1 + Delta2 + delta;
% onsite_B2 = - 2 * Delta2;
% 
% 
% Delta2 = (onsite_A1 + onsite_A3) / 2;
% Delta1 = (onsite_A1 - onsite_A3) / 2;
% delta = onsite_A1 + onsite_A2 + onsite_A3;
% 
% Delta2_ = -onsite_B2 / 2;
% delta_ = (onsite_B1 + onsite_B2 + onsite_B3) / 2;
% Delta1_ = (onsite_B1 - onsite_B3) / 2;

params = [Delta1, Delta2, delta];
type = "delta2onsite";
result = helper_convert_trilayer_diagonal(type, params);

params = result;
type = "onsite2delta";
result2 = helper_convert_trilayer_diagonal(type, params);