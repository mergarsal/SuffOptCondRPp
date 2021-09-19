function Ps = createScaleMatrices()

% create Ps 
Ps = zeros(6, 5, 5); 

% general matrix
P = zeros(6); 
% P(1,1) = 2; 
% P(3,3) = 2; 
% P(5,5) = 2; 
% P(2,2) = 0.5; 
% P(4,4) = 0.5; 
% P(6,6) = 0.5; 
% P(1, 3) = 2; 
% P(1, 5) = 2; 
% P(3, 5) = 2; 

P(1,1) = 2; % lambda 2
P(4,4) = 2; % lambda 5
P(6,6) = 2; % lambda 7
P(2,2) = 0.5;  % lambda 3
P(3,3) = 0.5;  % lambda 4
P(5,5) = 0.5;  % lambda 6
P(1, 4) = 2; 
P(1, 6) = 2; 
P(4, 6) = 2; 

P = 0.5 * (P + P'); 

% we discard lambda 2
idx = [2:6]; 
Pi = P(idx, idx); 
[up, dp]= eig(Pi); 
PPi = sqrt(dp) * up'; 
Ps(1, :, :) = pinv(PPi); 


% we discard lambda 3
clear Pi idx up dp PPi; 
idx = [1, 3:6]; 
Pi = P(idx, idx); 
[up, dp]= eig(Pi); 
PPi = sqrt(dp) * up';
Ps(2, :, :) = pinv(PPi); 


% we discard lambda 4
clear Pi idx up dp PPi; 
idx = [1, 2, 4:6]; 
Pi = P(idx, idx); 
[up, dp]= eig(Pi); 
PPi = sqrt(dp) * up';
Ps(3, :, :) = pinv(PPi);


% we discard lambda 5
clear Pi idx up dp PPi; 
idx = [1:3, 5:6]; 
Pi = P(idx, idx); 
[up, dp]= eig(Pi); 
PPi = sqrt(dp) * up';
Ps(4, :, :) = pinv(PPi);


% we discard lambda 6
clear Pi idx up dp PPi; 
idx = [1:4, 6]; 
Pi = P(idx, idx); 
[up, dp]= eig(Pi); 
PPi = sqrt(dp) * up';
Ps(5, :, :) = pinv(PPi);

% we discard lambda 7
clear Pi idx up dp PPi; 
idx = [1:5]; 
Pi = P(idx, idx); 
[up, dp]= eig(Pi); 
PPi = sqrt(dp) * up';
Ps(6, :, :) = pinv(PPi);

end
