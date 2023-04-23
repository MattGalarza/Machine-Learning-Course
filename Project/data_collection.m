% Number of samples
n = 250;

% Create the sample array
A = lhsdesign(n,8,'smooth','off'); 

% Upper and lower bounds
lbwt = 9e-6; % lower bound top electrode width
ubwt = 14e-6; % upper bound top electrode width
lbwb = 15e-6; % lower bound bottom electrode width
ubwb = 30e-6; % upper bound bottom electrode width
lblf = 300e-6; % lower bound electrode overlap length
ublf = 600e-6; % upper bound electrode overlap length
lbm = 1.25e-6; % lower bound shuttle mass
ubm = 2.25e-6; % upper bound shuttle mass
lbtfd = 25e-6; % lower bound device thickness
ubtfd = 50e-6; % upper bound device thickness
lbg0 = 10e-6; % lower bound initial gap
ubg0 = 20e-6; % upper bound initial gap
lbtp = 100e-9; % lower bound parylene layer thickness
ubtp = 200e-9; % upper bound parylene layer thickness
lbg = 0.1; % lower bound g-force
ubg = 1.5; % upper bound g-force

% Scale the sample array
A(:,1) = A(:,1).*(ubwt-lbwt)+lbwt;
A(:,2) = A(:,2).*(ubwb-lbwb)+lbwb;
A(:,3) = A(:,3).*(ublf-lblf)+lblf;
A(:,4) = A(:,4).*(ubm-lbm)+lbm;
A(:,5) = A(:,5).*(ubtfd-lbtfd)+lbtfd;
A(:,6) = A(:,6).*(ubg0-lbg0)+lbg0;
A(:,7) = A(:,7).*(ubtp-lbtp)+lbtp;
A(:,8) = A(:,8).*(ubg-lbg)+lbg;

% Run the simulation for the samples
iter = 0;
variables = 8;
var = 2*variables;
k = size(A,1);
B = zeros(k,1);
C = zeros(k,1);
D = zeros(k,var);
for i = 1:variables
    if i == 1
        disp('First set')
        for q = 1:k
            wt = A(q,1);
            wb = 30e-6;
            lf = 400e-6;
            m = 2.0933e-6;
            tfd = 25e-6;
            g0 = 14e-6;
            tp = 120e-9;
            g = 1;
            f = 20;
            
            runtime = 0.3; % simulation time
            sim('model',runtime)
            % Determines the max Vout and Vrms of the simulation
            tspan = 1/(2*f); % time for half cycle
            N = length(Vout);
            n = round(tspan/runtime*N); % number of datapoints in each cycle
            cycle = Vout((N-n):end); % last output cycle during simulation
            
            % Determine peaks from last output cycle
            peaks = zeros(10,1);
            tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
            for j = 1:size(tf,1)
                if tf(j) > 0
                    peaks(j) = tf(j);
                end
            end
            
            vrms = 0;
            numpeaks = 0;
            for t = 1:size(peaks,1)
                if peaks(t) > 0
                    vrms = vrms + (peaks(t))^2;
                    numpeaks = numpeaks + 1;
                end
            end
            X = sqrt(vrms/numpeaks);
            
            B(q) = X; % Vrms for given parameters
            C(q) = max(peaks); % Max Vout for given parameters
            D(2*q-q,i) = A(q,1);
            D(2*q-q,i+1) = X;
        end
    elseif i == 2
        disp('Second set')
        for q = 1:k
            wt = 9e-6;
            wb = A(q,2);
            lf = 400e-6;
            m = 2.0933e-6;
            tfd = 25e-6;
            g0 = 14e-6;
            tp = 120e-9;
            g = 1;
            f = 20;
            
            runtime = 0.3; % simulation time
            sim('model',runtime)
            % Determines the max Vout and Vrms of the simulation
            tspan = 1/(2*f); % time for half cycle
            N = length(Vout);
            n = round(tspan/runtime*N); % number of datapoints in each cycle
            cycle = Vout((N-n):end); % last output cycle during simulation
            
            % Determine peaks from last output cycle
            peaks = zeros(10,1);
            tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
            for j = 1:size(tf,1)
                if tf(j) > 0
                    peaks(j) = tf(j);
                end
            end
            
            vrms = 0;
            numpeaks = 0;
            for t = 1:size(peaks,1)
                if peaks(t) > 0
                    vrms = vrms + (peaks(t))^2;
                    numpeaks = numpeaks + 1;
                end
            end
            X = sqrt(vrms/numpeaks);
            
            B(q) = -X; % Vrms for given parameters
            C(q) = -max(peaks); % Max Vout for given parameters
            D(2*q-q,i+1) = A(q,2);
            D(2*q-q,i+2) = X;
        end
    elseif i == 3
        disp('Third set')
        for q = 1:k
            wt = 9e-6;
            wb = 30e-6;
            lf = A(q,3);
            m = 2.0933e-6;
            tfd = 25e-6;
            g0 = 14e-6;
            tp = 120e-9;
            g = 1;
            f = 20;
            
            runtime = 0.3; % simulation time
            sim('model',runtime)
            % Determines the max Vout and Vrms of the simulation
            tspan = 1/(2*f); % time for half cycle
            N = length(Vout);
            n = round(tspan/runtime*N); % number of datapoints in each cycle
            cycle = Vout((N-n):end); % last output cycle during simulation
            
            % Determine peaks from last output cycle
            peaks = zeros(10,1);
            tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
            for j = 1:size(tf,1)
                if tf(j) > 0
                    peaks(j) = tf(j);
                end
            end
            
            vrms = 0;
            numpeaks = 0;
            for t = 1:size(peaks,1)
                if peaks(t) > 0
                    vrms = vrms + (peaks(t))^2;
                    numpeaks = numpeaks + 1;
                end
            end
            X = sqrt(vrms/numpeaks);
            
            B(q) = -X; % Vrms for given parameters
            C(q) = -max(peaks); % Max Vout for given parameters
            D(2*q-q,i+2) = A(q,3);
            D(2*q-q,i+3) = X;
        end
    elseif i == 4
        for q = 1:k
            wt = 9e-6;
            wb = 30e-6;
            lf = 400e-6;
            m = A(q,4);
            tfd = 25e-6;
            g0 = 14e-6;
            tp = 120e-9;
            g = 1;
            f = 20;
            
            runtime = 0.3; % simulation time
            sim('model',runtime)
            % Determines the max Vout and Vrms of the simulation
            tspan = 1/(2*f); % time for half cycle
            N = length(Vout);
            n = round(tspan/runtime*N); % number of datapoints in each cycle
            cycle = Vout((N-n):end); % last output cycle during simulation
            
            % Determine peaks from last output cycle
            peaks = zeros(10,1);
            tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
            for j = 1:size(tf,1)
                if tf(j) > 0
                    peaks(j) = tf(j);
                end
            end
            
            vrms = 0;
            numpeaks = 0;
            for t = 1:size(peaks,1)
                if peaks(t) > 0
                    vrms = vrms + (peaks(t))^2;
                    numpeaks = numpeaks + 1;
                end
            end
            X = sqrt(vrms/numpeaks);
            
            B(q) = -X; % Vrms for given parameters
            C(q) = -max(peaks); % Max Vout for given parameters
            D(2*q-q,i+3) = A(q,4);
            D(2*q-q,i+4) = X;
        end
    elseif i == 5
        for q = 1:k
            wt = 9e-6;
            wb = 30e-6;
            lf = 400e-6;
            m = 2.0933e-6;
            tfd = A(q,5);
            g0 = 14e-6;
            tp = 120e-9;
            g = 1;
            f = 20;
            
            runtime = 0.3; % simulation time
            sim('model',runtime)
            % Determines the max Vout and Vrms of the simulation
            tspan = 1/(2*f); % time for half cycle
            N = length(Vout);
            n = round(tspan/runtime*N); % number of datapoints in each cycle
            cycle = Vout((N-n):end); % last output cycle during simulation
            
            % Determine peaks from last output cycle
            peaks = zeros(10,1);
            tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
            for j = 1:size(tf,1)
                if tf(j) > 0
                    peaks(j) = tf(j);
                end
            end
            
            vrms = 0;
            numpeaks = 0;
            for t = 1:size(peaks,1)
                if peaks(t) > 0
                    vrms = vrms + (peaks(t))^2;
                    numpeaks = numpeaks + 1;
                end
            end
            X = sqrt(vrms/numpeaks);
            
            B(q) = -X; % Vrms for given parameters
            C(q) = -max(peaks); % Max Vout for given parameters
            D(2*q-q,i+4) = A(q,5);
            D(2*q-q,i+5) = X;
        end
    elseif i == 6
        for q = 1:k
            wt = 9e-6;
            wb = 30e-6;
            lf = 400e-6;
            m = 2.0933e-6;
            tfd = 25e-6;
            g0 = A(q,6);
            tp = 120e-9;
            g = 1;
            f = 20;
            
            runtime = 0.3; % simulation time
            sim('model',runtime)
            % Determines the max Vout and Vrms of the simulation
            tspan = 1/(2*f); % time for half cycle
            N = length(Vout);
            n = round(tspan/runtime*N); % number of datapoints in each cycle
            cycle = Vout((N-n):end); % last output cycle during simulation
            
            % Determine peaks from last output cycle
            peaks = zeros(10,1);
            tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
            for j = 1:size(tf,1)
                if tf(j) > 0
                    peaks(j) = tf(j);
                end
            end
            
            vrms = 0;
            numpeaks = 0;
            for t = 1:size(peaks,1)
                if peaks(t) > 0
                    vrms = vrms + (peaks(t))^2;
                    numpeaks = numpeaks + 1;
                end
            end
            X = sqrt(vrms/numpeaks);
            
            B(q) = -X; % Vrms for given parameters
            C(q) = -max(peaks); % Max Vout for given parameters
            D(2*q-q,i+5) = A(q,6);
            D(2*q-q,i+6) = X;
        end
    elseif i == 7
        for q = 1:k
            wt = 9e-6;
            wb = 30e-6;
            lf = 400e-6;
            m = 2.0933e-6;
            tfd = 25e-6;
            g0 = 14e-6;
            tp = A(q,7);
            g = 1;
            f = 20;
            
            runtime = 0.3; % simulation time
            sim('model',runtime)
            % Determines the max Vout and Vrms of the simulation
            tspan = 1/(2*f); % time for half cycle
            N = length(Vout);
            n = round(tspan/runtime*N); % number of datapoints in each cycle
            cycle = Vout((N-n):end); % last output cycle during simulation
            
            % Determine peaks from last output cycle
            peaks = zeros(10,1);
            tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
            for j = 1:size(tf,1)
                if tf(j) > 0
                    peaks(j) = tf(j);
                end
            end
            
            vrms = 0;
            numpeaks = 0;
            for t = 1:size(peaks,1)
                if peaks(t) > 0
                    vrms = vrms + (peaks(t))^2;
                    numpeaks = numpeaks + 1;
                end
            end
            X = sqrt(vrms/numpeaks);
            
            B(q) = -X; % Vrms for given parameters
            C(q) = -max(peaks); % Max Vout for given parameters
            D(2*q-q,i+6) = A(q,7);
            D(2*q-q,i+7) = X;
        end
    else
        disp('Eigth set')
        for q = 1:k
            wt = 9e-6;
            wb = 30e-6;
            lf = 400e-6;
            m = 2.0933e-6;
            tfd = 25e-6;
            g0 = 14e-6;
            tp = 120e-9;
            g = A(q,8);
            f = 20;
            
            runtime = 0.3; % simulation time
            sim('model',runtime)
            % Determines the max Vout and Vrms of the simulation
            tspan = 1/(2*f); % time for half cycle
            N = length(Vout);
            n = round(tspan/runtime*N); % number of datapoints in each cycle
            cycle = Vout((N-n):end); % last output cycle during simulation
            
            % Determine peaks from last output cycle
            peaks = zeros(10,1);
            tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
            for j = 1:size(tf,1)
                if tf(j) > 0
                    peaks(j) = tf(j);
                end
            end
            
            vrms = 0;
            numpeaks = 0;
            for t = 1:size(peaks,1)
                if peaks(t) > 0
                    vrms = vrms + (peaks(t))^2;
                    numpeaks = numpeaks + 1;
                end
            end
            X = sqrt(vrms/numpeaks);
            
            B(q) = -X; % Vrms for given parameters
            C(q) = -max(peaks); % Max Vout for given parameters
            D(2*q-q,i+7) = A(q,8);
            D(2*q-q,i+8) = X;
        end
    
    end
%     for i = 1:k
%         wt = A(i,1);
%         wb = A(i,2);
%         lf = A(i,3);
%         m = A(i,4);
%         tfd = A(i,5);
%         lss = A(i,6);
%         g0 = A(i,7);
%         tp = A(i,8);
%         g = A(i,9);
%         f = A(i,10);
%         
%         runtime = 0.3; % simulation time
%         sim('model',runtime)
%         % Determines the max Vout and Vrms of the simulation
%         tspan = 1/(2*f); % time for half cycle
%         N = length(Vout);
%         n = round(tspan/runtime*N); % number of datapoints in each cycle
%         cycle = Vout((N-n):end); % last output cycle during simulation
%         
%         % Determine peaks from last output cycle
%         peaks = zeros(10,1);
%         tf = findpeaks(smooth(cycle)); % smooths output and finds peaks
%         for j = 1:size(tf,1)
%             if tf(j) > 0
%                 peaks(j) = tf(j);
%             end
%         end
%         
%         vrms = 0;
%         numpeaks = 0;
%         for k = 1:size(peaks,1)
%             if peaks(k) > 0
%                 vrms = vrms + (peaks(k))^2;
%                 numpeaks = numpeaks + 1;
%             end
%         end
%         X = sqrt(vrms/numpeaks);
%         
%         B(i) = -X; % Vrms for given parameters
%         C(i) = -max(peaks); % Max Vout for given parameters
%     end
end

% Save data to excel
filename = 'data2.xlsx';
writematrix(D,filename,'Sheet',1,'Range','A1:P250')

% B = -1*B;
% pos = 0;
% maxV = 0;
% m = size(B,1);
% for q = 1:m
%    if B(q) > maxV
%        maxV = B(q);
%        pos = pos + 1;
%    end
% end
% 
% wb = A(pos,1);
% wt = A(pos,2);
% Le = A(pos,3);
% disp(['wb = ', num2str(wb),' m'])
% disp(['wt = ', num2str(wt),' m'])
% disp(['Le = ', num2str(Le),' m'])
% 
% runtime = 0.3; % simulation time
% sim('opt_study',runtime)
% figure()
% plot(time,Vout)
% 
% figure()
% yyaxis left
% plot(Vout)
% xlabel('Iteration')
% ylabel('Vout (V)')
% title('Voltage vs. Shuttle Displacement')
% xlim([400000,452500])
% ylim([-.008 .01])
% hold on
% yyaxis right
% plot(y)
% ylabel('Displacement (m)')
% %ylim([-15,15])
% hold off
