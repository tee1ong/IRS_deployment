clc
clear

% Signal Parameters 
m = 4; %input('Enter Modulation Order: '); % Modulation order - VARIABLE!!!
k = log2(m); % Bits per Symbol
snr_vec = 0:30; % SNR vector - !
fc = 2.4e9; % Carrier frequency [Hz] - VARIABLE!!
lambda = 3e8/fc; % Wavelength
k_w = (2*pi)/lambda; % Wavenumber
n_sym = 3000; % number of test symbols - !
ts = 1e-5; % Symbol period - VARIABLE!
r_sym = 1/ts; % Symbol rate
b = r_sym*k; % Signal bandwidth

% Ask for number of users
n_users = 2;% input('Enter number of users: ');

% Antenna Layout and Parameters
p_in = 1; % Input Power [W] - !
n_t = 4; %input('Enter Number of Transmitter Antenna elements: '); % Number of transmitter antennas
n_r = 4; %input('Enter Number of Receiver Antenna elements: '); % Number of receiver antennas
da = 0.5*lambda; % Antenna spacing

antpuser = n_r/n_users;
for ntr = 1:n_t
    z_tr = (ntr-1)*da;
    tx_pos(ntr,:) = [30 10 5+z_tr];
end
for user_no = 1:n_users

    for nre = 1:antpuser
        z_re = (nre-1)*da;
        rx_pos{user_no}(nre,:) = [-30 -10 -5+z_re];
    end
end

%RIS Dimensions
n_rows = 4; 
n_columns = 4;

n_ris = n_rows*n_columns;
dris = 0.5*lambda;
% RIS Position
ris_pos = zeros(n_ris,3); % Cause 3D
idx = 1;
for row = 1:n_rows 
    for col = 1:n_columns
        y = (col-1)*dris;
        z = (row-1)*dris;
        ris_pos(idx,:) = [0,y,z];
        idx = idx + 1;
    end
end

for user_no = 1:n_users
    for ris_no = 1:n_ris
        for rx_no = 1:n_r
            dist_risue{user_no}(ris_no,rx_no) = norm((rx_pos{user_no}(n_r)-ris_pos(ris_no,:))); % Euclidean distance between the user andthe n-th RIS element 
            for rxchan = 1:size(rx_pos{user_no},1)
            v{user_no}(ris_no,rx_no) = acos(rx_pos{user_no}(rxchan,1)/dist_risue{user_no}(ris_no,rx_no)); % Angle of departure
            end
        end
        for tx_no = 1:n_t
            dist_txris(ris_no,tx_no) = norm(ris_pos(ris_no,:)-tx_pos(n_r));
        end
        dist_ris(ris_no) = norm((ris_pos(ris_no,:)-(-0.1)));
        z_imp(ris_no) = (sinc(k_w*(dist_ris(ris_no)-dist_ris(1))));
    end
end

% Path loss Calculation
for user_no = 1:n_users
    dist_tot{user_no} = dist_txris+dist_risue{user_no};
    p_l{user_no} = p_in*(((lambda/4)^2)./(4*pi*dist_tot{user_no}.^2));
end

%init
len = length(snr_vec);
ber_mu_ris = zeros(1,len);

for snr_value = 1:len
    snr = snr_vec(snr_value);
    sigma = sqrt((1)/(2*(10^(snr/10))));

    for i = 1:n_sym
        tx_data = randi([0 m-1],n_r,1);
        tx_sig = pskmod(tx_data,m);

        for user_no = 1:n_users
           % H{user_no} = sqrt(0.5)*randn(n_r/antpuser,n_t)+1j*randn(n_r/antpuser,n_t);
           H{user_no} = sqrt(0.5)*channel_ris(n_r/antpuser,n_t,n_ris);
        end
        [U1,S1,V1] = svd(H{1});
        W2 = V1(:,3:4);
        [U2,S2,V2] = svd(H{2});
        W1 = V2(:,3:4);

        n = (sqrt(0.5/(10^(snr/10)))) * (randn(min(n_r,n_t),1)+1j*randn(min(n_t,n_r),1)); % This is not

        tdata = W1*tx_sig(1:2)+W2*tx_sig(3:4);
        Rx1 = H{1}*tdata;
        Rx2 = H{2}*tdata;

        % Compute equalizers using block diagonalizing matrices
        W1_H1 = H{1}*W1;
        EQ1 = W1_H1'*inv(W1_H1*W1_H1');

        W2_H2 = H{2}*W2;
        EQ2 = W2_H2'*inv(W2_H2*W2_H2');

        % Equalize receivrd data and estimate symbols
        y = [EQ1*Rx1;EQ2*Rx2] + n;

        rx_data = pskdemod(y,m);
        [numerr errate] = biterr(rx_data,tx_data);
        
        ber_mu_ris(1,snr_value) = ber_mu_ris(1,snr_value)+errate;

%         n = (sqrt(0.5/(10^(snr/10))))*(randn(min(n_r,n_t),1)+1j*randn(min(n_t,n_r),1));

    end
end
ber_mu_ris = ber_mu_ris/n_sym;