clear
% close all

%% Enter Inputs
name_EC = 'b03.xlsx'; %name of pdb file
N = 360; %number of residues
m1 = 1; %first mode
m2 = 10; %last mode (make m1 = m2 to evaluate the individual modes)
tau_0 = 6;
threshold = 0.80;

%% Prepare Kirchhoff matrix
EC = xlsread(name_EC);

K = zeros(N,N); %Kirchhoff matrix
for i = 1:length(EC) %Coupled contacts
    if EC(i,4) > threshold
        K(EC(i,1),EC(i,2)) = -EC(i,3);
        K(EC(i,2),EC(i,1)) = -EC(i,3);
    end
end
for i = 1:N %Chain connectivity
    for j = 1:N
        if i ~= j && abs(i-j) < 4
            K(i,j) = K(i,j) -1;
        end
    end
end
K(1:size(K,1)+1:end) = -sum(K); %diagonal elements

%% Matrix Decomposition 
[U,S,V] = svd(K); %standart value decomposition of K
S = diag(S);
S = flipud(S); %eigenvalues in ascending order
U = fliplr(U); %eigenvectors

%% Find Maximum Information Tau
lower = 0.01; %lower bound for tau
inc = 20; %incremantation of tau
upper = 4000; %upper bound for tau
Tau = zeros(N,N); %max information tau for each residue pair
for i = [73,75:76,79,125,148,149,161,193,197,246,274,293,296] %number of pair can be adjusted here
    for j = 5:10:N
        index = 0; %index of Tij
        Tij = zeros(1,length(lower:inc:upper)); %information transfer for every tau
        for tau = lower:inc:upper
            index = index + 1;
            sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0; sum5 = 0; 
            %Sum number matchs with order of appearance of distinct
            %summation terms in equation 14 of Supporting Information text
            %of https://doi.org/10.1002/prot.25272
            for k = m1+1:m2+1
                sum1 = sum1 + U(j,k) * (1 / S(k)) * U(j,k);
                sum2 = sum2 + U(j,k) * (1 / S(k)) * U(j,k) * exp(-S(k) * tau / tau_0);
                sum3 = sum3 + U(i,k) * (1 / S(k)) * U(i,k);
                sum4 = sum4 + U(i,k) * (1 / S(k)) * U(j,k);
                sum5 = sum5 + U(i,k) * (1 / S(k)) * U(j,k) * exp(-S(k) * tau / tau_0);
            end
            
            %Term numbers match with order of appearrance of distinct
            %natural logarithm terms in equation 14 in the same text.
            term1 = sum1^2 - sum2^2;
            term2 = sum3 * sum1^2 + 2 * sum4 * sum2 * sum5 - (sum5^2 + sum4^2) * sum1 - sum2^2 * sum3;
            term3 = sum1;
            term4 = sum3 * sum1 - sum4^2;
            
            %Necessary to avoid compuational crush
            if abs(term1)<9*10^-4 && term1<0
                term1 = abs(term1);
            end
            if abs(term2)<9*10^-4 && term2<0
                term2 = abs(term2);
            end
            if abs(term3)<9*10^-4 && term3<0
                term3 = abs(term3);
            end
            if abs(term4)<9*10^-4 && term4<0
                term4 = abs(term4);
            end
            
            %Final forms of each term in equation 14
            term1 = 0.5 * reallog(term1);
            term2 = 0.5 * reallog(term2);
            term3 = 0.5 * reallog(term3);
            term4 = 0.5 * reallog(term4);
            
            Tij(index) = term1 - term2 - term3 + term4; %equation 14
        end

        [pks, locs] = findpeaks(Tij,lower:inc:upper); %Detect tau giving max information transfer
        if isempty(locs) == 0
            Tau(i,j) = locs(1);
        end
    end
end

taus = nonzeros(Tau);
tau = mean(taus)*3; %Value of time delay to be used in rest of the calculation
%% Calculate Transfer Entropy
% Same calculation as previous section is run for determined tau.
Tij = zeros(N,N);
for i = 1:N
    for j = 1:N
        if i == j
            Tij(i,j) = 0;
        else
            sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0; sum5 = 0; 

            for k = m1+1:m2+1
                sum1 = sum1 + U(j,k) * (1 / S(k)) * U(j,k);
                sum2 = sum2 + U(j,k) * (1 / S(k)) * U(j,k) * exp(-S(k) * tau / tau_0);
                sum3 = sum3 + U(i,k) * (1 / S(k)) * U(i,k);
                sum4 = sum4 + U(i,k) * (1 / S(k)) * U(j,k);
                sum5 = sum5 + U(i,k) * (1 / S(k)) * U(j,k) * exp(-S(k) * tau / tau_0);
            end

            term1 = sum1^2 - sum2^2;
            term2 = sum3 * sum1^2 + 2 * sum4 * sum2 * sum5 - (sum5^2 + sum4^2) * sum1 - sum2^2 * sum3;
            term3 = sum1;
            term4 = sum3 * sum1 - sum4^2;

            if abs(term1)<9*10^-4 && term1<0
                term1 = abs(term1);
            end
            if abs(term2)<9*10^-4 && term2<0
                term2 = abs(term2);
            end
            if abs(term3)<9*10^-4 && term3<0
                term3 = abs(term3);
            end
            if abs(term4)<9*10^-4 && term4<0
                term4 = abs(term4);
            end

            term1 = 0.5 * reallog(term1);
            term2 = 0.5 * reallog(term2);
            term3 = 0.5 * reallog(term3);
            term4 = 0.5 * reallog(term4);

            Tij(i,j) = term1 - term2 - term3 + term4;
        end
    end
end

%% Calculate Net Transfer Entropy
netTE = zeros(N,N); %Net transfer entropy
for i = 1:N
    for j = 1:N
        netTE(i,j) = Tij(i,j) - Tij(j,i);
    end
end

%% Figure
figure;
imagesc(netTE(1:length(netTE),1:length(netTE)));
set(gca,'YDir','normal');
set(gca, 'FontSize',16);
axis square;
colorbar;
clim([-0.01 0.01])
colormap(jet)
xlabel('Residue Index (Affected)')
ylabel('Residue Index (Effector)')
if m1 ~= m2
    name = strcat(num2str(m1),'-',num2str(m2),32,'Modes - Net Transfer Entropy');
else
    name = strcat('Mode',32,num2str(m1),32,'- Net Transfer Entropy.fig');
end
title(name)

% saveas(gcf,strcat(name,'.fig'))
% saveas(gcf,strcat(name,'.png'))

%% Find Collectivity
net_pos = max(netTE, 0);
net_sq = net_pos.^2;
sum_sq = sum(net_sq, 2); %sum along columns for each residue

alfa1 = 1 ./ sum_sq;

P = net_sq .* alfa1;
P_logP = zeros(size(P));
pos_mask = P > 0;
P_logP(pos_mask) = P(pos_mask) .* log(P(pos_mask));

insidesum = sum(P_logP, 2)';
col = exp(-insidesum) / N; %collectivity

