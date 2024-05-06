MonteCarlo = 1e2;%仿真次数，蒙特卡洛
K = 2000;%活跃设备数
N = 40000;%总设备
M_set = [2,4,8,16,32,64,128];%天线数
L = 300;%导频长度？


N_md = K*MonteCarlo*ones(length(M_set),1);%活跃设备对应的实验总数
N_fa = (N-K)*MonteCarlo*ones(length(M_set),1);失效设备对应的实验总数
P_md = zeros(length(M_set),1);%每种天线数的漏检概率
P_fa = zeros(length(M_set),1);%每种天线数的错检概率

for i = 1:length(M_set)
    display(strcat('M_idx=',num2str(i)));
    M = M_set(i);%索引M
    D_channel = zeros(N,MonteCarlo);%N×Mon全零矩阵，每个mon实验中的信道
    D_act = false(N,MonteCarlo);%N×Mon全逻辑零矩阵，激活？
    D_signal = zeros(N,M,MonteCarlo);%三维数组，每个蒙特卡罗模拟中的每个设备和每个天线配置的接收信号强度？
    D_absxmm = zeros(N,MonteCarlo);%该矩阵可以表示对接收信号的一些绝对值运算
    D_tau = zeros(MonteCarlo,1);
    imc = 1;%蒙特卡洛的索引
    while imc <= MonteCarlo
        display(strcat('Mc_idx=',num2str(imc)));
        
        A = (randn(L,N) + sqrt(-1)*randn(L,N))*sqrt(1/2.*1/L);%L×N矩阵，表示导频序列，零均值、1/L方差复高斯序列
        supp = randperm(N);%N个设备的随机排列
        D_act(supp(1:K),imc) = true;把设备随机排列后，选择前K个设备作为活跃设备
      
        distance = zeros(N,1);
        l_user = ((rand(N,2) - 0.5*ones(N,2))*2*100*sqrt(2) + 400*sqrt(2)).*randsrc(N,2);%N×2，矩阵每行作为设备的坐标（x，y）
        for n = 1:N
            distance(n) = sqrt((l_user(n,1))^2 + (l_user(n,2))^2);%每个设备到参考点（基站）的矩阵
        end%距离和原文不一样，原文设置的是[0.05km,1km]，这个设置的是[0.6km,1km]
        

            

        path_loss = zeros(N,1);
        for n=1:N
            path_loss(n) = 10^((-128.1 - 37.6*log10(distance(n)*1e-3))/10);%路径损耗，单位是mW？
        end
        
        power = 10^(1.3)*10^(-3);%发射功率，单位是W，论文里面是23dbm
        noise_power = 10^(-16.9)*10^(-3);
        B = 1e7;信道带宽
        noise = noise_power*B;

        
        D_channel(:,imc) = path_loss;

        h=zeros(N,M);
        for n=1:N
            h(n,:) = sqrt(1/2)*(randn(1,M) + sqrt(-1)*randn(1,M))*sqrt(path_loss(n));%信道建模，CN(0,path_lossI)
        end


        x = zeros(N,M);
        for m=1:M
            x(supp(1:K),:) = h(supp(1:K),:);
        end
        
        w = (randn(L,M) + sqrt(-1)*randn(L,M))*sqrt(1/2);%噪声建模
        

        
        noise_r = noise/power/L;
        sigma_w = sqrt(noise_r);
        

        
        y = A*x + w*sigma_w;
        
       [xnoise,xhat,mse,tau_real,tau_est] = noisyCAMPmmseforKLS(A,N,M,L,y,x,50,K/N,path_loss,sigma_w);
       
        D_tau(imc) = tau_est(end);
        
        D_signal(:,:,imc) = xnoise;
        for n = 1:N
            D_absxmm(n,imc) = norm(D_signal(n,:,imc));
        end
        imc = imc + 1;
        
       

    end
    for imc = 1:MonteCarlo 
        for n = 1:N
            t = M*log(1+path_loss(n)/D_tau(imc)^2)/(1/D_tau(imc)^2-1/(path_loss(n)+D_tau(imc)^2));
            if D_absxmm(n,imc)^2 >= t && D_act(n,imc) == 1
                N_md(i) = N_md(i) - 1;
            end
            if D_absxmm(n,imc)^2 <= t && D_act(n,imc) == 0
                N_fa(i) = N_fa(i) - 1;
            end
        end
    end
    

    P_md(i) = N_md(i)/(K*MonteCarlo);
    P_fa(i) = N_fa(i)/((N-K)*MonteCarlo);

    
end


M = M_set;
figure
semilogy(M,P_md,'ko-',M,P_fa,'k+-');
grid on
xlabel('Number of BS Antennas: {\it M}');
ylabel('Activity Detection Error Probability');
legend('{\it P}^{MD}: MMSE Denoiser','{\it P}^{FA}: MMSE Denoiser');

