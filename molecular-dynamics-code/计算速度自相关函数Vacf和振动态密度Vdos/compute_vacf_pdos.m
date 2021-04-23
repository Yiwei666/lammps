v_all=load('1dump-file2');   %加载速度文件
N=1400;                      %原子个数
dt=0.001;                    %ps
Nc=500;                      %关联点数
Nt=500;
Ns=size(v_all,1)/N;          %帧数

for i=1:Ns
    D(:,:,i)=v_all((i-1)*N+1:i*N,:);
end
M=Ns-Nc;               % number of time origins for time average 

%function [vacf]=find_vacf(v_all,M,Nt)
% v_all(:,:,:): all the velocity data
% M: number of "time origins" used in "time average"
% Nt: number of correlation time points
% vacf(:): VACF data
vacf=zeros(Nc,1);
for nt=0:Nc-1
    for m=1:M
        vacf(nt+1,:)=vacf(nt+1,:)+sum(sum(D(:,:,m+0).*D(:,:,m+nt)));
    end
end
vacf=vacf/vacf(1);
vacf_use=vacf;     

figure(1);
t=(0:Nc-1)*dt;
plot(t,vacf_use);
xlabel('Time (ps)','fontsize',20);
ylabel('VACF (normalized)','fontsize',20);
set(gca,'fontsize',20,'linewidth',1.5);
set(gca,'ticklength',get(gca,'ticklength')*2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vacf=load('vacf4.txt');
i=1;                                  %选取vacf.txt文件的第i列作为输入
omega=0:0.05:90*pi;         %ps-1=1Thz 角频率
%0.05:0.05:20

vacf=vacf(:,i);
%Nt=500;
delta_t=0.001;        %ps

%function [pdos]=find_pdos(vacf,omega,Nt,delta_t)
% vacf(:): VACF
% omega: phonon angular frequency in units of ps^{-1}
% Nt: number of correlation time points
% delta_t: time interval between two measurements, in units of ps
% pdos(:): PDOS
vacf=vacf.'.*(cos(pi*(0:Nt-1)/Nt)+1)*0.5; % apply a window function
vacf=vacf.*[1,2*ones(1,Nt-1)]/pi; % C(t)=C(-t) and there is only one C(0)
pdos=zeros(length(omega),1); % phonon density of states
for n=1:length(omega) % Discrete cosine transformation
    pdos(n)=delta_t*sum(vacf.*cos(omega(n)*(0:Nt-1)*delta_t));
end

%作图

figure(2)
plot(omega,pdos);

xlabel('\omega (1/ps)','fontsize',20);
ylabel('PDOS (ps)','fontsize',20);
%ylim([0,0.01]);
set(gca,'fontsize',20,'linewidth',1.5);
set(gca,'ticklength',get(gca,'ticklength')*2);



% An important check is that PDOS should be normalized to about 1
normalization_of_pdos=trapz(omega,mean(pdos,2))
omega=omega'/(3*pi)*50;     %单位cm-1



