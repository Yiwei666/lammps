% parameters setting of MD simulation
a=29.256558;             
b=29.256558;
c=29.256558;   
type2=1;                                 %配位原子为data文件里编号为5的原子即O原子
type_1=[4,3];                            %中心原子为data文件里编号为1,2的原子即Si,Al原子
N_atoms=1878;                            %原子数
rcut=[2.12,2.45];                          %依次为Si-O,Al-0的rcut

box=[a,0,0;0,b,0;0,0,c];                 % box matrix
A=importdata('Mg15.txt');                        %加载多帧的坐标文件r0.txt
pbc=[1,1,1].';                           % boundary conditions
Ns=size(A,1)/N_atoms;                    % 帧数
Ns=single(Ns);
A1=ones(1,Ns)*N_atoms;                   
B=mat2cell(A,A1,4);                      %分割矩阵 分成Ns份 
w_total=zeros(N_atoms,N_atoms);
for i=1:Ns
    B1=cell2mat(B(i));
    type=B1(:,1).';
    w=zeros(N_atoms,N_atoms);
    for j=1:2
        type1=type_1(j);       
        w=w+find_matrix(type1,type2,type,B1,box,pbc,rcut(j));
    end   
   	w_total=w_total+w;
end
w_ave=w_total/Ns;
w_ave=round(w_ave);                       %强制取整转化,极其重要的一点               

%统计不同种类氧的百分比
%桥氧 非桥氧 自由氧 氧簇
w_hang=sum(w_ave);
n_count=[0,0,0,0,0];                      %初始化计数1位为自由氧，2位非桥氧，3位桥氧，4位氧簇即三键氧,5位其他氧
position_BO=zeros(1,N_atoms);             %用来记录桥氧的位置
w_coum=sum(w_ave,2);
for k=1:N_atoms
    if type(k)~=type2
        continue;
    end
    if w_hang(k)==0                      %自由氧
        n_count(1)=n_count(1)+1;
		continue;
    end
    if w_hang(k)==1                      %非桥氧
        n_count(2)=n_count(2)+1;
		continue;
    end
    if w_hang(k)==2                      %桥氧
        position_BO(k)=position_BO(k)+1; %记录下桥氧位置
        n_count(3)=n_count(3)+1;
		continue;
    end
    if w_hang(k)==3                      %三键氧
        position_BO(k)=position_BO(k)+1; %将三键氧当做桥氧计算Qn
        n_count(4)=n_count(4)+1;
		continue;
    end
    n_count(5)=n_count(5)+1;             %其他种类的氧
end
n1_count=n_count(1:4);
O_percent=100*n1_count/sum(n1_count);    %四种氧的百分比例
fprintf('四种氧的百分比例%f\n',O_percent);               

%统计Qi的百分比    Q为Si或者Al原子
BO_lie1=zeros(N_atoms,1);                 %列初值化
[mi,ni]=find(position_BO==1);
mi=mi+1111111;                            %使用mi避免程序报错
for ibo=ni
    BO_lie1=BO_lie1+w_ave(:,ibo);         %桥氧列求和
end                                        
%统计Si或Al原子链接桥氧个数分别为0,1,2,3,4个
n_Si=[0,0,0,0,0];                    
n_Al=[0,0,0,0,0];

for kk=1:N_atoms
    if ~ismember(type(kk),type_1)                 
        continue;
    end
    if type(kk)==type_1(1)                        %当粒子为Si原子的时候
        if BO_lie1(kk)==0                 %0个桥氧
            n_Si(1)=n_Si(1)+1;
			continue;
        end
        if BO_lie1(kk)==1                 %1个桥氧
            n_Si(2)=n_Si(2)+1;
			continue;
        end
        if BO_lie1(kk)==2                 %2个桥氧
            n_Si(3)=n_Si(3)+1;
			continue;
        end
        if BO_lie1(kk)==3                 %3个桥氧                 
            n_Si(4)=n_Si(4)+1;
			continue;
        end
        if BO_lie1(kk)==4                 %4个桥氧                 
            n_Si(5)=n_Si(5)+1;
			continue;
        end
    end
    %当粒子为Al的时候,应该注意到Al5 A16的存在，并且含量不可忽略
	if type(kk)==type_1(2) && (w_coum(kk)==4)                       
		if BO_lie1(kk)==0                 %0个桥氧
			n_Al(1)=n_Al(1)+1;
			continue;
		end
		if BO_lie1(kk)==1                 %1个桥氧
			n_Al(2)=n_Al(2)+1;
			continue;
		end
		if BO_lie1(kk)==2                 %2个桥氧
			n_Al(3)=n_Al(3)+1;
			continue;
		end
		if BO_lie1(kk)==3                 %3个桥氧                 
			n_Al(4)=n_Al(4)+1;
			continue;
		end
		if BO_lie1(kk)==4                 %4个桥氧                 
			n_Al(5)=n_Al(5)+1;
			continue;
		end		
  end
end
fprintf('不同Sin个数输出%f\n',n_Si);                 %不同Qi个数输出
fprintf('不同Aln个数输出%f\n',n_Al);                 %个数Qi个数输出

%2019  体系Si-Al-Ca-0  txt文件数据格式：原子类型 x y z  %多帧求时间平均  
function[d]=find_matrix(type1,type2,type,r,box,pbc,rcut)
% determine some parameters      type is a vector pbc=[1,1,1] 
N=size(r,1);                     % number of particles
d=zeros(N,N);                    %creat a N*N matrix d
for n1=1:N
    if type(n1)~=type1           % type1 is the center atom
        continue;
    end
    for n2=1:N
        if type(n2)~=type2 || n1==n2     % type2 is the other atom
            continue;
        end
        r12=r(n2,2:4)-r(n1,2:4);
        r12=r12.'; 
        r12=box\r12;
        r12=r12-pbc.*round(r12); 
        r12=box*r12; 
        d12=sqrt(sum(r12.*r12));           % distance
        if d12>rcut                        
            continue;
        end
        d(n1,n2)=d(n1,n2)+1;
    end
end
end