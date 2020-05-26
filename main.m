
[xdata,textdata]=xlsread('11.xls'); %����10�����е����ݣ����ݰ��ձ���е�λ�ñ�����Excel�ļ�11.xls��
x_label=xdata(:,2); %�ڶ���Ϊ������
y_label=xdata(:,3); %������Ϊ������
C=[x_label y_label];      %�������
n=size(C,1);  %n��ʾ���и���
D=zeros(n,n); %D��ʾ��ȫͼ�ĸ�Ȩ�ڽӾ��󣬼��������D��ʼ��
for i=1:n
   for j=1:n
       if i~=j
           D(i,j)=((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5; %����������֮��ľ���
       else
           D(i,j)=0;   %i=j, �����Ϊ0��
       end
   end
end
 %%==================��Ⱥ�㷨ʵ�ֹ���======================================================
 %%============== ��һ�� ������ʼ��==============
iter_max=100;   %����������
m=30;           % ���ϸ���
Alpha=1;        % ������Ϣ����Ҫ�̶ȵĲ���
Beta=5;         % ��������ʽ������Ҫ�̶ȵĲ���
Rho=0.8;        % ��Ϣ������ϵ��
Q=10;           % ��Ϣ������ǿ��ϵ��
Eta=1./D;          % EtaΪ�ܼ���������������Ϊ����ĵ���
Tau=ones(n,n);     % TauΪ��Ϣ�ؾ��󣬳�ʼ��ȫΪ1
Tabu=zeros(m,n);   % �洢����¼·��������
nC=1;              % ����������
R_best=zeros(iter_max,n);   %�������·�ߣ���Ϊ��������������Ϊ���и���
L_best=inf.*ones(iter_max,1);%%�������·�ߵĳ��ȣ�infΪ�����
L_ave=zeros(iter_max,1);     % ����ƽ��·�߳���
 
 %%============== �ڶ��� ��mֻ���Ϸŵ�������==============
while nC<=iter_max    %ֹͣ����֮һ���ﵽ���������� 
    Randpos=[];
    for i=1:(ceil(m/n))       %ceil��ʾ�������ȡ��
        Randpos=[Randpos,randperm(n)]; %randperm(n)����ʾ�������һ����������
    end
 Tabu(:,1)=(Randpos(1,1:m))'; %ÿֻ���ϣ�mֻ������Ӧ��һ��λ�ã�Tabu(:,1)Ϊÿֻ�����߹��ĵ�һ������
  
%% ============== ������ mֻ���ϰ����ʺ���ѡ����һ�����У���ɸ��Ե�����==============
  for j=2:n       %���дӵڶ�����ʼ
     for  i=1:m
        visited=Tabu(i,1:(j-1));      %�ѷ��ʵĳ���
        J=zeros(1,(n-j+1));           %�����ʵĳ���
        P=J;                          %�����ʳ��е�ѡ����ʷֲ�����ʼ����
        Jc=1;                         %ѭ���±�
            
       for k=1:n     %����ѭ���������ʳ��У������k�����в������ѷ��ʳ��У�����Ϊ�����ʳ���
          if  length(find(visited==k))==0
            J(Jc)=k;
            Jc=Jc+1;   %�±��1��������һ���洢�����ʵĳ���
          end
       end
      
       for k=1:length(J)   % �����������ʳ��еĸ��ʷֲ���length(J)��ʾ�����ʳ��и���
         P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta); %���ʼ��㹫ʽ�еķ���
       end
         P=P/(sum(P));     %���ʷֲ�������Ϊ�����ʳ��и���
         Pcum=cumsum(P);   %���ۻ����ʺͣ�cumsum��[1 2 3])=1 3 6,Ŀ������ʹ��Pcum��ֵ���д���rand����
         Select=find(Pcum>=rand);  %������ѡȡ��һ�����У����ۻ����ʺʹ��ڸ��������������ѡ����ͱ����ϵ����һ��������Ϊ�������ʵĳ���
       if  isempty(Select)    %��ѡ�����Ϊ�ռ������������һ���м�����ɱ���
         Tabu(i,j)=round(1+(n-1)*rand);
       else
         next_visit=J(Select(1));   %next_visit��ʾ�������ʵĳ���
         Tabu(i,j)=next_visit;      %�����ʹ��ĳ��м�����ɱ���
       end
     end
  end
    
    if nC>=2;Tabu(1,:)=R_best(nC-1,:);end  %�������������ڵ���2������һ�ε��������·�ߴ��뵽Tabu�ĵ�һ����
 
%% ==============���Ĳ� ��¼���ε������·��==============
 L=zeros(m,1);
  for i=1:m;
      R=Tabu(i,:);
    for j=1:(n-1)
      L(i)=L(i)+D(R(j),R(j+1));  %��·������
    end
      L(i)=L(i)+D(R(1),R(n));    %�������һ���������һ������֮��ľ���
  end
  L_best(nC)=min(L);            %����·��Ϊ������̵�·��
  pos=find(L==L_best(nC));      %�ҳ�����·����Ӧ��λ�ã���Ϊ��ֻ����
  R_best(nC,:)=Tabu(pos(1),:);  %ȷ������·����Ӧ�ĳ���˳��
  L_ave(nC)=mean(L);            %���k�ε�����ƽ������
  nC=nC+1;
   
%% ==============���岽 ������Ϣ�أ��˴�����ϵͳ==============
 Delta_Tau=zeros(n,n);  %Delta_Tau(i,j)��ʾ�����������ڵ�i�����е���j������·���ϵ���Ϣ������
   for i=1:m
      for j=1:(n-1)     %����������·�������ͷ���Ϣ��
        Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
      end
        Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
   end
Tau=(1-Rho).*Tau+Delta_Tau;   %��Ϣ�ظ��¹�ʽ
 
%% ==============������ ���ɱ�����==============
Tabu=zeros(m,n);
end
 
%% ==============���߲� ������==============
Pos=find(L_best==min(L_best));     %�ҵ�L_best����Сֵ���ڵ�λ��
Shortest_Route=R_best(Pos(1),:)   %��ȡ���·��
Shortest_Length=L_best(Pos(1))    %��ȡ���·������
 
%% ==============��ͼ==============
figure(1)   %��������������ͼ
x=linspace(0,iter_max,iter_max);
y=L_best(:,1);
plot(x,y,'-','LineWidth',2);
xlabel('��������'); ylabel('���·������');
 
figure(2)   %�����·��ͼ
Shortest_Route=[Shortest_Route Shortest_Route(1)];
plot([C(Shortest_Route,1)],[C(Shortest_Route,2)],'o-');
grid on
for i = 1:size(C,1)
    text(C(i,1),C(i,2),['   ' num2str(i)]);
end
xlabel('���к�����'); ylabel('����������'); 