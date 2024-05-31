function dy = sir_ComDemografia(t,y)
% ============== Declarando o valor de cada Par�metro ==================
% ===== Esses s�o os valores que ser�o considerados no modelo ==========
n=50000; %50 mil indiv�duos
beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
mu=1/(80*365); %(expectativa de 80 anos de vida)
%mu=0;
%d=0;
nu=1/60; %(a taxa � 1 dividido pelo tempo em que se fica infectado)

%alpha=0.001; %( mortalidade pela doen�a de 1 a cada 1000)
alpha=0;
%% ================ Equa��es do Modelo =================================
dy = zeros(3,1);

% %  ======= vers�o com demografia e mortalidade induzida pela doen�a ===
dy(1)=mu*n-beta*y(1)*y(2)-mu*y(1);              %Suscet�veis      
dy(2)=beta*y(1)*y(2)-(mu+nu+alpha)*y(2);        %Infectados   
dy(3)=nu*y(2)-mu*y(3);                          %Recuperados
% % ============== vers�o do modelo SIR sem demografia =================
% dy(1)=-beta*y(1)*y(2);              %Suscet�veis      
% dy(2)=beta*y(1)*y(2)-nu*y(2);       %Infectados   
% dy(3)=nu*y(2);                      %Recuperados
