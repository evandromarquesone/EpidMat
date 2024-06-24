% Material Complementar do Livro: 
% Introdução à Epidemiologia Matemática: Métodos em Estudos Transversais

% = Outros Materiais estão disponíveis em https://linktr.ee/livroepidmat =

function dy = sir_ComDemografia(t,y)
% ============== Declarando o valor de cada Parâmetro ==================
% ===== Esses são os valores que serão considerados no modelo ==========
n=50000; %50 mil indivíduos
beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
mu=1/(80*365); %(expectativa de 80 anos de vida)
%mu=0;
%d=0;
nu=1/60; %(a taxa é 1 dividido pelo tempo em que se fica infectado)

%alpha=0.001; %( mortalidade pela doença de 1 a cada 1000)
alpha=0;
%% ================ Equações do Modelo =================================
dy = zeros(3,1);

% %  ======= versão com demografia e mortalidade induzida pela doença ===
dy(1)=mu*n-beta*y(1)*y(2)-mu*y(1);              %Suscetíveis      
dy(2)=beta*y(1)*y(2)-(mu+nu+alpha)*y(2);        %Infectados   
dy(3)=nu*y(2)-mu*y(3);                          %Recuperados
% % ============== versão do modelo SIR sem demografia =================
% dy(1)=-beta*y(1)*y(2);              %Suscetíveis      
% dy(2)=beta*y(1)*y(2)-nu*y(2);       %Infectados   
% dy(3)=nu*y(2);                      %Recuperados
