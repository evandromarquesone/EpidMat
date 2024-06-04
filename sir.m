% Material Complementar do Livro: 
% T�picos Matem�ticos Aplicados na Modelagem em Epidemiologia - Estudos Transversais

% = Outros Materiais est�o dispon�veis em https://linktr.ee/livroepidmat =

function dy = sir(t,y)
% ============== Declarando o valor de cada Par�metro ==================
% Esses ser�o os valores utilizados para computar as colu��es num�ricas de cada popula��o do modelo 
n=50000; %50 mil indiv�duos
beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
nu=1/60; %(a taxa � 1 dividido pelo tempo em que se fica infectado)
%% ================ Equa��es do Modelo =================================
dy = zeros(3,1); % declarando um vetor nulo, para evitar que o vetor criado na simula��o seja n�o nulo.

dy(1)=-beta*y(1)*y(2);              %Suscet�veis      
dy(2)=beta*y(1)*y(2)-nu*y(2);       %Infectados   
dy(3)=nu*y(2);                      %Recuperados