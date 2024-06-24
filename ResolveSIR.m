%% Material Complementar do Livro: 
% Introdução à Epidemiologia Matemática: Métodos em Estudos Transversais

% = Outros Materiais estão disponíveis em https://linktr.ee/livroepidmat =

% ======== Programa para Simular o Modelo epidemiológico SIR ============

%comandos para não levar nenhum resíduo computacional da simulação anterior
clc;  clear all;  close all; 
% =========== Declarando os parâmetros (por dia) ========================
% Esses parâmetros são utilizados aqui caso calcule o valor de Ro ou das
% coordenadas do ponto de equilíbrio endêmico, por exemplo.

%% ======= Calculando Ro e o Ponto de Equilíbrio Endêmico ===============
%%============= Comente caso não queria calcular ========================
% Esses valores são utilizados somente para o cálculo de Ro, ou do ponto de equilíbrio (quando escrito). 
% Para simular o modelo, é necessário declarar o valor dentro da função "sir".  
n=50000; %50 mil indivíduos
beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
nu=1/60; %(a taxa é 1 dividido pelo tempo em que se fica infectado)

Ro=(beta*n)/nu
% Seq=n/Ro     
% Ieq=(Ro-1)*(mu/beta)   
% Ieq2=((mu*n)/(nu+mu))-(mu/beta)   
% Req=(Ro-1)*nu/beta
% 
% Discriminante=(mu*mu)*(Ro*Ro)-(4*mu*(mu+nu)*(Ro-1))

%% =======================================================================
% % =================  Resolvendo o Sistema de EDO ========================
t=600; %tempo em dias

%aplicando o pacote ode45 no arquivo função sir, que contém as equações
% 0 vetor são as condições iniciais de suscetíveis, infectados e
% recuperados, respectivamente

[T,Y]=ode45('sir',[0 t],[49999 1 0]);
% Aqui T é o tempo em dias e Y a matriz com a solução do modelo. O pacote
% ode45 está considerando a função "sir", que contém as equações do modelo,
% no intervalo de tempo inicial "0" e tempo final "t" (que é o valor
% definido acima). Já nos últimos argumentos, temos como "49999 a
% quantidade de indivíduos Suscetíveis, 1 Infectado e 0 Recuperados. Note
% que a soma da quantidade de indivíduos é exatamente igual ao N (a
% população total considerada no modelo).

%% ================ Plotando as Soluções ====================
figure(1) %plotando as 3 soluções do modelo no mesmo gráfico
plot(T,Y(:,1),'k',T,Y(:,2),'--k',T,Y(:,3),'-.k')
xlabel('Tempo (dias)'),
ylabel('População Total'),
legend('Suscetíveis', 'Infectados', 'Recuperados')

figure(2) %plotando somente a solução de infectados do modelo
plot(T,Y(:,2),'--k')
xlabel('Tempo (dias)'),
ylabel('População Total de Infectados'),
legend( 'Infectados')

