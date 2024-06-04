% Material Complementar do Livro: 
% T�picos Matem�ticos Aplicados na Modelagem em Epidemiologia - Estudos Transversais

% = Outros Materiais est�o dispon�veis em https://linktr.ee/livroepidmat =

% ======== Programa para Simular o Modelo epidemiol�gico SIR ============

%comandos para n�o levar nenhum res�duo computacional da simula��o anterior
clc;  clear all;  close all; 
% =========== Declarando os par�metros (por dia) ========================
% Esses par�metros s�o utilizados aqui caso calcule o valor de Ro ou das
% coordenadas do ponto de equil�brio end�mico, por exemplo.

%% ======= Calculando Ro e o Ponto de Equil�brio End�mico ===============
%%============= Comente caso n�o queria calcular ========================
% Esses valores s�o utilizados somente para o c�lculo de Ro, ou do ponto de equil�brio (quando escrito). 
% Para simular o modelo, � necess�rio declarar o valor dentro da fun��o "sir".  
n=50000; %50 mil indiv�duos
beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
nu=1/60; %(a taxa � 1 dividido pelo tempo em que se fica infectado)

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

%aplicando o pacote ode45 no arquivo fun��o sir, que cont�m as equa��es
% 0 vetor s�o as condi��es iniciais de suscet�veis, infectados e
% recuperados, respectivamente

[T,Y]=ode45('sir',[0 t],[49999 1 0]);
% Aqui T � o tempo em dias e Y a matriz com a solu��o do modelo. O pacote
% ode45 est� considerando a fun��o "sir", que cont�m as equa��es do modelo,
% no intervalo de tempo inicial "0" e tempo final "t" (que � o valor
% definido acima). J� nos �ltimos argumentos, temos como "49999 a
% quantidade de indiv�duos Suscet�veis, 1 Infectado e 0 Recuperados. Note
% que a soma da quantidade de indiv�duos � exatamente igual ao N (a
% popula��o total considerada no modelo).

%% ================ Plotando as Solu��es ====================
figure(1) %plotando as 3 solu��es do modelo no mesmo gr�fico
plot(T,Y(:,1),'k',T,Y(:,2),'--k',T,Y(:,3),'-.k')
xlabel('Tempo (dias)'),
ylabel('Popula��o Total'),
legend('Suscet�veis', 'Infectados', 'Recuperados')

figure(2) %plotando somente a solu��o de infectados do modelo
plot(T,Y(:,2),'--k')
xlabel('Tempo (dias)'),
ylabel('Popula��o Total de Infectados'),
legend( 'Infectados')

