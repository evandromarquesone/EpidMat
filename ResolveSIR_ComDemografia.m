%% Material Complementar do Livro: 
% Introdução à Epidemiologia Matemática: Métodos em Estudos Transversais

% = Outros Materiais estão disponíveis em https://linktr.ee/livroepidmat =

% ======== Programa para Simular o modelo epidemiológico SIR ============

%comandos para não levar nenhum resíduo computacional da simulação anterior
clc;  clear all;  close all; 
%% ======= Calculando Ro e o Ponto de Equilíbrio Endêmico ===============
% Esses valores de parâmetros são utilizados aqui para o cálculo de Ro e/ou 
% das coordenadas do ponto de equilíbrio endêmico. Para o modelo, alterar
% os valores no arquivo sir_ComDemografia.m. NÃO SE ESQUEÇA DISSO!
%%============= Comente caso não queria calcular ========================
n=5000; %50 mil indivíduos
beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
mu=1/(80*365); %(expectativa de 80 anos de vida)
%d=0;
nu=1/60; %(a taxa é 1 dividido pelo tempo em que se fica infectado)

%alpha=0.001; %( mortalidade pela doença de 1 a cada 1000 dias)

Ro=(beta*n)/(nu+mu)
Seq=n/Ro     
Ieq=(Ro-1)*(mu/beta)   
Ieq2=((mu*n)/(nu+mu))-(mu/beta)   
Req=(Ro-1)*nu/beta

%% ===== Calculando Dois Autovalores da Matriz Jacobiana do Sistema ====
%%============= Comente caso não queria calcular ========================
Discriminante=(mu*mu)*(Ro*Ro)-(4*mu*(mu+nu)*(Ro-1))
Autovalor1=0.5*((-mu*Ro)+sqrt(Discriminante))
Autovalor2=0.5*((-mu*Ro)-sqrt(Discriminante))
%% =======================================================================
% =================  Resolvendo o Sistema de EDO ========================
% =======================================================================
t=600; %tempo em dias

%aplicando o pacote ode45 no arquivo função sir_ComDemografia, que contém 
%as equações. O vetor [49999 1 0] são as condições iniciais de suscetíveis,  
% infectados e recuperados, respectivamente.

[T,Y]=ode45('sir_ComDemografia',[0 t],[49999 1 0]);

%% ================ Plotando as Soluções ====================
figure(1)
plot(T,Y(:,1),'k',T,Y(:,2),'--k',T,Y(:,3),'-.k')
xlabel('Tempo (dias)'),
ylabel('População Total'),
legend('Suscetíveis', 'Infectados', 'Recuperados')

figure(2)
plot(T,Y(:,2),'k')
xlabel('Tempo (dias)'),
ylabel('População Total de Infectados'),
legend( 'Infectados')
