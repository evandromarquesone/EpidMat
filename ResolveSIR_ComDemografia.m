% ======== Programa para Simular o modelo epidemiol�gico SIR ============

%comandos para n�o levar nenhum res�duo computacional da simula��o anterior
clc;  clear all;  close all; 
%% ======= Calculando Ro e o Ponto de Equil�brio End�mico ===============
% Esses valores de par�metros s�o utilizados aqui para o c�lculo de Ro e/ou 
% das coordenadas do ponto de equil�brio end�mico. Para o modelo, alterar
% os valores no arquivo sir_ComDemografia.m. N�O SE ESQUE�A DISSO!
%%============= Comente caso n�o queria calcular ========================
n=5000; %50 mil indiv�duos
beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
mu=1/(80*365); %(expectativa de 80 anos de vida)
%d=0;
nu=1/60; %(a taxa � 1 dividido pelo tempo em que se fica infectado)

%alpha=0.001; %( mortalidade pela doen�a de 1 a cada 1000 dias)

Ro=(beta*n)/(nu+mu)
Seq=n/Ro     
Ieq=(Ro-1)*(mu/beta)   
Ieq2=((mu*n)/(nu+mu))-(mu/beta)   
Req=(Ro-1)*nu/beta

%% ===== Calculando Dois Autovalores da Matriz Jacobiana do Sistema ====
%%============= Comente caso n�o queria calcular ========================
Discriminante=(mu*mu)*(Ro*Ro)-(4*mu*(mu+nu)*(Ro-1))
Autovalor1=0.5*((-mu*Ro)+sqrt(Discriminante))
Autovalor2=0.5*((-mu*Ro)-sqrt(Discriminante))
%% =======================================================================
% =================  Resolvendo o Sistema de EDO ========================
% =======================================================================
t=600; %tempo em dias

%aplicando o pacote ode45 no arquivo fun��o sir_ComDemografia, que cont�m 
%as equa��es. O vetor [49999 1 0] s�o as condi��es iniciais de suscet�veis,  
% infectados e recuperados, respectivamente.

[T,Y]=ode45('sir_ComDemografia',[0 t],[49999 1 0]);

%% ================ Plotando as Solu��es ====================
figure(1)
plot(T,Y(:,1),'k',T,Y(:,2),'--k',T,Y(:,3),'-.k')
xlabel('Tempo (dias)'),
ylabel('Popula��o Total'),
legend('Suscet�veis', 'Infectados', 'Recuperados')

figure(2)
plot(T,Y(:,2),'k')
xlabel('Tempo (dias)'),
ylabel('Popula��o Total de Infectados'),
legend( 'Infectados')