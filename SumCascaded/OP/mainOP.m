clear all
close all

% A piece of code to calculate the Outage Probability
% Pedro Henrique Dornelas Almeida - 03/11/2023

functions_path = "functions";
addpath(functions_path);

verify_python();

% Se o caminho do Python foi encontrado, continue com o código
disp('Continuing with code...');

% ms = 2;
% alpha = 2;
% mu = 1;
% z = 1.5;

% % FoxH "arguments" pre computable assuming constant parameters
% % Ai = (1 - ms, 1/alpha)
% Ai1 = [1-ms, 1/alpha]
% Ai2 = [1-ms, 1/alpha]

% % Bi = (((z_PE**2/alpha) + 1), 1/alpha)
% Bi1 = (((z_PE^2/alpha) + 1), 1/alpha)
% Bi2 = (((z_PE^2/alpha) + 1), 1/alpha)

% % Ci = (mu, 1 / alpha)  % a tuple with a pair 
% Ci1 = (mu, 1 / alpha)   % a tuple with a pair 
% Ci2 = (mu, 1 / alpha)   % a tuple with a pair 

% % Di = (z_PE**2/alpha, 1/alpha)    % a list of tuples of pairs
% Di1 = (z_PE^2/alpha, 1/alpha)    % a list of tuples of pairs
% Di2 = (z_PE^2/alpha, 1/alpha)    % a list of tuples of pairs

% Epsilon1 = (1, 1)
% Epsilon2 = (0, 1/2)
% Epsilon3 = (1, 1/2)

% % compute analytic Outage Probability
% % a -> top right
% % b -> bot right
% % c -> top left
% % d -> bot left

% % Parameters for FoxH
% a = [[(1,1), Ai1, Ai2, Bi1, Bi2]] * 1
% b = [(Ci1, Ci2, Di1, Di2)] * 1
% c = [tuple([1] + [1/2] * 1)]                             % Epsilon3
% d = [(tuple([1] + [1] * 1)), (tuple([0] + [1/2] * 1))]   % Epsilon1 and Epsilon2

% mn = [(0, 1)] + [(4, 3)] * 1
% pq = [(1, 2)] + [(5, 4)] * 1
param = z, mn, pq, c, d, a, b

pyModule = py.importlib.import_module('multiFoxH');
H = pyModule.compMultiFoxH(param, nsubdivisions=35, boundaryTol=1e-5)