clear, clc, clf

%% Question 1.a)

beta = 5.2;                                                                 % Set values from the exercise
C = 1.5;
c_L = 6.78;
c_eta = 0.4;
p_0 = 2;

Re_lambda = [30, 300, 1500];                                                % 3 chosen values for Re_lambda
Re_L = (3/20)*Re_lambda.^2;                                                 % Calculating Re_L from Re_lambda
K = logspace(-6,1,500);                                                     % Creating a list of values for K between the minimum and maximum x-axis value
K1 = logspace(-6,1,500);                                                    % Creating a list of values for K between the minimum and maximum x-axis value
E_11 = zeros(1,500);                                                        % Creating an empty list to store E_11 values, same length as K1
figure(1)
axes('XScale', 'log', 'YScale', 'log')                                      % Forcing the loglog scale of the plot
axis([10^(-6) 10 10^(-6) 10^8])                                             % Setting axis limits
xlabel('$\tilde{\kappa}, \tilde{\kappa}_1$', 'Interpreter', 'latex')                % Axis lables
ylabel('$\tilde{E}(\tilde{\kappa}), \tilde{E}_{11}(\tilde{\kappa}_1)$', 'Interpreter','latex')

TKE = zeros(1,length(Re_lambda));                                           % Initialization for the variables used in 1.b)
fun2 = zeros(1,length(Re_lambda));
eps = zeros(1,length(Re_lambda));
Re_L2 = zeros(1,length(Re_lambda));
Re_lambda2 = zeros(1,length(Re_lambda));

for i = 1:length(Re_L)                                                      % For-loop to go through the 3 chosen Reynolds numbers
    E = @(K) C.*K.^(-5/3).*((K.*Re_L(i)^(3/4)./(((K.^2.*Re_L(i).^(3/2))+c_L).^(1/2)))...
       .^((5/3)+p_0)).*exp(-beta.*(((K.^4+c_eta.^4).^(1/4))-c_eta));        % Determining the value of E(K) as a function of K

    for j = 1:length(K1)
        fun = @(K) (E(K)./K).*(1-(K1(j).^2./K.^2));                         % Function for E_11
        E_11(j) = integral(fun,K1(j),Inf);                                  % Determining the value of E_11 by integrating over the range for a changing K1
    end

    hold on                                                                 % Plotting everything
    loglog(K,E(K))
    loglog(K1,E_11, '--')
    hold off
    
    %% Question 1.b)
    TKE(i) = integral(E,0,inf);                                             % Calculating turbulent kinetic energy
                                                                  
    fun2 = @(K) K.^2.*E(K);                                                 % Function in the integral for the TKE dissipation
    eps(i) = integral(fun2,0,inf);                                          % Calculating the TKE dissipation                                              
                                         
    Re_L2(i) = TKE(i)^2;                                                    % Calculating Re_L from the numerical integration                                                 
    Re_lambda2(i) = sqrt(20/3)*TKE(i);                                      % Calculating Re_lambda from the numerical integration    

end

h = get(gca, 'Children');                                                   % Setting colours for more clarity (could probably have been done easier though)
set([h(1), h(2)], 'Color', 'r');                                            % Colouring for Re_lambda = 1500
set([h(3), h(4)], 'Color', 'g');                                            % Colouring for Re_lambda = 300
set([h(5), h(6)], 'Color', 'b');                                            % Colouring for Re_lambda = 30

legend(['$\tilde{E}(\tilde{\kappa}): Re_{\lambda} = $', num2str(Re_lambda(1))], ...
    ['$\tilde{E}_{11}(\tilde{\kappa}_1): Re_{\lambda} = $', num2str(Re_lambda(1))], ...     % Creation of legend
    ['$\tilde{E}(\tilde{\kappa}): Re_{\lambda} = $', num2str(Re_lambda(2))], ...
    ['$\tilde{E}_{11}(\tilde{\kappa}_1): Re_{\lambda} = $', num2str(Re_lambda(2))], ...
    ['$\tilde{E}(\tilde{\kappa}): Re_{\lambda} = $', num2str(Re_lambda(3))], ...
    ['$\tilde{E}_{11}(\tilde{\kappa}_1): Re_{\lambda} = $', num2str(Re_lambda(3))], 'Interpreter', 'latex')
