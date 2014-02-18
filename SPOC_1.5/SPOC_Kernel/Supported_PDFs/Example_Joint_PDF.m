function PDF = Example_Joint_PDF(meshed_discretization_values)
% All joint pdfs must be written as function of a matrix of
% meshed discretization values. Rows of this matrix are realizations of
% values of parameter space.
% PDF is a column vector of joint pdf values for list of meshed discretization values

global CONSTANTS

number_of_meshed_values = size(meshed_discretization_values,1);
number_of_params = size(meshed_discretization_values,2);

PDF = zeros(number_of_meshed_values,1);


%This example will be a trivial joint pdf which just returns a
% multivariate beta distribution with parameters alpha=beta=2.
alpha = 2;  
beta = 2;

for i=1:number_of_meshed_values
    joint_values = meshed_discretization_values(i,:);
    
    temp_pdf=1;
    for k = 1:number_of_params
        w0 = CONSTANTS.ParameterSpace.W0(k);
        wf = CONSTANTS.ParameterSpace.WF(k);
        temp_pdf = temp_pdf*betapdf((joint_values(k)-w0)/(wf-w0),alpha,beta)/(wf-w0);
    end
    PDF(i,1)=temp_pdf;
end




