clc

% Represent the given data
figure;
mesh(id.X{1}, id.X{2}, id.Y);
colormap("cool");
title('Initial data');
xlabel('Input data 1');
ylabel('Input data 2');
zlabel('Output value');

%%
clc

% Declare given var for ease of use
X1_id = id.X{1};
X2_id = id.X{2};
Y_id = id.Y;

X1_val = val.X{1};
X2_val = val.X{2};
Y_val = val.Y;

% Initialize the MSE vectors
MSE_id_all = zeros(26,1);
MSE_val_all = zeros(26,1);

m = 26; 
offset = -150;

% Outher for = degree, Inner fors = indices
% We use an outher for because we need all the degrees' MSE for the final plot
for d = 1:m

    terms = ((d + 1) * (d + 2)) / 2;

    phi_id = zeros(length(X1_id)^2, terms);

    phi_id(:,1) = ones(length(X1_id)^2, 1);

    c = 2;

    for i = 1 : d
        for j = 0 : i
            x1_temp = X1_id.^(i - j);
            x2_temp = X2_id.^j;

            col = zeros(length(X1_id) * length(X2_id), 1);

            index = 1;
            for i1 = 1:length(x1_temp)
                for j2 = 1:length(x2_temp)
                    col(index) = x1_temp(i1) * x2_temp(j2);
                    index = index + 1;
                end
            end

            phi_id(:,c) = col;
            c = c + 1;
        end
    end

    Y_id_flat = reshape(Y_id, [], 1);

    theta = phi_id \ Y_id_flat;

    Y_id_approx = phi_id * theta;

    Y_id_approx_matrix = reshape(Y_id_approx, 41, 41);

    % Plot only for the optimal degree
    if(d == 6)
        figure;
        sgtitle('Id data vs Approximator');
        subplot(121);
        mesh(X1_id, X2_id, Y_id);
        colormap("cool");
        hold on;
        mesh(X1_id, X2_id, Y_id_approx_matrix);
        xlabel('Input data 1');
        ylabel('Input data 2');
        zlabel('Output value');

        subplot(122);
        mesh(X1_id, X2_id, Y_id);
        colormap("cool");
        hold on;
        subplot(122);
        mesh(X1_id, X2_id, Y_id_approx_matrix + offset);
        xlabel('Input data 1');
        ylabel('Input data 2');
        zlabel('Output value');
    end


    phi_val = zeros(length(X1_val)^2, terms);

    phi_val(:,1) = ones(length(X1_val)^2, 1);

    c = 2;

    for i = 1 : d
        for j = 0 : i
            x1_temp = X1_val.^(i - j);
            x2_temp = X2_val.^j;

            col = zeros(length(X1_val) * length(X2_val), 1);

            index = 1;
            for i1 = 1:length(x1_temp)
                for j2 = 1:length(x2_temp)
                    col(index) = x1_temp(i1) * x2_temp(j2);
                    index = index + 1;
                end
            end

            phi_val(:,c) = col;
            c = c + 1;
        end
    end

    Y_val_flat = reshape(Y_val, [], 1);

    Y_val_approx = phi_val * theta;

    Y_val_approx_matrix = reshape(Y_val_approx, 31, 31);

    % Creating the MSEs
    MSE_id = sum(1/length(Y_id) * (Y_id_flat - Y_id_approx).^2);
    MSE_val = sum(1/length(Y_val) * (Y_val_flat - Y_val_approx).^2);

    MSE_id_all(d) = MSE_id;
    MSE_val_all(d) = MSE_val;

    % Plot only for the optimal degree
    if(d == 6)
        figure;
        sgtitle('Val data vs Approximator');
        subplot(121);
        mesh(X1_val, X2_val, Y_val);
        colormap("cool");
        hold on;
        mesh(X1_val, X2_val, Y_val_approx_matrix);
        xlabel('Input data 1');
        ylabel('Input data 2');
        zlabel('Output value');

        subplot(122)
        mesh(X1_val, X2_val, Y_val);
        colormap("cool");
        hold on;
        mesh(X1_val, X2_val, Y_val_approx_matrix + offset);
        xlabel('Input data 1');
        ylabel('Input data 2');
        zlabel('Output value');
    end
end

% Plot of the MSEs
figure;
n = 1 : m;
plot(n, MSE_val_all, '-c', 'DisplayName', 'MSE VAL');
hold on;
plot(n, MSE_id_all, '-m', 'DisplayName', 'MSE ID');
grid;
hold on;
plot(6, 387.684, '*k', 'DisplayName', 'Optimal degree');
text(6, 387.684, 'm=6', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
title('MSE comparison');
xlabel('Degree m');
ylabel('MSE');
legend show;