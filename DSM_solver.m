% Extensional Stiffness
EA1 = 21000000;
EA2 = 10000000;
EA = [EA1, EA1, EA1, EA1, EA2];

% Truss Lengths
L = [60, 60, sqrt(45^2 + 15^2), sqrt(45^2 + 15^2), sqrt(2) * 45];

% Angles
alpha = [-pi/2, 0, -atan(45/15), -atan(15/45), pi/4];

% Problem Data
n_ele = 5; % Number of elements

% ID Matrix of the System
ID_t = [
    1, 2, 3, 4;
    3, 4, 7, 8;
    5, 6, 7, 8;
    1, 2, 5, 6;
    3, 4, 5, 6
];

ID = transpose(ID_t);      %Real ID matrix of the system  

% Element Stiffness Matrices
Q = cell(n_ele, 1); % Transformation matrices
k = cell(n_ele, 1); % Element stiffness matrices

for i = 1:n_ele
    % Local stiffness matrix
    k_bar = (EA(i) / L(i)) * [1, 0, -1, 0; 0, 0, 0, 0; -1, 0, 1, 0; 0, 0, 0, 0];
    
    % Transformation matrix
    c = cos(alpha(i));
    s = sin(alpha(i));
    Q{i} = [c, s, 0, 0; -s, c, 0, 0; 0, 0, c, s; 0, 0, -s, c];
    
    % Global element stiffness matrix
    k{i} = Q{i}' * k_bar * Q{i};
end

% Global Stiffness Matrix Assembly
K = zeros(8); % Initialize global stiffness matrix

for e = 1:n_ele
    for i = 1:4
        for j = 1:4
            K(ID(i, e), ID(j, e)) = K(ID(i, e), ID(j, e)) + k{e}(i, j);
        end
    end
end


% Applying Boundary Conditions
% No deformation in directions 1, 2, 3 (D1 = D2 = D3 = 0)
K_red = K(4:8, 4:8);

% Reduced prescribed forces
F_red = [0; 0; 0; 0; 200];

% Solving for reduced global displacements
D_red = K_red \ F_red;

% Constructing full displacement vector
D = [0; 0; 0; D_red];

% Compute system forces (prescribed and reactions)
F = K * D;

% Recover member forces
d = cell(n_ele, 1); % Local displacement vectors
d_bar = cell(n_ele, 1); % Member local displacements
N = zeros(n_ele, 1); % Member forces

for i = 1:n_ele
    d{i} = [D(ID(1, i)); D(ID(2, i)); D(ID(3, i)); D(ID(4, i))];
    d_bar{i} = Q{i} * d{i};
    N(i) = EA(i) * (d_bar{i}(3) - d_bar{i}(1)) / L(i);
end

% Display Results
disp('Nodal Displacements:');
disp(D);

disp('Support Reactions:');
disp(F);

disp('Member Forces:');
disp(N);
