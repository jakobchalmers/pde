function A0 = BdyStiffnessMatrix(nodes)
    % Edge vectors of triangle using 3rd node as origin
    e1 = nodes(:, 1) - nodes(:, 3); % (p1 - p3)
    e2 = nodes(:, 2) - nodes(:, 3); % (p2 - p3)

    % Solve dot(e_i, grad phi_j) = 1, i=j and dot(e_i, grad phi_j) = 0, i!=j
    basis = [e1, e2];
    e1
    e2
    basis
    dualbasis = inv(basis');

    % [grad phi_1, grad phi_2, grad phi_3 = - (grad phi_2 + grad phi_3]
    grads = [dualbasis(:, 1), dualbasis(:, 2), - dualbasis(:, 1) - dualbasis(:,2)]; 

    area = det(basis) / 2; % triangle area
    % integral of dot(grad phi_i, grad phi_j) for all combinations i, j
    A0 = grads' * grads * area;
end
