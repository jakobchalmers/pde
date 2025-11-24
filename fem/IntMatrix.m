function A0 = IntMatrix(nodes)
    e1 = nodes(:, 1) - nodes(:, 3); % (p1 - p3) choose 3rd node as origin
    e2 = nodes(:, 2) - nodes(:, 3); % (p2 - p3)
    basis = [e1, e2];
    dualbasis = inv(basis');
    grads = [dualbasis(:, 1), dualbasis(:, 2), - dualbasis(:, 1) - dualbasis(:,2)]; % [grad1, grad2, grad3]
    area = det(basis) / 2; % half a parallellogram
    A0 = grads' * grads * area;
end
