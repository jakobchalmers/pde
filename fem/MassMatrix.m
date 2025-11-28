function B0 = MassMatrix(nodes)
    % Edge vectors of triangle using 3rd node as origin
    e1 = nodes(:, 1) - nodes(:, 3); % (p1 - p3)
    e2 = nodes(:, 2) - nodes(:, 3); % (p2 - p3)
    basis = [e1, e2];
    area = det(basis) / 2; % area of triangle
    diagonal_element = area / 6; % integral for i=j
    other_element = area / 12; % integral for i not equal to j
    B0 = [[diagonal_element other_element other_element];
        [other_element, diagonal_element, other_element];
        [other_element, other_element, diagonal_element]
        ]; % integral for all i,j combinations
end