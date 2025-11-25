function B0 = MassMatrix(nodes)
    e1 = nodes(:, 1) - nodes(:, 3); % (p1 - p3) choose 3rd node as origin
    e2 = nodes(:, 2) - nodes(:, 3); % (p2 - p3)
    basis = [e1, e2];
    J = det(basis);
    diagonal_element = J / 12;
    other_element = J / 24;
    B0 = [[diagonal_element other_element other_element];
        [other_element, diagonal_element, other_element];
        [other_element, other_element, diagonal_element]
        ];
end