function B1 = BdyMassMatrix(nodes)
    segment_length = norm(nodes(:, 1) - nodes(:, 2));
    B1 = segment_length * [[1/3 1/6]; [1/6 1/3]];
end