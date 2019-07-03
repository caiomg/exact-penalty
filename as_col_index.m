function cols_indices = as_col_index(are_cols, cols_to_select)

    assert(isempty(find(cols_to_select & ~are_cols, 1)))
    cols_indices = cumsum(are_cols);
    cols_indices = cols_indices(cols_to_select);
end