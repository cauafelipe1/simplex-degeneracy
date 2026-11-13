function [A, b, c, raw_mps] = mps_parser(mps_filepath, output_mat_path, to_standard_form)
% mps_parser
%  parses a Mathematical Programming System (MPS) file into sparse Octave matrices.
%   [A, b, c, raw_mps] = mps_parser(mps_filepath, output_mat_path, to_standard_form)
%   reads the MPS sections (NAME, ROWS, COLUMNS, RHS, RANGES, BOUNDS, ENDATA),
%   maps alphanumeric row and column identifiers to numeric matrix indices in O(1),
%   assembles the technological matrix A using sparse triplet allocation, and saves
%   the resulting matrices (A, b, c, raw_mps) in binary .mat format.
%
%   hint: incremental insertion into a sparse matrix A(i, j) = v inside a loop
%   triggers O(nnz^2) Compressed Sparse Column (CSC) reallocations. To ensure high
%   benchmark performance, all non-zero coefficients are accumulated in preallocated
%   triplet vectors (row_idx, col_idx, val_vec) and assembled in a single call to
%   sparse(row_idx, col_idx, val_vec, m, n) in O(nnz * log(nnz)) time.
%   When to_standard_form is true (default), the raw MPS problem:
%    min c_raw^T x   s.t.   L/G/E constraints on A_raw * x and b_raw,  lb <= x <= ub
%   is algebraically transformed into Vanderbei's canonical standard form:
%    max c^T x       s.t.   A * x <= b,  x >= 0
%   via:
%    1. Variable bound shifting: x_j' = x_j - lb_j for finite lb_j != 0
%    2. Free variable splitting: x_j = x_j^+ - x_j^- for lb_j = -Inf
%    3. Inequality normalization:
%       - 'L' rows (a_i^T x <= b_i) -> kept as a_i^T x <= b_i
%       - 'G' rows (a_i^T x >= b_i) -> negated to -a_i^T x <= -b_i
%       - 'E' rows (a_i^T x == b_i) -> split into a_i^T x <= b_i and -a_i^T x <= -b_i
%    4. Finite upper bounds: appended as sparse rows x_j' <= ub_j - lb_j
%    5. Objective negation: c = -c_raw (converting minimization to maximization).

    % 1. Defensive Programming and Input Validation
    if nargin < 1 || isempty(mps_filepath)
        error("mps_parser: mps_filepath argument is required and cannot be empty");
    endif

    if ~ischar(mps_filepath)
        error("mps_parser: mps_filepath must be a character string");
    endif

    if exist(mps_filepath, "file") != 2
        error(sprintf("mps_parser: file not found at '%s'", mps_filepath));
    endif

    if nargin < 2 || isempty(output_mat_path)
        [~, base_name, ~] = fileparts(mps_filepath);
        output_mat_path = fullfile("data", "processed", [base_name, ".mat"]);
    endif

    if ~ischar(output_mat_path)
        error("mps_parser: output_mat_path must be a character string");
    endif

    if nargin < 3 || isempty(to_standard_form)
        to_standard_form = true;
    endif

    fid = fopen(mps_filepath, "r");
    if fid == -1
        error(sprintf("mps_parser: unable to open file '%s' for reading", mps_filepath));
    endif

    % 2. Data Structures and Preallocated Triplet Buffers
    model_name = "";
    section = "";
    obj_row_name = "";
    active_rhs_name = "";
    active_range_name = "";
    active_bound_name = "";

    row_map = containers.Map("KeyType", "char", "ValueType", "double");
    col_map = containers.Map("KeyType", "char", "ValueType", "double");

    row_names = {};
    col_names = {};
    row_types = "";

    m_raw = 0;
    n_raw = 0;

    triplet_cap = 65536;
    nnz_count = 0;
    row_idx_buf = zeros(triplet_cap, 1);
    col_idx_buf = zeros(triplet_cap, 1);
    val_buf = zeros(triplet_cap, 1);

    c_cap = 4096;
    c_raw = zeros(c_cap, 1);

    rhs_entries_row = [];
    rhs_entries_val = [];

    range_entries_row = [];
    range_entries_val = [];

    bound_types = {};
    bound_cols = [];
    bound_vals = [];

    % 3. Sequential Line-by-Line State Machine
    while ~feof(fid)
        raw_line = fgetl(fid);
        if ~ischar(raw_line)
            break;
        endif

        % Strip carriage return (\r) if present, preserving leading spaces
        clean_line = regexprep(raw_line, "[\r\n]+$", "");
        if isempty(strtrim(clean_line))
            continue;
        endif

        % Skip full-line comments starting with '*' in column 1 or '$'
        first_char = clean_line(1);
        if first_char == "*" || first_char == "$"
            continue;
        endif

        % Section Header Detection (non-whitespace in column 1)
        if ~isspace(first_char)
            header_tokens = strsplit(strtrim(clean_line));
            sec_keyword = upper(header_tokens{1});

            if strcmp(sec_keyword, "NAME")
                if numel(header_tokens) >= 2
                    model_name = header_tokens{2};
                else
                    model_name = "UNNAMED";
                endif
                section = "NAME";
            elseif strcmp(sec_keyword, "ROWS")
                section = "ROWS";
            elseif strcmp(sec_keyword, "COLUMNS")
                section = "COLUMNS";
            elseif strcmp(sec_keyword, "RHS")
                section = "RHS";
            elseif strcmp(sec_keyword, "RANGES")
                section = "RANGES";
            elseif strcmp(sec_keyword, "BOUNDS")
                section = "BOUNDS";
            elseif strcmp(sec_keyword, "ENDATA")
                section = "ENDATA";
                break;
            else
                section = sec_keyword;
            endif
            continue;
        endif

        % Tokenize data record
        tokens = strsplit(strtrim(clean_line));
        num_tokens = numel(tokens);
        if num_tokens == 0
            continue;
        endif

        % Dispatch according to active MPS section
        if strcmp(section, "ROWS")
            if num_tokens < 2
                continue;
            endif
            r_type = upper(tokens{1}(1));
            r_name = tokens{2};

            if r_type == "N"
                if isempty(obj_row_name)
                    obj_row_name = r_name;
                endif
            elseif r_type == "L" || r_type == "G" || r_type == "E"
                if ~isKey(row_map, r_name)
                    m_raw = m_raw + 1;
                    row_map(r_name) = m_raw;
                    row_names{m_raw, 1} = r_name;
                    row_types(m_raw, 1) = r_type;
                endif
            else
                error(sprintf("mps_parser: unsupported row type '%s' for row '%s'", r_type, r_name));
            endif

        elseif strcmp(section, "COLUMNS")
            % Ignore integer 'MARKER' lines ('INTORG' / 'INTEND')
            if num_tokens >= 2 && strcmp(tokens{2}, "'MARKER'")
                continue;
            endif

            c_name = tokens{1};
            if isKey(col_map, c_name)
                j_col = col_map(c_name);
            else
                n_raw = n_raw + 1;
                j_col = n_raw;
                col_map(c_name) = j_col;
                col_names{j_col, 1} = c_name;
                if j_col > c_cap
                    c_raw = [c_raw; zeros(c_cap, 1)];
                    c_cap = c_cap * 2;
                endif
            endif

            % Process up to two (row_name, value) pairs on the line
            tok_pos = 2;
            while tok_pos + 1 <= num_tokens
                r_name = tokens{tok_pos};
                if r_name(1) == "$"
                    break;
                endif
                val_num = str2double(tokens{tok_pos + 1});
                if isnan(val_num)
                    error(sprintf("mps_parser: invalid numeric value '%s' in COLUMNS for (%s, %s)", ...
                                  tokens{tok_pos + 1}, c_name, r_name));
                endif

                if strcmp(r_name, obj_row_name)
                    c_raw(j_col) = c_raw(j_col) + val_num;
                elseif isKey(row_map, r_name)
                    i_row = row_map(r_name);
                    nnz_count = nnz_count + 1;
                    if nnz_count > triplet_cap
                        row_idx_buf = [row_idx_buf; zeros(triplet_cap, 1)];
                        col_idx_buf = [col_idx_buf; zeros(triplet_cap, 1)];
                        val_buf = [val_buf; zeros(triplet_cap, 1)];
                        triplet_cap = triplet_cap * 2;
                    endif
                    row_idx_buf(nnz_count) = i_row;
                    col_idx_buf(nnz_count) = j_col;
                    val_buf(nnz_count) = val_num;
                endif

                tok_pos = tok_pos + 2;
            endwhile

        elseif strcmp(section, "RHS")
            % Handle both named RHS records (3 or 5 tokens) and unnamed RHS records (2 or 4 tokens)
            if mod(num_tokens, 2) == 1
                rhs_vec_name = tokens{1};
                start_pos = 2;
            else
                rhs_vec_name = "__DEFAULT_RHS__";
                start_pos = 1;
            endif

            if isempty(active_rhs_name)
                active_rhs_name = rhs_vec_name;
            endif

            % Only parse the first RHS vector encountered in the file
            if strcmp(rhs_vec_name, active_rhs_name)
                tok_pos = start_pos;
                while tok_pos + 1 <= num_tokens
                    r_name = tokens{tok_pos};
                    if r_name(1) == "$"
                        break;
                    endif
                    val_num = str2double(tokens{tok_pos + 1});
                    if isKey(row_map, r_name)
                        i_row = row_map(r_name);
                        rhs_entries_row(end + 1, 1) = i_row;
                        rhs_entries_val(end + 1, 1) = val_num;
                    endif
                    tok_pos = tok_pos + 2;
                endwhile
            endif

        elseif strcmp(section, "RANGES")
            if mod(num_tokens, 2) == 1
                rng_vec_name = tokens{1};
                start_pos = 2;
            else
                rng_vec_name = "__DEFAULT_RANGE__";
                start_pos = 1;
            endif

            if isempty(active_range_name)
                active_range_name = rng_vec_name;
            endif

            if strcmp(rng_vec_name, active_range_name)
                tok_pos = start_pos;
                while tok_pos + 1 <= num_tokens
                    r_name = tokens{tok_pos};
                    if r_name(1) == "$"
                        break;
                    endif
                    val_num = str2double(tokens{tok_pos + 1});
                    if isKey(row_map, r_name)
                        i_row = row_map(r_name);
                        range_entries_row(end + 1, 1) = i_row;
                        range_entries_val(end + 1, 1) = val_num;
                    endif
                    tok_pos = tok_pos + 2;
                endwhile
            endif

        elseif strcmp(section, "BOUNDS")
            if num_tokens < 2
                continue;
            endif
            b_type = upper(tokens{1});

            % Bound types without numeric value: FR (free), MI (-Inf), PL (+Inf), BV (binary 0-1)
            if strcmp(b_type, "FR") || strcmp(b_type, "MI") || strcmp(b_type, "PL") || strcmp(b_type, "BV")
                if num_tokens >= 3
                    bnd_name = tokens{2};
                    c_name = tokens{3};
                else
                    bnd_name = "__DEFAULT_BND__";
                    c_name = tokens{2};
                endif
                val_num = 0.0;
            else
                % Bound types with numeric value: LO, UP, FX, LI, UI
                if num_tokens >= 4
                    bnd_name = tokens{2};
                    c_name = tokens{3};
                    val_num = str2double(tokens{4});
                elseif num_tokens == 3
                    bnd_name = "__DEFAULT_BND__";
                    c_name = tokens{2};
                    val_num = str2double(tokens{3});
                else
                    continue;
                endif
            endif

            if isempty(active_bound_name)
                active_bound_name = bnd_name;
            endif

            if strcmp(bnd_name, active_bound_name) && isKey(col_map, c_name)
                j_col = col_map(c_name);
                bound_types{end + 1, 1} = b_type;
                bound_cols(end + 1, 1) = j_col;
                bound_vals(end + 1, 1) = val_num;
            endif
        endif
    endwhile

    fclose(fid);

    % 4. Validate Extracted Dimensions and Build Sparse Raw Matrix
    if m_raw == 0 || n_raw == 0
        error(sprintf("mps_parser: invalid or empty MPS model in '%s' (m=%d, n=%d)", ...
                      mps_filepath, m_raw, n_raw));
    endif

    c_raw = c_raw(1:n_raw);
    b_raw = zeros(m_raw, 1);
    for k = 1:numel(rhs_entries_row)
        b_raw(rhs_entries_row(k)) = rhs_entries_val(k);
    endfor

    % Single sparse assembly call in O(nnz * log(nnz))
    A_raw = sparse(row_idx_buf(1:nnz_count), col_idx_buf(1:nnz_count), ...
                   val_buf(1:nnz_count), m_raw, n_raw);

    % Construct variable bounds vectors (default MPS bounds: 0 <= x < +Inf)
    lb = zeros(n_raw, 1);
    ub = inf(n_raw, 1);

    for k = 1:numel(bound_cols)
        j_col = bound_cols(k);
        b_type = bound_types{k};
        b_val = bound_vals(k);

        if strcmp(b_type, "LO") || strcmp(b_type, "LI")
            lb(j_col) = b_val;
        elseif strcmp(b_type, "UP") || strcmp(b_type, "UI")
            ub(j_col) = b_val;
        elseif strcmp(b_type, "FX")
            lb(j_col) = b_val;
            ub(j_col) = b_val;
        elseif strcmp(b_type, "FR")
            lb(j_col) = -inf;
            ub(j_col) = inf;
        elseif strcmp(b_type, "MI")
            lb(j_col) = -inf;
            % Note: in MPS convention, MI leaves upper bound at 0 unless UP/PL was specified
            if isinf(ub(j_col))
                ub(j_col) = 0.0;
            endif
        elseif strcmp(b_type, "PL")
            ub(j_col) = inf;
        elseif strcmp(b_type, "BV")
            lb(j_col) = 0.0;
            ub(j_col) = 1.0;
        endif
    endfor

    % Populate raw_mps metadata struct
    raw_mps = struct();
    raw_mps.name = model_name;
    raw_mps.obj_row_name = obj_row_name;
    raw_mps.A = A_raw;
    raw_mps.b = b_raw;
    raw_mps.c = c_raw;
    raw_mps.row_types = row_types;
    raw_mps.row_names = row_names;
    raw_mps.col_names = col_names;
    raw_mps.lb = lb;
    raw_mps.ub = ub;

    % 5. Canonical Conversion to Vanderbei Standard Form (max c^T x s.t. A x <= b, x >= 0)
    if to_standard_form
        A_work = A_raw;
        b_work = b_raw;
        c_work = c_raw;
        lb_work = lb;
        ub_work = ub;
        rtypes_work = row_types;

        % 5.1 Apply RANGES if any were present in the MPS file
        if ~isempty(range_entries_row)
            extra_A_rows = {};
            extra_b_vals = [];
            for k = 1:numel(range_entries_row)
                i_row = range_entries_row(k);
                r_val = range_entries_val(k);
                r_type = rtypes_work(i_row);
                if r_type == "G"
                    % b_i <= a_i^T x <= b_i + |r_i|
                    extra_A_rows{end + 1, 1} = A_work(i_row, :);
                    extra_b_vals(end + 1, 1) = b_work(i_row) + abs(r_val);
                elseif r_type == "L"
                    % b_i - |r_i| <= a_i^T x <= b_i
                    extra_A_rows{end + 1, 1} = -A_work(i_row, :);
                    extra_b_vals(end + 1, 1) = -(b_work(i_row) - abs(r_val));
                elseif r_type == "E"
                    if r_val > 0
                        rtypes_work(i_row) = "G";
                        extra_A_rows{end + 1, 1} = A_work(i_row, :);
                        extra_b_vals(end + 1, 1) = b_work(i_row) + r_val;
                    elseif r_val < 0
                        rtypes_work(i_row) = "L";
                        extra_A_rows{end + 1, 1} = -A_work(i_row, :);
                        extra_b_vals(end + 1, 1) = -(b_work(i_row) + r_val);
                    endif
                endif
            endfor
            if ~isempty(extra_A_rows)
                A_extra_rng = vertcat(extra_A_rows{:});
                A_work = [A_work; A_extra_rng];
                b_work = [b_work; extra_b_vals];
                rtypes_work = [rtypes_work; repmat("L", numel(extra_b_vals), 1)];
            endif
        endif

        % 5.2 Handle variables with lb = -Inf and finite ub (x_j <= ub_j) via negation x_j' = ub_j - x_j >= 0
        neg_inf_bounded = find(isinf(lb_work) & (lb_work < 0) & ~isinf(ub_work));
        if ~isempty(neg_inf_bounded)
            b_work = b_work - A_work(:, neg_inf_bounded) * ub_work(neg_inf_bounded);
            A_work(:, neg_inf_bounded) = -A_work(:, neg_inf_bounded);
            c_work(neg_inf_bounded) = -c_work(neg_inf_bounded);
            lb_work(neg_inf_bounded) = 0.0;
            ub_work(neg_inf_bounded) = inf;
        endif

        % 5.3 Shift finite non-zero lower bounds: x_j' = x_j - lb_j >= 0
        finite_lb_mask = ~isinf(lb_work) & (lb_work != 0);
        if any(finite_lb_mask)
            shift_vec = zeros(size(lb_work));
            shift_vec(finite_lb_mask) = lb_work(finite_lb_mask);
            b_work = b_work - A_work * shift_vec;
            raw_mps.obj_const_shift = c_work' * shift_vec;
            ub_work(finite_lb_mask) = ub_work(finite_lb_mask) - lb_work(finite_lb_mask);
            lb_work(finite_lb_mask) = 0.0;
        else
            raw_mps.obj_const_shift = 0.0;
        endif

        % 5.4 Split remaining free variables (lb = -Inf, ub = +Inf): x_j = x_j^+ - x_j^-
        free_cols = find(isinf(lb_work) & (lb_work < 0) & isinf(ub_work));
        if ~isempty(free_cols)
            A_work = [A_work, -A_work(:, free_cols)];
            c_work = [c_work; -c_work(free_cols)];
            ub_work = [ub_work; inf(numel(free_cols), 1)];
            lb_work(free_cols) = 0.0;
            lb_work = [lb_work; zeros(numel(free_cols), 1)];
        endif

        % 5.5 Convert all row constraints ('L', 'G', 'E') into canonical 'A * x <= b'
        idx_L = find(rtypes_work == "L");
        idx_G = find(rtypes_work == "G");
        idx_E = find(rtypes_work == "E");

        A_ineq = [A_work(idx_L, :); ...
                  -A_work(idx_G, :); ...
                  A_work(idx_E, :); ...
                  -A_work(idx_E, :)];

        b_ineq = [b_work(idx_L); ...
                  -b_work(idx_G); ...
                  b_work(idx_E); ...
                  -b_work(idx_E)];

        % 5.6 Append finite upper bounds as sparse rows: x_j' <= ub_j'
        finite_ub_cols = find(~isinf(ub_work));
        if ~isempty(finite_ub_cols)
            num_ub = numel(finite_ub_cols);
            n_curr = size(A_ineq, 2);
            A_ub = sparse(1:num_ub, finite_ub_cols, ones(num_ub, 1), num_ub, n_curr);
            A_ineq = [A_ineq; A_ub];
            b_ineq = [b_ineq; ub_work(finite_ub_cols)];
        endif

        % 5.7 Convert objective from MPS minimization to Vanderbei maximization (max c^T x)
        A = A_ineq;
        b = b_ineq;
        c = -c_work;
    else
        A = A_raw;
        b = b_raw;
        c = c_raw;
    endif

    % 6. Final Defensive Dimension and Sparsity Checks
    if ~issparse(A)
        A = sparse(A);
    endif

    [m_final, n_final] = size(A);
    if length(b) != m_final
        error(sprintf("mps_parser: dimension mismatch between A (%d rows) and b (%d elements)", ...
                      m_final, length(b)));
    endif

    if length(c) != n_final
        error(sprintf("mps_parser: dimension mismatch between A (%d cols) and c (%d elements)", ...
                      n_final, length(c)));
    endif

    % 7. Save Binary .mat File
    [out_dir, ~, ~] = fileparts(output_mat_path);
    if ~isempty(out_dir) && exist(out_dir, "dir") != 7
        mkdir(out_dir);
    endif

    save("-mat7-binary", output_mat_path, "A", "b", "c", "raw_mps");
endfunction
