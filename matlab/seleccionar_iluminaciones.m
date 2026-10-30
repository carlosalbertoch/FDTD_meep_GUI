% FDTD_meep v1.0.0 | 2026-10-02 | Carlos Chacón
function seleccionar_iluminaciones(config_json)

if nargin < 1
    config_json = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'inputs', 'config', 'phantom_4pi.json');
end
raiz = fileparts(fileparts(mfilename('fullpath')));
cfg = jsondecode(fileread(config_json));
n_sel = 32; tukey = 0.10;
if isfield(cfg, 'selection')
    if isfield(cfg.selection, 'n'), n_sel = cfg.selection.n; end
    if isfield(cfg.selection, 'tukey'), tukey = cfg.selection.tukey; end
end
salida = ruta(raiz, cfg.output);
vistas = fieldnames(cfg.views);

for v = 1:numel(vistas)
    vista = vistas{v};
    carpeta = ruta(raiz, cfg.views.(vista));
    out = fullfile(salida, vista);
    sel_dir = fullfile(out, 'selection');
    if ~exist(sel_dir, 'dir'), mkdir(sel_dir); end

    c = leer_config(fullfile(carpeta, 'config_manip.txt'));
    Gtot = (c.F_TUBE / c.F_OBJ) / c.RF;
    FOV = c.DIM_ROI * c.TPCAM / Gtot;
    REwald = (2*pi*c.N0/c.LAMBDA) / (2*pi/FOV);
    dimHolo = 2 * ceil(REwald * c.NA / c.N0);

    files = dir(fullfile(carpeta, 'i*.pgm'));
    files = files(~cellfun(@isempty, regexp({files.name}, '^i\d+\.pgm$', 'once')));
    idx = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), {files.name});
    [idx, o] = sort(idx); files = files(o);
    N = numel(files);

    img0 = double(imread(fullfile(carpeta, files(1).name)));
    [Ny, Nx] = size(img0);
    [csx, csy] = coord_shift(c.CIRCLE_CX, c.CIRCLE_CY, Nx, Ny);
    rows = mod((csy - dimHolo/2 + 1 : csy + dimHolo/2) - 1, Ny) + 1;
    cols = mod((csx - dimHolo/2 + 1 : csx + dimHolo/2) - 1, Nx) + 1;
    pc = dimHolo/2 + 1;
    w1 = tukeywin_local(Nx, tukey); W = sqrt(w1(:) * w1(:)');

    T = table();
    T.index = idx(:); T.filename = string({files.name})';
    [T.peak_row, T.peak_col, T.peak_row_subpix, T.peak_col_subpix, ...
     T.drow_subpix, T.dcol_subpix, T.peak_quality, T.peak_snr] = deal(zeros(N, 1));
    for k = 1:N
        F = fft2(W .* double(imread(fullfile(carpeta, files(k).name))));
        P = abs(F(rows, cols)).^2;
        [pk, id] = max(P(:)); [pr, pcc] = ind2sub(size(P), id);
        m = centroide(P, pr, pcc, 3);
        T.peak_row(k) = pr; T.peak_col(k) = pcc;
        T.peak_row_subpix(k) = m.r; T.peak_col_subpix(k) = m.c;
        T.drow_subpix(k) = m.r - pc; T.dcol_subpix(k) = m.c - pc;
        T.peak_quality(k) = pk / (mean(P(:)) + eps);
        T.peak_snr(k) = pk / (m.bg + eps);
    end
    T.kx_norm_raw = T.dcol_subpix / REwald;
    T.ky_norm_raw = T.drow_subpix / REwald;
    T.rho_norm = hypot(T.kx_norm_raw, T.ky_norm_raw);
    T.theta_deg = asind(min(T.rho_norm, 1));
    T.phi_deg = atan2d(T.ky_norm_raw, T.kx_norm_raw);

    K = [T.kx_norm_raw, T.ky_norm_raw];
    orden = zeros(N, 1); [~, orden(1)] = min(T.rho_norm);
    dmin = hypot(K(:,1) - K(orden(1),1), K(:,2) - K(orden(1),2));
    for j = 2:N
        dmin(orden(1:j-1)) = -Inf;
        [~, orden(j)] = max(dmin);
        dmin = min(dmin, hypot(K(:,1) - K(orden(j),1), K(:,2) - K(orden(j),2)));
    end
    T.fps_rank = zeros(N, 1); T.fps_rank(orden) = (1:N)';
    writetable(T, fullfile(sel_dir, 'illumination_positions_all.csv'));

    S = T(orden(1:n_sel), :);
    S.new_index = (0:n_sel-1)';
    S.new_filename = compose("i%03d.pgm", S.new_index);
    writetable(S, fullfile(sel_dir, 'illumination_positions.csv'));
    illum = table2struct(S, 'ToScalar', true);
    illum.filename = cellstr(S.filename); illum.new_filename = cellstr(S.new_filename);
    meta = struct('view', vista, 'folder', carpeta, 'n_selected', n_sel, 'n_total', N, ...
                  'tukey', tukey, 'k_axis_convention', 'RAW: +dcol -> +kx, +drow -> +ky');
    save(fullfile(sel_dir, 'illumination_positions.mat'), 'illum', 'meta', '-v7');

    fid = fopen(fullfile(sel_dir, 'seleccion.txt'), 'w');
    [~, nombre] = fileparts(carpeta);
    fprintf(fid, '# %s: %d de %d hologramas, orden FPS en (kx,ky) RAW\n', nombre, n_sel, N);
    fprintf(fid, '# nuevo_nombre  archivo_original  kx_norm_raw  ky_norm_raw  theta_deg\n');
    for j = 1:n_sel
        fprintf(fid, '%s  %s  %+.5f  %+.5f  %6.2f\n', S.new_filename(j), S.filename(j), ...
                S.kx_norm_raw(j), S.ky_norm_raw(j), S.theta_deg(j));
    end
    fclose(fid);

    rs = fullfile(out, 'real_subset');
    if exist(rs, 'dir'), rmdir(rs, 's'); end
    mkdir(rs);
    for j = 1:n_sel
        copyfile(fullfile(carpeta, S.filename(j)), fullfile(rs, S.new_filename(j)));
    end
    for f = {'Intensite_ref.pgm', 'recon.txt'}
        if isfile(fullfile(carpeta, f{1})), copyfile(fullfile(carpeta, f{1}), rs); end
    end
    txt = regexprep(fileread(fullfile(carpeta, 'config_manip.txt')), '(?m)^NB_HOLO\s+\d+', ...
                    sprintf('NB_HOLO %d', n_sel));
    fid = fopen(fullfile(rs, 'config_manip.txt'), 'w'); fwrite(fid, txt); fclose(fid);
    copyfile(fullfile(sel_dir, 'seleccion.txt'), rs);
    fprintf('%s: %d de %d iluminaciones (theta %.1f-%.1f deg) -> %s\n', vista, n_sel, N, ...
            min(S.theta_deg), max(S.theta_deg), out);
end
end

function p = ruta(raiz, p)
if ~startsWith(p, filesep), p = fullfile(raiz, p); end
end

function cfg = leer_config(archivo)
cfg = struct();
for linea = splitlines(string(fileread(archivo)))'
    p = split(strtrim(extractBefore(linea + "#", "#")));
    if numel(p) >= 2 && ~isnan(str2double(p(2)))
        cfg.(matlab.lang.makeValidName(upper(p(1)))) = str2double(p(2));
    end
end
end

function [sx, sy] = coord_shift(cx, cy, nx, ny)
if cx - nx/2 > 0, sx = cx - nx/2; else, sx = nx/2 + cx; end
if cy - ny/2 > 0, sy = cy - ny/2; else, sy = ny/2 + cy; end
sx = round(sx); sy = round(sy);
end

function w = tukeywin_local(N, r)
x = (0:N-1).' / (N-1); w = ones(N, 1);
i1 = x < r/2;   w(i1) = 0.5 * (1 + cos(2*pi/r * (x(i1) - r/2)));
i2 = x > 1-r/2; w(i2) = 0.5 * (1 + cos(2*pi/r * (x(i2) - 1 + r/2)));
end

function M = centroide(P, pr, pc, R)
[Nr, Nc] = size(P);
r1 = max(1, pr-R); r2 = min(Nr, pr+R); c1 = max(1, pc-R); c2 = min(Nc, pc+R);
L = P(r1:r2, c1:c2); [C, Rr] = meshgrid(c1:c2, r1:r2);
b = false(size(L)); b([1 end], :) = true; b(:, [1 end]) = true;
M.bg = median(L(b)); Wt = max(L - M.bg, 0); s = sum(Wt(:));
if s <= eps, M.r = pr; M.c = pc;
else, M.r = sum(Rr(:) .* Wt(:)) / s; M.c = sum(C(:) .* Wt(:)) / s; end
end
