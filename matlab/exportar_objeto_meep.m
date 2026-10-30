% FDTD_meep v1.0.0 | 2026-10-02 | Carlos Chacón
function info = exportar_objeto_meep(V_yxz, voxel_um, n_medium, h5file, es_delta_n)

if nargin < 5, es_delta_n = true; end
assert(ndims(V_yxz) == 3 && all(isfinite(V_yxz(:))), 'V_yxz debe ser un volumen 3D finito V(y,x,z).');
if isscalar(voxel_um), voxel_um = voxel_um * [1 1 1]; end
if es_delta_n, n = n_medium + double(V_yxz); else, n = double(V_yxz); end
assert(all(n(:) > 0), 'El índice debe ser positivo.');

[Ny, Nx, Nz] = size(n);
borde = [reshape(n([1 end], :, :), [], 1); reshape(n(:, [1 end], :), [], 1); reshape(n(:, :, [1 end]), [], 1)];
salto = max(abs(borde - n_medium));
if salto > 1e-3
    warning('El borde del volumen difiere del medio en %.4g: añade relleno.', salto);
end

if isfile(h5file), delete(h5file); end
h5create(h5file, '/n', [Nx Ny Nz], 'Datatype', 'single', 'ChunkSize', min([Nx Ny Nz], [64 64 16]), 'Deflate', 4);
h5write(h5file, '/n', single(permute(n, [2 1 3])));
eje = {'x', Nx, voxel_um(1); 'y', Ny, voxel_um(2); 'z', Nz, voxel_um(3)};
for i = 1:3
    nombre = sprintf('/coords/%s_um', eje{i, 1});
    h5create(h5file, nombre, eje{i, 2});
    h5write(h5file, nombre, ((0:eje{i, 2}-1) - (eje{i, 2}-1)/2) * eje{i, 3});
end
h5writeatt(h5file, '/', 'axis_order', 'zyx');
h5writeatt(h5file, '/', 'n_medium', n_medium);
h5writeatt(h5file, '/', 'voxel_um_xyz', voxel_um(:)');

info = struct('h5file', h5file, 'size_yxz', [Ny Nx Nz], ...
              'extent_um_xyz', [Nx Ny Nz] .* voxel_um(:)', 'n_min', min(n(:)), 'n_max', max(n(:)));
fprintf('%s: %dx%dx%d voxeles (y,x,z), %.2fx%.2fx%.2f um, n=[%.4f, %.4f]\n', h5file, Ny, Nx, Nz, ...
        info.extent_um_xyz, info.n_min, info.n_max);
end
