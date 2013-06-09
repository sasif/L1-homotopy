% script to update support
if length(gamma_xh)==M && flag == 1
    % swap an element... 
    Ag = A(:,gamma_xh)'*A(:,idelta);
    xg = xk_1(gamma_xh);
    zx = -sign(pk(idelta)+delta*dk(idelta));
    
    % delx = -AtAgx\(Ag*zx);
    % delx = -AtAgx\(Ag*zx); % WHY?????????
    switch delx_mode
        case 'mil'
            delx = -iAtA*(Ag*zx);
        case 'qr'
            delx = -R\(R'\(Ag*zx));
    end
    
    delta_out = -xg./delx;
    idelta_out = find(delta_out>0);
    if ~isempty(idelta_out)
        [delta_t idelta_tmp] = min(delta_out(idelta_out));
        idelta_t = gamma_xh(idelta_out(idelta_tmp));
        
        x_old = xk_1;
        xk_1(gamma_xh) = xg+delta_t*delx;
        xk_1(idelta) = delta_t*zx;
        xk_1(idelta_t) = 0;
        
        % Remember remember the Simplex pivoting
        % disp(sprintf(['iter = %d, delta = %3.4g, delta_t = %3.4g, idelta = %d, idelta_t = %d, flag = %d'], iter, delta, delta_t, idelta, idelta_t, -2));
        
        outx_index = find(gamma_xh==idelta_t);
        gamma_xh = [gamma_xh(1:outx_index-1); gamma_xh(outx_index+1:end); gamma_xh(outx_index)];
        g_old = gamma_xh;
        
        gamma_xh(end) = idelta;
        
        % itr_history = [itr_history; idelta_t delta_t 2];
        
        new_order = [1:outx_index-1, outx_index+1:length(g_old), outx_index]';
        
        in_delx.new_order = new_order;
        in_delx.out_i = outx_index;
        in_delx.add_rem = 2;
        in_delx.Gamma = g_old;
        in_delx.new_x = idelta;
        idelta = idelta_t; % To set its constraints.
    else
        disp('Fatal error');
        break;
    end
else
    if flag == 1
        % add an element
        g_old = gamma_xh;
        gamma_xh = [gamma_xh; idelta];
        xk_1(idelta) = 0;
        
        gamma_xc(gamma_xc==idelta) = [];
        
        in_delx.new_x = idelta;
        in_delx.add_rem = 1;
        in_delx.Gamma = g_old;
    else
        % remove an element
        outx_index = find(gamma_xh==idelta);
        len_gamma = length(gamma_xh);
        if len_gamma == 1
            break;
        end
        gamma_xh = [gamma_xh(1:outx_index-1); gamma_xh(outx_index+1:end); gamma_xh(outx_index)];
        g_old = gamma_xh;
        
        gamma_xh = gamma_xh(1:len_gamma-1);
        
        gamma_xc = [gamma_xc; idelta];
        
        xk_1(idelta) = 0;
        
        new_order = [1:outx_index-1, outx_index+1:length(g_old), outx_index]';
        
        in_delx.new_order = new_order;
        in_delx.out_i = outx_index;
        in_delx.add_rem = 0;
        in_delx.Gamma = g_old;
    end
end

% temp_gamma = zeros(N,1);
% temp_gamma(gamma_xh) = gamma_xh;
% gamma_xc = find([1:N]' ~= temp_gamma);
