% script to ubpdate delx
switch update_mode
    case 'init0'
        switch delx_mode
            case 'mil'
                iAtA = pinv(A(:,gamma_xh)'*A(:,gamma_xh));
                delx = iAtA*rhs(gamma_xh);
            case {'qr','chol'}
                [Q R] = qr(A(:,gamma_xh),0);
                delx = R\(R'\rhs(gamma_xh));
            case 'qrM'
                [Q0 R0] = qr(A(:,gamma_xh));
                delx = R0\(R0'\rhs(gamma_xh));            
        end
    case 'init1'
        switch delx_mode
            case 'mil'
                % AtAgx = in.AtA;
                iAtA = opts.iAtA;
                delx = iAtA*rhs(gamma_xh);
                % iAtAgx = opts.iAtAgx;
                % delx_update = iAtAgx*(ds(gamma_xh)-d(gamma_xh));
            case 'qr'
                Q = opts.Q;
                R = opts.R;
                delx = R\(R'\rhs(gamma_xh));
                % delx_update = R\(R'\(ds(gamma_xh)-d(gamma_xh)));
            case 'qrM'
                Q0 = opts.Q0;
                R0 = opts.R0;
                delx = R0\(R0'\rhs(gamma_xh));            
                % delx_update = R0\(R0'\(ds(gamma_xh)-d(gamma_xh)));
            case 'chol'
                R = opts.R;
                delx = R\(R'\rhs(gamma_xh));
        end
        % delx = opts.delx;
    case 'update'
        % script to update delx
        in_delx.rhs = rhs;
        in_delx.itr = iter;
        in_delx.M = M;
        switch delx_mode
            case 'mil'
                in_delx.iAtA = iAtA;
                in_delx.delx = delx;
                in_delx.type = delx_mode;
                out = compute_delx(in_delx,A);
                delx = out.delx;
                iAtA = out.iAtA;
            case 'qr'
                in_delx.R = R;
                in_delx.Q = Q;
                in_delx.delx = delx;
                in_delx.type = delx_mode;
                in_delx.nre = 2;
                out = compute_delx(in_delx,A);
                R = out.R;
                Q = out.Q;
                delx = out.delx;
            case 'qrM'
                in_delx.R0 = R0;
                in_delx.Q0 = Q0;
                in_delx.delx = gamma_xh;
                in_delx.type = 'qrM';
                out = compute_delx(in_delx,A);
                R0 = out.R0;
                Q0 = out.Q0;
                delx = out.delx;
            case 'chol'
                in_delx.R = R;
                in_delx.delx = delx;
                in_delx.type = delx_mode;
                out = compute_delx(in_delx,A);
                R = out.R;
                delx = out.delx;
            case 'cg'
                % NOT TESTED ... 
                % delx_orig = (A(:,gamma_xh)'*A(:,gamma_xh))\(ds(gamma_xh)-d(gamma_xh));
                in_delx.x0 = xk_1;
                in_delx.Gamma = gamma_xh;
                in_delx.delx = [];
                in_delx.W = in.W_new(gamma_xh);
                in_delx.type = 'cg';
                delta = theta;
                in_delx.rhs = cg_rhs-(epsilon+delta)*d-(1-epsilon-delta)*ds+tau*sign(pk_old);
                out = compute_delx(in_delx,A);
                delx = (out.delx-xk_1(gamma_xh))/delta;
                stp = 1;
        end
    case 'recompute'
        switch delx_mode
            case 'mil'
                delx = iAtA*rhs(gamma_xh);
                % iAtAgx = in.iAtAgx;
                % delx_update = iAtAgx*(ds(gamma_xh)-d(gamma_xh));
            case {'qr','chol'}
                delx = R\(R'\rhs(gamma_xh));
            case 'qrM'
                delx = R0\(R0'\rhs(gamma_xh));
        end
end