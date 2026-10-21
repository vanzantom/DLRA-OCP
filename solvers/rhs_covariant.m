function rhs = rhs_covariant(Abasis, chi, M, Xi,f,dt,Y,Yd,z,i)


if isvector(Xi), w = full(Xi(:)); else, w = full(diag(Xi)); end
        
N=length(w);

if isstruct(Y)==0
    [U,S,V]=svd(Y);
    clear Y;
    Y.U=U; Y.S=S; Y.V=V;
end
GPhi = [Yd, -Y.U*Y.S];
GPsi=[ones(N,1),Y.V]; 

if isempty(f)
    gstate.Phi = z;
else
    gstate.Phi= z + f((i-1)*dt);
end
gstate.Psi= ones(N,1);


rhs = @(Q) sum_tan(rhs1(Q,Y,Abasis,chi,GPhi,GPsi,M,Xi),rhs2(Q,Y,gstate,Abasis,chi,M,Xi));


end

function out=rhs1(Q,Y,Abasis,chi,GPhi,GPsi,M,Xi)

AQ=applyA_tangent(Q, Abasis, chi);
AQ.U=-(M\AQ.U);
Z.U = [AQ.U, GPhi];
Z.V = [AQ.V, GPsi];
Z.S = blkdiag(AQ.S, eye(size(GPhi,2)));

outtemp=projection(Y, Z, M, Xi);
out.dU=outtemp.Up;
out.dS=outtemp.M;
out.dV=outtemp.Vp;

end


function output=sum_tan(rhs1,rhs2)

output.dS=rhs1.dS;
output.dU=rhs1.dU+rhs2.dU;
output.dV=rhs1.dV+rhs2.dV;

end

function output=rhs2(Q,Y,Gstate,Abasis,chi,M,Xi)
S1=Y.S;
dU1=Q.dU;
dV1=Q.dV;
w=diag(Xi);
AY = applyA(Y, Abasis, chi);
                                             % approximate ambient F at Y as factors (use negative M\applyA plus forcing part
UF_at = -(M \ AY.U);
SF_at = AY.S;
VF_at = AY.V;
UF = [UF_at, Gstate.Phi]; VF = [VF_at, Gstate.Psi]; SF = blkdiag(SF_at, eye(size(Gstate.Phi,2)));
                                            % first term is F\perp \Xi Vp_tilde S^(-t)
                                            % compute W = Xi * dV1 (recall Vp_tilde is dV1)
W = bsxfun(@times, dV1, w);
                                            % F_Y * W via low-r factors
F_Y_W = UF * (SF * (VF' * W));              % projection of F_at_Y onto tangent at Y(i+1)
ZF.U = UF; ZF.S = SF; ZF.V = VF;
Zproj_tmp = projection(Y, ZF, M, Xi);    
PYF_W = Y.U*(Zproj_tmp.M*(Y.V'*W))+ Zproj_tmp.Up*(Y.V'*W) + Y.U*(Zproj_tmp.Vp'*W);
FperpW = F_Y_W - PYF_W;                    % F_perp * W    


                                            % regularized S^{-T} from S1
%rcond_thr = getfield_or(opts, 'rcond', 1e-10);
rcond_thr=1e-10;
[Us,Ss,Vs] = svd(S1);
svals = diag(Ss);
cutoff = rcond_thr * max(svals);
inv_s = zeros(size(svals));
inv_s(svals > cutoff) = 1 ./ svals(svals > cutoff);
S_invT = Us * diag(inv_s) * Vs';            % regularized S^{-⊤}
                                            % first term: Uhat = (F_perp * (Xi*dV1)) * S^{-T}
output.dU = FperpW * S_invT;
                                            % second term is
                                            % US^(-t)U_p^\top M F_perp
                                            % compute dU1' * (M * F_perp) via low-r products
                                            % dU1'*(M*F_Y) = (dU1'*(M*UF)) * (SF * VF')                                     
dU_M_FY = (dU1'*(M*UF))*(SF*VF');                              
dU_M_PYF=(dU1'*(M*Y.U))*(Zproj_tmp.M*Y.V')+(dU1'*(M*Zproj_tmp.Up))*(Y.V')+(dU1'*M*Y.U)*(Zproj_tmp.Vp');
dUMF = dU_M_FY - dU_M_PYF;                       

output.dV = (S_invT * dUMF)';

end

function Astruct = applyA_tangent(Q, Abasis, chi)
    % Apply A to  Q = U*dS*V' + dU*V' + U*dV'
    % A(X) = A0*X + sum_k Ak*X*Dk, Dk = diag(chi(:,k))
    % Returns Astruct such that
    % A(Q) = Astruct.U * Astruct.S * Astruct.V'

    U  = Q.Ubase;
    dU = Q.dU;
    dS = Q.dS;
    V  = Q.Vbase;
    dV = Q.dV;

    p = size(chi,2);

    % -------------------------------------------------------------
    % First contribution: U*dS*V'
    % -------------------------------------------------------------

    Ublocks1 = cell(p+1,1);
    Vblocks1 = cell(p+1,1);

    Ublocks1{1} = Abasis{1} * U;
    Vblocks1{1} = V;

    for kk = 1:p
        Ublocks1{kk+1} = Abasis{kk+1} * U;
        Vblocks1{kk+1} = bsxfun(@times,V,chi(:,kk));
    end

    % -------------------------------------------------------------
    % Second contribution: dU*V'
    % -------------------------------------------------------------

    Ublocks2 = cell(p+1,1);
    Vblocks2 = cell(p+1,1);

    Ublocks2{1} = Abasis{1} * dU;
    Vblocks2{1} = V;

    for kk = 1:p
        Ublocks2{kk+1} = Abasis{kk+1} * dU;
        Vblocks2{kk+1} = bsxfun(@times,V,chi(:,kk));
    end

    % -------------------------------------------------------------
    % Third contribution: U*dV'
    % -------------------------------------------------------------

    Ublocks3 = cell(p+1,1);
    Vblocks3 = cell(p+1,1);

    Ublocks3{1} = Abasis{1} * U;
    Vblocks3{1} = dV;

    for kk = 1:p
        Ublocks3{kk+1} = Abasis{kk+1} * U;
        Vblocks3{kk+1} = bsxfun(@times,dV,chi(:,kk));
    end

    % -------------------------------------------------------------
    % Combine the three contributions
    % -------------------------------------------------------------

    Uall = [Ublocks1{:}, Ublocks2{:}, Ublocks3{:}];
    Vall = [Vblocks1{:}, Vblocks2{:}, Vblocks3{:}];

    r = size(U,2);

    Sall = blkdiag( ...
        kron(eye(p+1), dS), ...
        kron(eye(p+1), eye(r)), ...
        kron(eye(p+1), eye(r)) );

    Astruct.U = Uall;
    Astruct.S = Sall;
    Astruct.V = Vall;
end
    
function Astruct = applyA(X, Abasis, chi)
        % Return low-r factors for applyA(X) = A0*X + sum_k A_k * X * Dk
        Ux = X.U; Sx = X.S; Vx = X.V;
        p = size(chi,2);
        Ublocks = cell(p+1,1);
        Vblocks = cell(p+1,1);
        Ublocks{1} = Abasis{1} * Ux;
        Vblocks{1} = Vx;
        for kk = 1:p
            Ublocks{kk+1} = Abasis{kk+1} * Ux;
            Vblocks{kk+1} = bsxfun(@times, Vx, chi(:,kk));
        end
        Uall = [Ublocks{:}];
        Vall = [Vblocks{:}];
        Sall = kron(eye(p+1), Sx);
        Astruct.U = Uall;
        Astruct.S = Sall;
        Astruct.V = Vall;
    end