function Zproj = projection(X,Z,K,D)
%PROJECTION computes the orthogonal projection of Z on the tangent space at X.

    ZDV = apply_ambient(Z, D*X.V);
    UtKZDV = X.U'*(K*ZDV);
    ZtKU = apply_ambient_transpose(Z, K*X.U);

    Zproj.M = UtKZDV;
    Zproj.Up = ZDV  - X.U*UtKZDV;
    Zproj.Vp = ZtKU - X.V*UtKZDV';

end

 function ZW = apply_ambient(Z, W)
        if ~isstruct(Z)
            ZW = Z*W;
        else
            ZW = Z.U*(Z.S*(Z.V'*W));
        end
 end

 function ZtW = apply_ambient_transpose(Z, W)
        if ~isstruct(Z)
            ZtW = Z'*W;
        else
            ZtW = Z.V*(Z.S'*(Z.U'*W));
        end
    end