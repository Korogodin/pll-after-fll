classdef CKalmanEqMesConstK < handle
    %CKALMAN 
    
    properties
        %Scalar
        Band
        nx
        T
        
        %Vectors
        Xest
        Xextr
        Koeff
        
        % Matrix
        F
    end
    
    methods
        function K = CKalmanEqMesConstK(Xest0, F, T)
            K.Xest = Xest0;
            K.F = F; 
            K.nx = length(Xest0);
            K.T = T;
        end
        
        function Extrapolate(K)
            K.Xextr = K.F*K.Xest;
        end
        
        function Estimate(K, u_norm)
            K.Xest = K.Xextr + K.Koeff*u_norm;
        end
        
        function setKoeffFromBand(K, Band)
            K.Band = Band;
            if K.nx == 2
                Kn2 = 32 / 9 * Band^2;
                Kn1 = sqrt(2*Kn2);
                K.Koeff = [Kn1; Kn2] * K.T;
            else
                disp('Need function: setKoeffFromBand');
            end
        end
        
        function setKoeffFromSpectr(K, Sx, Smeas)
            if K.nx == 2
                Kn2 = sqrt(Sx / Smeas);
                Kn1 = sqrt(2*Kn2);
                K.Koeff = [Kn1; Kn2] * K.T;
            elseif K.nx == 3 
                Kn3 = sqrt(Sx / Smeas);
                Kn2 = 2*(Kn3)^(2/3);
                Kn1 = 2*(Kn3)^(1/3);
                K.Koeff = [Kn1; Kn2; Kn3] * K.T;
            else
                disp('Need function: setKoeffFromSpectr');
            end
        end        
        
        function setX(K, X)
            K.Xest = X;
        end
        
        function setF(K, F)
            K.F = F;
        end
        
        function setKoeff(K, Koeff)
            K.Koeff = Koeff;
        end
        
        function X = getXest(K)
            X = K.Xest;
        end
        
        function calcBand(K)
            T = K.T;
            
            K1 = K.Koeff(1);
            K2 = K.Koeff(2);
            if K.nx == 3
                K3 = K.Koeff(3);
            else
                K3 = 0;
            end
            K.Band = (4*K1^2*K2*T - 6*K1*K2^2*T^2 + 7*K1*K2*K3*T^3 - 2*K1*K3^2*T^4 - 4*K1*K3*T^2 + 4*K2^2*T^2 - 4*K2*K3*T^3 + K3^2*T^4) / ...
                            ( (K1*K2*T - K3*T^2) * (K3*T^2 - 2*K2*T - 4*K1 + 8) * 2 * T );
        end
    end
end

