classdef Plato < handle

    properties
        nro
        z_i
        x_i
        y_i
        v_i
        l_i
        T
        P
        V
        L
        K
        num_sust
        aliment = Corriente.empty(0,0) % tipo: 'Corriente'
        salidaV = 0;
        salidaL = 0;
        pumparoundOut
        pumparoundIn
    end
    
    methods 
        function self = Plato(nro, T, P, v_i, l_i, x_i, y_i, V, L, K, alim, saleV, saleL, pumpout, pumpin)
            if nargin > 0 && ~isempty(nro)
                self.nro = nro;
            end
            if nargin > 1 && ~isempty(T)
                self.T = T;
            end
            if nargin > 2 && ~isempty(P)
                self.P = P;
            end
            if nargin > 3 && ~isempty(v_i)
                self.v_i = v_i;
            end
            if nargin > 4 && ~isempty(l_i)
                self.l_i = l_i;
            end
            if nargin > 5 && ~isempty(x_i) 
                self.x_i = x_i;
            end
            if nargin > 6 && ~isempty(y_i)
                self.y_i = y_i;
            end
            
            if nargin > 10 && ~isempty(alim)
                self.aliment = alim;
            end
            if nargin > 11 && ~isempty(saleV)
                self.salidaV = saleV;
            end
            if nargin > 12 && ~isempty(saleL)
                self.salidaL = saleL;
            end
            if nargin > 7 && ~isempty(V)
                self.V = V;
                if isempty(v_i) && ~isempty(y_i)
                    self.v_i = self.V.*self.y_i;
                end
            end
            if nargin > 8 && ~isempty(L)
                self.L = L;
                if isempty(l_i) && ~isempty(x_i)
                    self.l_i = self.L.*self.x_i;
                end
            end
            if nargin > 9 && ~isempty(K)
                self.K = K;
            end
            if nargin > 13 && ~isempty(pumpout)
                self.pumparoundOut = pumpout;
            end
            if nargin > 14 && ~isempty(pumpin)
                self.pumparoundIn = pumpin;
            end
            try
                self.mixtured();
            catch
            end
        end
        function self = setT(self, T)
            if nargin > 1
                if ~isempty(T)
                    self.T = T;
                end
            end
        end
        function self = setP(self,P)
            if nargin > 1 
                if ~isempty(P)     
                    self.P = P;
                end
            end
        end
        function self = setvi(self, vi)
            if nargin > 1
                if ~isempty(vi)
                    self.v_i = vi;
                    self.y_i = vi./sum(vi);
                    self.V = sum(vi);
                end  
                if isempty(self.num_sust)
                    [self.num_sust] = length(vi);
                end
                try
                    self.mixtured();
                catch
                end
            end
        end
        function self = setli(self, li)
            if nargin > 1
                if ~isempty(li)
                    self.l_i = li;
                    self.x_i = li/sum(li);
                    self.L = sum(li);
                
                end
                if isempty(self.num_sust)
                    [self.num_sust] = length(li);
                end
                try 
                    self.mixtured();
                catch
                end
            end
        end
        function self = setxi(self, xi)
            if nargin > 1
                if ~isempty(xi)
                    self.x_i = xi;
                    self.l_i = self.L.*xi;
                end
                try
                    self.mixtured();
                catch
                end
            end
        end
        function self = setyi(self, yi)
            if nargin > 1
                if ~isempty(yi)
                    self.y_i = yi;
                    self.v_i = yi.*self.V;
                end
                try
                    self.mixtured();
                catch
                end
            end
        end
        function self = setV(self, V)
            if nargin > 1
                if ~isempty(V)
                    self.V = V;
                    self.v_i = self.y_i.*self.V;
                end
                try 
                    self.mixtured();
                catch
                end
            end
        end
        function self = setL(self, L)
            if nargin > 1
                if ~isempty(L)
                    self.L = L;
                    self.l_i = self.x_i.*self.L;
                end
                try
                    self.mixtured();
                catch
                end
            end
        end
        function self = mixtured(self)
            self.z_i = (self.L.*self.x_i + self.V.* self.y_i)/(self.V+self.L);
        end
    end
end