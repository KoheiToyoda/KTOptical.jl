module KTOptical
    using Parameters
# Write your package code here.

    # l,p: LG のモードパラメータ
    # x,y: 座標(使用する際はブロードキャストする)
    # ω,λ: モードフィールド径、波長

    # ラゲール倍多項式
    @with_kw struct struct_beamparam
        ω::Float64
        z::Float64
        λ::Float64
        k::Float64
        zr::Float64
        ωz::Float64
        Rz::Float64
    end

    export LG_E
    export LG_I
    export HG_E
    export HG_I

    function setParam(ω_,z_,λ_)
        # 計算結果を使ってパラメータセットするときは、
        # 引数内の宣言内のスコープが引き継がれない。
        # よってω、ｚ、λを使った計算はk,zr,wz,rzに入れられない。
        # ω_、ｚ_、λ_を使う必要がある。
        k = 2π / λ_
        zr = k * ω_^2 / 2
        ωz = ω_ * √(1 + (z_^2 / zr^2))
        Rz = k * ω_^2 / 2

        global bp = struct_beamparam(ω_,z_,λ_,k,zr,ωz,Rz)
    end

    function __Llp__(l,p,u)
        function Llp_setinit(l,p,u,m)
            if m < 0 #再起関数の終了条件
                return 1
            else
                return Llp_setinit(l,p,u,m-1) + (-1)^m*(2p+abs(l)-m)* u^m / factorial(m)
            end
        end
        return Llp_setinit(l,p,u,p)
    end

    function LG_E(l,p,x,y)
        R = √(x^2 + y^2)
        Θ = atan(y/x)
        E = @. √(2 * factorial(p) / (π * factorial(p + abs(l)))) *
            (√2R / bp.ω)^abs(l) *
            __Llp__(l, p, 2 * R^2 / bp.ω^2) *
            (bp.ω / bp.ωz) *
            exp(-R^2 / bp.ω^2) *
            exp(1im * l * Θ) *
            exp(-1im * (1 + 2p + abs(l)) * atan(bp.z / bp.zr)) *
            exp(-1im * bp.k * R^2 / (2bp.Rz))
        return E
    end

    function LG_E(l,p,x)
        R = √(x^2)
        Θ = atan(0)
        E = @. √(2 * factorial(p) / (π * factorial(p + abs(l)))) *
            (√2R / bp.ω)^abs(l) *
            __Llp__(l, p, 2 * R^2 / bp.ω^2) *
            (bp.ω / bp.ωz) *
            exp(-R^2 / bp.ω^2) *
            exp(1im * l * Θ) *
            exp(-1im * (1 + 2p + abs(l)) * atan(bp.z / bp.zr)) *
            exp(-1im * bp.k * R^2 / (2bp.Rz))
        return E
    end

    function LG_I(l,p,x,y)
        return abs(LG_E(l,p,x,y))^2
    end

# m,n: HG のモードパラメータ
# x,y: 座標(使用する際はブロードキャストする)

# ラゲール倍多項式
    function __Hn__(n,u)
        function Hn_setinit(n,u,m)
            if m > floor(n/2) #再起関数の終了条件
                return 0
            else
                return Hn_setinit(n, u, m + 1) + (-1)^m/(factorial(m)*factorial(n-2m))*(2u)^(n-2m)
            end
        end
        return factorial(n) * Hn_setinit(n,u,0)
    end

    function HG_E(m,n,x,y)
        E = __Hn__(n, √2 * y / bp.ω) *
        __Hn__(m, √2 * x / bp.ω) *
        exp(-(x^2 + y^2) / bp.ω^2) *
        exp(-1im * (1 + m + n) * atan(bp.z / bp.zr)) *
        exp(-1im * bp.k * (x^2 + y^2) / (2bp.Rz))
        return E
    end

    function HG_E(m,n,x)
        E = __Hn__(m, √2 * x / bp.ω) *
        exp(-(x^2) / bp.ω^2) *
        exp(-1im * (1 + m + n) * atan(bp.z / bp.zr)) *
        exp(-1im * bp.k * (x^2) / (2bp.Rz))
        return E
    end

    function HG_I(m,n,x,y)
        return abs(HG_E(m,n,x,y))^2
    end

    function HG_I(m,n,x)
        return abs(HG_E(m,n,x))^2
    end

end

