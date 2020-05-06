
d_0           = 1e-4         #
L             = 4e-4         # Domain radial length
N_p           = 64           #
carrier_subst = "air"        #
fuel_subst    = "acetone"    #
p_atm         = 101325       #
T_0_L         = 294.92       #
T_infty       = 700          #
Y_F_infty     = 0            #
Fo_crit       = 0.4          # Fo_crit<0.5
CFL_crit      = 0.8          # CFL_crit<1.0

function INIT(d_0,L,N_p,carrier_subst,fuel_subst,p_atm,T_0_L,T_infty,Y_F_infty)
end

function molarMass(subst)
#Tabulates the Molar Mass of given substance
    if subst=="air"
        Mw = 28.9647
    elseif subst=="acetone"
        Mw = 58.08
    end
    return Mw
end

function boilingTemperature(fuel_subst)
#Tabulates the Boiling Temperature of given fuel_substance
    if fuel_subst=="acetone"
        T_b = 329
    end
    return T_b
end

function heatVaporization(method,fuel_subst)
#Computes the heat of Vaporization h_vap of given fuel_substance
    if method=="CC"
    elseif method=="tabulated"
        if fuel_subst=="acetone"
            h_vap = 5.18e5
        end
    end
    return h_vap
end

function density(method,subst)
    if method=="NASA"
    elseif method=="tabulated"
        if subst=="air"
            ρ = 1
        elseif subst=="acetone"
            ρ = 700
    end
    return ρ
end

function heatConductivity(method,subst)
    if method=="NASA"
    elseif method=="tabulated"
        if subst=="air"
            λ = 5.2e-2
        elseif subst=="acetone"
            λ = 1.61e-1
        end
    end
    return λ
end

function binaryDiffusionCoef(method,subst_1,subst_2)
    if method=="NASA"
    elseif method=="tabulated"
        if subst_1=="air" && subst_2=="acetone"
            D = 5.2e-5
        end
    end
    return D
end

function specificHeat(method,subst)
    if method=="NASA"
    elseif method=="tabulated"
        if subst=="air"
            c_p = 1000
        elseif subst=="acetone"
            c_p = 2000
        end
    end
    return c_p
end


α_G      = λ_G/(ρ_G*c_p_G)
α_L      = λ_L/(ρ_L*c_p_L)
r_int    = d_0/2
r_int_0  = r_int
dr       = L/(N_p-1)
dt1      = Fo_crit*dr^2/α_L
dt2      = Fo_crit*dr^2/α_G
dt3      = Fo_crit*dr^2/D
dt       = min(dt1,dt2,dt3)
N_c      = N_p-1
p_grid   = collect(0:dr:N_p)
c_grid   = collect(0:dr/2:N_c)
u_r      = zeros(1,N_p)
N_p_L    = Int(floor(r_int/dr)+1)
N_c_L    = N_p_L-1
theta    = ((N_p_L+1)*dr-r_int)/dr

λ_field   = Array{Float64}(undef,N_c)
ρ_field   = Array{Float64}(undef,N_c)
c_ρ_field = Array{Float64}(undef,N_c)

mutable struct paramStruct
    λ_L,λ_cell,λ_G::Float64
    ρ_L,ρ_cell,ρ_G::Float64
    c_p_L,c_p_cell,c_p_G::Float64
    D::Float64
    N_c_L_old,N_c_L,N_p_L,N_p::Int64
    theta::Float64
    dr::Float64
end

#Terminar de fazer o mapeamento de paramStruct vers os parametros
param = paramStruct()
param.λ_cell   = cellVolAvg(λ_L,λ_G,(N_p_L-1)*dr,(N_p_L-theta)*dr,N_p_L*dr)
param.ρ_cell   = cellVolAvg(ρ_L,ρ_G,(N_p_L-1)*dr,(N_p_L-theta)*dr,N_p_L*dr)
param.c_p_cell = cellVolAvg(c_p_L,c_p_G,(N_p_L-1)*dr,(N_p_L-theta)*dr,N_p_L*dr)

for j in 1:N_c
    if j<=N_c_L
        λ_field[j]   = λ_L
        ρ_field[j]   = ρ_L
        c_p_field[j] = c_p_L
    elseif j==N_c_L+1
        λ_field[j]   = λ_cell
        ρ_field[j]   = ρ_cell
        c_p_field[j] = c_p_cell
    else
        λ_field[j]   = λ_G
        ρ_field[j]   = ρ_G
        c_p_field[j] = c_p_G
    end
end

function cellVolAvg(qt_L,qt_G,r1,r2,r3)
    cellVolAvg = (4/3*π*(r3^3-r2^3)*qt_G+4/3*π*(r2^3-r1^3)*qt_L)/(4/3*π*(r3^3-r1^3))
end

function propFieldUpdate!(λ_field::Array{Float64},ρ_field::Array{Float64},c_p_field::Array{Float64},param::paramStruct)
    if N_c_L_old ~= N_c_L
        λ_field[N_c_L+2]   = param.λ_G
        ρ_field[N_c_L+2]   = param.ρ_G
        c_p_field[N_c_L+2] = param.c_p_G
    end
    λ_field[N_c_L+1]   = param.λ_cell
    ρ_field[N_c_L+1]   = param.ρ_cell
    c_p_field[N_c_L+1] = param.c_p_cell

    return λ_field, ρ_field, c_p_field
end

λ_field, ρ_field, c_p_field = propFieldUpdate(λ_L,λ_cell,ρ_L,c_p_L,λ_G,ρ_G,c_p_G,D,N_c_L,N_p_L,theta,dr,N_p)
