
d_0           = 1e-4
L             = 4e-4         # Domain radial length
N_p           = 64           #
carrier_subst = "air"        #
fuel_subst    = "water"      #
p_atm         = 101325       #
T_0_L         = 294.92       #
T_0_G         = 700          #
T_infty       = 700          #
Y_F_infty     = 0            #

W_F      = 58
W_G      = 29
T_b      = 329
h_vap    = 5.18e5
ρ_G      = 1
ρ_L      = 700
λ_L      = 1.61e-1
λ_G      = 5.2e-2
D        = 5.2e-5
c_p_L    = 2000
c_p_G    = 1000
α_G      = λ_G/(ρ_G*c_p_G)
α_L      = λ_L/(ρ_L*c_p_L)

r_int    = d_0/2
r_int_0  = r_int
dr       = L/(N_p-1)
Fo_crit  = 0.4
CFL_crit = 0.8
dt1      = Fo_crit*dr^2/α_L
dt2      = Fo_crit*dr^2/α_G
dt3      = Fo_crit*dr^2/D
dt       = min(dt1,dt2,dt3)
x = collect(0:dx:1)
N_c      = N_p-1
p_grid   = [dr*i for i in 0:N_p-1]
c_grid   = [dr/2*j for j in 1:N_c]
u_r      = zeros(1,N_p)
N_p_L    = Int(floor(r_int/dr)+1)
N_c_L    = N_p_L-1
theta    = ((N_p_L+1)*dr-r_int)/dr

function volAvg(qt_L,qt_G,r1,r2,r3)
    volAvg = (4/3*π*(r3^3-r2^3)*qt_G+4/3*π*(r2^3-r1^3)*qt_L)/(4/3*π*(r3^3-r1^3))
end

λ_cell   = volAvg(λ_L,λ_G,(N_p_L-1)*dr,(N_p_L-theta)*dr,N_p_L*dr)
ρ_cell   = volAvg(ρ_L,ρ_G,(N_p_L-1)*dr,(N_p_L-theta)*dr,N_p_L*dr)
c_p_cell = volAvg(c_p_L,c_p_G,(N_p_L-1)*dr,(N_p_L-theta)*dr,N_p_L*dr)

λ_field  = ones(N_c)
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
#Abaixo como modificar o struct
param.λ_cell = volAvg(λ_L,λ_G,(N_p_L-1)*dr,(N_p_L-theta)*dr,N_p_L*dr)

λ_field   = Array{Float64}(undef,N_c)
ρ_field   = Array{Float64}(undef,N_c)
c_ρ_field = Array{Float64}(undef,N_c)
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

function propFieldUpdate!(λ_field,ρ_field::Array{Float64},c_p_field::Array{Float64},hello::paramStruct)
    λ_field::Array{Float64}
    if N_c_L_old ~= N_c_L
        λ_field[N_c_L+2]   = λ_G
        ρ_field[N_c_L+2]   = ρ_G
        c_p_field[N_c_L+2] = c_p_G
    end
    λ_field[N_c_L+1]   = λ_cell
    ρ_field[N_c_L+1]   = ρ_cell
    c_p_field[N_c_L+1] = c_p_cell

    return λ_field, ρ_field, c_p_field
end

λ_field, ρ_field, c_p_field = propFieldUpdate(λ_L,λ_cell,ρ_L,c_p_L,λ_G,ρ_G,c_p_G,D,N_c_L,N_p_L,theta,dr,N_p)
