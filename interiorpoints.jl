#
# Interior Points
#

# Guilherme Pereira Freire Machado

function interior_points(A::Array{Float64}, b::Array{Float64}, c::Array{Float64})
    err = 1e-3
    α = 0.9

    # 1) initialization
    x, s, p, A, b, c = initialization(A, b, c)
    n = size(A)[2]

    maxit = 100
    k=0
    for i in 1:maxit
        println(i)
        # 2) 1st test for convergence
        if optimality_test(s, x, err)
            # convergiu
            return "gg" # todo
        end
        
        # 3) computation of newton directions
        μ = 1 *  x' * s / n # ρ *
        dx, dp, ds = compute_directions(x, s, p, μ, A, b, c)
        
        # 4) and 5) update variables
        x, s, p = update_variables(x, s, p, dx, dp, ds, α)

        k += 1
    end

end

function initialization(A::Array{Float64,2}, b::Array{Float64}, c::Array{Float64})
    # get maximum value
    mx = [maximum(A), maximum(b), maximum(c)]
    U = maximum(mx)
    M = U*1e5

    # in order to get the center of the polyhedron:

    # analytic center for P

    # analytic center for Q
    μ0 = 4 * (norm(c)^2 + M^2)^0.5

    m, nx = size(A)

    # x
    x1 = ones(nx)
    x2 = 1
    x3 = 1
    x = [x1; x2; x3] 
    
    # p
    p1 = zeros(size(A)[1])
    p2 = - μ0 
    p = [p1; p2]

    # s
    e = ones(length(c))
    s1 = c + μ0*e
    s2 = M + μ0 
    s3 = μ0
    s = [s1; s2; s3]

    # new c
    c_new = [c; M; 0]

    b_bar = (nx + 2) * b / (nx*(m*U)^m)
    
    # new A
    A_new = [A (b_bar - A*ones(nx)) zeros(m); ones(nx)' 1 1] 
    
    # new b
    b_new = [b_bar; nx + 2]

    return x , s, p, A_new, b_new, c_new
end

function interior_phase1(A, b, c)
    m, n = size(A)

    x0 = ones(n)

    v = b - A*x0

    nv = length(v)

    if all(v .== zeros(n))

    end

    # new opt problem
    c_1 = [zeros(n) ; 1]
    A_1 = [A ones(nv)]
    b_1 = b

    # initial guess
    x0 = [v; 1]

    interior_algorithm(A_1, b_1, c_1, x0, s, p)
    # r, z, status = Simplex([A_1 eye(m) ], b, [c_1; zeros(m)])

    # u value found
    u = r[n+1]

    # update v
    v = v * u

    return v
end 


function interior_bigM(A, b, c)
    m, n = size(A)

    x0 = ones(n)

    v = b - A*x0

    nv = length(v)

    if all(v .== zeros(n))

    end

    # new opt problem
    c_1 = [zeros(n) ; 1]
    A_1 = [A ones(nv)]
    b_1 = b

    # initial guess
    x0 = [v; 1]

    interior_algorithm(A_1, b_1, c_1, x0, s, p)
    # r, z, status = Simplex([A_1 eye(m) ], b, [c_1; zeros(m)])

    # u value found
    u = r[n+1]

    # update v
    v = v * u

    return v
end 

function interior_algorithm(A, b, c, x, s, p)
    err = 1e-3
    α = 0.9

    n = size(A)[2]

    maxit = 100
    k=0
    for i in 1:maxit
        println(i)
        # 2) 1st test for convergence
        if optimality_test(s, x, err)
            # convergiu
            println("Interior Points algorithm converged!")
            return x, s, p
        end
        
        # 3) computation of newton directions
        μ = 1 *  x' * s / n # ρ *
        dx, dp, ds = compute_directions(x, s, p, μ, A, b, c)
        
        # 4) and 5) update variables
        x, s, p = update_variables(x, s, p, dx, dp, ds, α)

        k += 1
    end
end

function optimality_test(s::Array{Float64}, x::Array{Float64}, err::Float64)
    return s'*x < err
end

function compute_directions(x::Array{Float64}, s::Array{Float64}, p::Array{Float64}, μ::Float64,  A::Array{Float64}, b::Array{Float64}, c::Array{Float64})
    # solve linear system:
    nx = length(x)
    np = length(p)
    ns = length(s)
    e = ones(nx)

    matrix = [A zeros(size(A)[1], size(A')[2]) zeros(size(A)[1], nx) ;
            zeros(size(A')[1], size(A)[2]) A' eye(size(A')[1]) ;
            diagm(s) zeros(length(x), size(A')[2]) diagm(x)]

    d = matrix \ - [A * x - b ;
                A' * p + s - c;
                diagm(x)*diagm(s)*e - μ * e]
    
    dx = d[1:nx]
    dp = d[(nx + 1) : (nx + np)]
    ds = d[(nx + np + 1) : end]

     return dx, dp, ds
end

function update_variables(x::Array{Float64}, s::Array{Float64}, p::Array{Float64}, dx::Array{Float64}, dp::Array{Float64}, ds::Array{Float64}, α::Float64)
    # compute step β for primal
    ratio_x = -x ./ dx
    ratio_x[dx .>= 0] = Inf
    βp = minimum([1, α*minimum(ratio_x)]) 
    
    # compute step β for dual
    ratio_s = -s ./ ds
    ratio_s[ds .>= 0] = Inf
    βd = minimum([1, α*minimum(ratio_s)]) 
    
    # update variables
    x = x + βp * dx
    s = s + βd * ds
    p = p + βd * dp
    return x, s, p
end

function open_log(A::Array{Float64,2}, b::Array{Float64,1}, c::Array{Float64,1})
    fname = "Simplex.log"
    if isfile(fname)
        stream = open(fname, "a")
        pwrite(stream, "=======================")
        pwrite(stream, "Comeco da Solucao do PL")
        pwrite(stream, "=======================")
        pwrite(stream, "Problema:")
        pwrite(stream, "A = $A")
        pwrite(stream, "b = $b")
        pwrite(stream, "c = $c")
        pwrite(stream, "")
        close(stream)
    else
        stream = open(fname, "w")
        pwrite(stream, "=======================")
        pwrite(stream, "Comeco da Solucao do PL")
        pwrite(stream, "=======================")
        pwrite(stream, "Problema:")
        pwrite(stream, "A = $A")
        pwrite(stream, "b = $b")
        pwrite(stream, "c = $c")
        pwrite(stream, "")
        close(stream)
    end
    nothing
end

function get_log(state::Int)
    fname = "SimplexFase2.log"
    stream = open(fname, "a")
    if state == 1
        pwrite(stream, "Simplex Fase 1")
        pwrite(stream, "--------------")
        
    else
        pwrite(stream, "Simplex Fase 2")
        pwrite(stream, "--------------")

    end
    return stream
end

function result_log(it::Int, x::Array{Float64,1}, bidx::Array{Int,1}, nidx::Array{Int,1}, z::Float64, status::Int, stream::IOStream)
    pwrite(stream, "iter $it:")
    pwrite(stream, "x = $x")
    pwrite(stream, "Base = $bidx")
    pwrite(stream, "Nbase = $nidx")
    pwrite(stream, "")
    
    if status == 1
        pwrite(stream, "| Solucao otima obtida:")
        pwrite(stream, "| ---------------------")
        pwrite(stream, "| x = $x")
        pwrite(stream, "| z = $z")
        pwrite(stream, "| status = $status")
        pwrite(stream, "")
    elseif status == -1
        pwrite(stream, "| Solucao ilimitada obtida:")
        pwrite(stream, "| -------------------------")
        pwrite(stream, "| de = $x")
        pwrite(stream, "| z = $z")
        pwrite(stream, "| status = $status")
        pwrite(stream, "")
    end
end

function pwrite(stream::IOStream, string::AbstractString)
    println(string)
    write(stream, string * "\n")
end

function problemas()
    cd(pwd())
    # 2)

    # a) Problema da Producao
    println("a) Problema da Producao")
    println("")
    A = float([2 1; 1 2])
    b = float([4 ; 4])
    c = float([4 ; 3])
    x,z,status = interior_points(A, b, c)
    
    # b) Prob 2
    println("b) Problema ilimitado")
    println("")
    A = float([0.5 -1 1 0; -4 1 0 1])
    b = float([0.5 ; 1])
    c = float([1 ; 1; 0; 0])
    x,z,status = interior_points(A, b, c)
    
    # c) Prob 3 - fase 1
    println("c) Problema fase 1")
    println("")
    A = float([2 1 1 0 0; 1 2 0 1 0; -1 -1 0 0 1])
    b = float([4 ; 4 ; -1])
    c = float([4 ; 3; 0; 0; 0])
    x,z,status = interior_points(A, b, c)
end