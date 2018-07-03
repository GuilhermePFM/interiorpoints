#
# Interior Points
#

# Guilherme Pereira Freire Machado

function interior_points(A::Array{Float64}, b::Array{Float64}, c::Array{Float64}, debug=true)

    open_log_i(A, b, c)
    stream = get_log_i()

    # initial guess
    m, n = size(A)

    x0 = ones(n)
    p0 = zeros(m)
    s0 = ones(n)
    
    x = ones(n)
    p = zeros(m)
    s = ones(n)

    x, s, p, status, it = interior_algorithm(A, b, c, x0, s0, p0, stream, debug)
    # x, p, s, status, it = interior_bigM(A, b, c, debug)

    return x, p, s, status, it
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

function interior_bigM(A, b, c, debug)
    stream = get_log_i()
    
    m, n = size(A)
    
    U = maximum(c)
    M = U*1e5
    
    # add bigM variable if necessary
    if any((b - A*ones(n)) .!= 0)
        c_1 = [c ; M]
        A_1 = [A (b - A*ones(n))]
        b_1 = b
        # initial feasible solution
        x0 = ones(n+1)
    else
        c_1 = c
        A_1 = A
        b_1 = b
        # initial feasible solution
        x0 = ones(n)
    end
    p = (c'[1:(m)]*A[:,1:(m)])' # ones(m)#
    s0 = (c_1' - p'A_1)'
    
    A=A_1
    c=c_1
    x=x0
    s=s0
    pwrite(stream, "Usando big M para comecar no ponto viavel x = e", debug)
    pwrite(stream, "Problema:", debug)
    pwrite(stream, "A = $A_1", debug)
    pwrite(stream, "b = $b_1", debug)
    pwrite(stream, "x0 = $x0", debug)
    pwrite(stream, "s0 = $s0", debug)
    pwrite(stream, "p = $p", debug)
    pwrite(stream, "")

    x, s, p, status, it = interior_algorithm(A_1, b_1, c_1, x0, s0, p, stream, debug)

    return x, p, s, status, it
end 

function interior_algorithm(A, b, c, x, s, p, stream, debug=true)
    err = 1e-3
    α = 0.9
    ρ = 0.5
    n = size(A)[2]
    μ = 0.0
    dx = 0.0

    maxit = 200
    it=0
    for i in 1:maxit
        pwrite(stream, "It: $i - ϵ = $(convergence_error(s, x, μ)) - x = $x", debug)

        # 2) 1st test for convergence
        if s'*x < err  # i > 40 && optimality_test(s, x, μ, err) && maximum(abs.(dx)) < err
            # convergiu
            status = 1

            if norm(x) < norm(p) # check for unbounded problem
                status = -1
                pwrite(stream, "The problem is unbounded!")
            else

                pwrite(stream, "Interior Points algorithm converged! ( s*x criteria)")
            end


            # if check_unbounded(A, b, c, x)
            #     status = -1
            # end

            # if check_infeasible(x)
            #     status = -2
            # end

            result_log_i(i, x,  (c'*x), status, stream)
            return x, s, p, status, it
        end
        
        # 3) computation of newton directions
        μ = ρ *  x' * s / n 
        dx, ds, dp = compute_directions_old(x, s, p, μ, A, b, c)
        ρ = update_rho(ρ, x + dx, x)

        # 4) and 5) update variables
        x, s, p, unbounded = update_variables(x, s, p, dx, ds, dp, α)

        if  convergence_error(s, x, μ) > 1e5 &&  all(x .> 0)
            status = -1
            result_log_i(i, x, (c'*x), status, stream, debug)
            return x, s, p, status, it
        end

        it += 1
    end

    pwrite(stream, "Maximum number of iterations($maxit) exceeded!")
    result_log_i(maxit, x, (c'*x), 0, stream)

    return x, s, p, 0, it
end

function optimality_test(s::Array{Float64}, x::Array{Float64}, μ::Float64, err::Float64)
    nx = length(x)
    op = convergence_error(s, x, μ)
    return op < err
end

function check_infeasible(x)
    if round(x[end],2) > 0
        println(x)
        return true
    end

    return false
end

function convergence_error(s::Array{Float64}, x::Array{Float64}, μ::Float64)
    nx = length(x)
    op = x'*s - (μ*ones(nx)')*ones(nx)
    return op
end

function compute_directions(x::Array{Float64}, s::Array{Float64}, p::Array{Float64}, μ::Float64,  A::Array{Float64}, b::Array{Float64}, c::Array{Float64})
    # solve linear system:
    nx = length(x)
    np = length(p)
    ns = length(s)
    e = ones(nx)

    inv_s = diagm(s) \ I 
    D2 = diagm(x) * inv_s
    D = D2^0.5

    parcial_1 = (A * D2 * A') \ (A * D)
    P = D * A' * parcial_1
    vμ = (diagm(x) \ D) * (μ * ones(nx) - diagm(x) * diagm(s) * ones(nx))
    

    dx = D * (I - P) * vμ
    dp = -((A * D2 * A') \ (A*D)) * vμ
    ds = (D \ P) * vμ

     return dx, ds, dp
end

function compute_directions_old(x::Array{Float64}, s::Array{Float64}, p::Array{Float64}, μ::Float64,  A::Array{Float64}, b::Array{Float64}, c::Array{Float64})
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

     return dx, ds, dp
end

function update_variables(x::Array{Float64}, s::Array{Float64}, p::Array{Float64}, dx::Array{Float64}, ds::Array{Float64}, dp::Array{Float64}, α::Float64)
    # compute step β for primal
    ratio_x = -x ./ dx
    ratio_x[ dx .>= 0] = Inf
    # ratio_x = [v for (i,v) in enumerate(ratio_x) if dx[i] < 0]
    βp = minimum([1, α*minimum(ratio_x)]) 
    
    # compute step β for dual slack
    ratio_s = -s ./ ds
    ratio_s[ds .>= 0] = Inf
    # ratio_s = [v for (i,v) in enumerate(ratio_s) if ds[i] < 0]
    βd = minimum([1, α*minimum(ratio_s)]) 
    
    # update variables
    x = x + βp * dx
    s = s + βd * ds
    p = p + βd * dp

    return x, s, p, false
end

function update_rho(ρ, x_new, x_old)
    if maximum(abs.(x_new - x_old)) < 1
        ρ *= 0.9

    else
        ρ = 0.8
    end
    return ρ
end 

function open_log_i(A::Array{Float64,2}, b::Array{Float64,1}, c::Array{Float64,1}, debug=true)
    fname = "Interior_Points.log"
    if isfile(fname)
        stream = open(fname, "a")
        pwrite(stream, "=======================", debug)
        pwrite(stream, "Comeco da Solucao do PL", debug)
        pwrite(stream, "=======================", debug)
        pwrite(stream, "Problema:", debug)
        pwrite(stream, "A = $A", debug)
        pwrite(stream, "b = $b", debug)
        pwrite(stream, "c = $c", debug)
        pwrite(stream, "", debug)
        close(stream)
    else
        stream = open(fname, "w", debug)
        pwrite(stream, "=======================", debug)
        pwrite(stream, "Comeco da Solucao do PL", debug)
        pwrite(stream, "=======================", debug)
        pwrite(stream, "Problema:", debug)
        pwrite(stream, "A = $A", debug)
        pwrite(stream, "b = $b", debug)
        pwrite(stream, "c = $c", debug)
        pwrite(stream, "", debug)
        close(stream)
    end
    nothing
end

function get_log_i()
    fname = "Interior_Points.log"
    stream = open(fname, "a")
    pwrite(stream, "Interior Points")
    pwrite(stream, "--------------")

    return stream
end

function result_log_i(it::Int, x::Array{Float64,1}, z::Float64, status::Int, stream::IOStream, debug=true)
    pwrite(stream, "iter $it:", debug)
    pwrite(stream, "x = $x", debug)
    pwrite(stream, "", debug)
    
    if status == 1
        pwrite(stream, "| Solucao otima obtida:", debug)
        pwrite(stream, "| ---------------------", debug)
        pwrite(stream, "| x = $x", debug)
        pwrite(stream, "| z = $z", debug)
        pwrite(stream, "| status = $status", debug)
        pwrite(stream, "")
    elseif status == -1
        pwrite(stream, "| Solucao ilimitada obtida:", debug)
        pwrite(stream, "| -------------------------", debug)
        pwrite(stream, "| x = $x", debug)
        pwrite(stream, "| z = $z", debug)
        pwrite(stream, "| status = $status", debug)
        pwrite(stream, "", debug)
    elseif status == 0 
        pwrite(stream, "| Solucao subotima obtida:", debug)
        pwrite(stream, "| -------------------------", debug)
        pwrite(stream, "| x = $x", debug)
        pwrite(stream, "| z = $z", debug)
        pwrite(stream, "| status = $status", debug)
        pwrite(stream, "")
    elseif status == -2 
        pwrite(stream, "| Problema e inviavel:", debug)
        pwrite(stream, "| -------------------------", debug)
        pwrite(stream, "| status = $status", debug)
        pwrite(stream, "", debug)
    end
end

function pwrite(stream::IOStream, string::AbstractString, debug=true)
    if debug
        println(string)
    end
    write(stream, string * "\n")
end

function problemas()
    cd(pwd())
    # 2)

    # a) Problema da Producao
    println("a) Problema da Producao")
    println("")
    A = float([2 1 1 0; 1 2 0 1])
    b = float([4 ; 4])
    c = -float([4 ; 3; 0; 0])
    x, p, s, status, it = interior_points(A, b, c)
    
    # b) Prob 2
    println("b) Problema ilimitado")
    println("")
    A = float([0.5 -1 1 0; -4 1 0 1])
    b = float([0.5 ; 1])
    c = -float([1 ; 1; 0; 0])
    x, p, s, status, it = interior_points(A, b, c)
     
    # c) Prob 3 - fase 1
    println("c) Problema fase 1")
    println("")
    A = float([2 1 1 0 0; 1 2 0 1 0; -1 -1 0 0 1])
    b = float([4 ; 4 ; -1])
    c = -float([4 ; 3; 0; 0; 0])
    x, p, s, status, it = interior_points(A, b, c)

    # inviavel
    A = float([1 3 1 0 0 ;3 2 0 1 0;-1 -3 0 0 1]); 
    b = float([8 ;12 ;-13]);
    c = -float([ 1; 1 ; 0; 0;0]);
    x, p, s, status, it = interior_points(A, b, c)

end