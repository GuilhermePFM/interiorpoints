#
# Interior Points
#

# Guilherme Pereira Freire Machado

function interior_points(A::Array{Float64}, b::Array{Float64}, c::Array{Float64})

    open_log_i(A, b, c)
    x, p, s, status = interior_bigM(A, b, c)

    return x, p, s, status
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
    pwrite(stream, "Usando big M para comecar no ponto viavel x = e")
    pwrite(stream, "Problema:")
    pwrite(stream, "A = $A_1")
    pwrite(stream, "b = $b_1")
    pwrite(stream, "x0 = $x0")
    pwrite(stream, "s0 = $s0")
    pwrite(stream, "p = $p")
    pwrite(stream, "")

    x, p, s, status = interior_algorithm(A_1, b_1, c_1, x0, s0, p, stream)

    return x, p, s, status
end 

function interior_algorithm(A, b, c, x, s, p, stream)
    err = 1e-5
    α = 0.9
    ρ = 1.0
    n = size(A)[2]
    μ = 0.0
    dx = 0.0

    maxit = 200
    k=0
    for i in 1:maxit
        pwrite(stream, "It: $i - ϵ = $(convergence_error(s, x, μ)) - x = $x")

        # 2) 1st test for convergence
        if i > 40 && optimality_test(s, x, μ, err) && maximum(abs.(dx)) < err
            # convergiu
            status = 1

            pwrite(stream, "Interior Points algorithm converged! ( s*x criteria)")

            if check_unbounded(A, b, c, x)
                status = -1
            end

            # if check_infeasible(x)
            #     status = -2
            # end

            result_log_i(i, x,  (c'*x), status, stream)
            return x, s, p, status
        end
        
        # 3) computation of newton directions
        μ = 0.9 *  x' * s / n 
        dx, dp, ds = compute_directions(x, s, p, μ, A, b, c)

        # 2nd test for convergence
        if maximum(abs.([dx ; dp; ds])) <= err
            pwrite(stream, "Interior Points algorithm converged! ( abs(d) criteria)")
            # convergiu
            status = 1
            println(p)
            println(s)
            if check_unbounded(A, b, c, x)
                status = -1
            end
            # if check_infeasible(x)
            #     status = -2
            # end

            result_log_i(i, x,  (c'*x), status, stream)
            return x, s, p, status
        end

        # 4) and 5) update variables
        x, s, p, unbounded = update_variables(x, s, p, dx, dp, ds, α)

        # update rho in case where convergence get stuck
        # ρ = update_rho(ρ, x_new, x)
        # x = x_new

        if unbounded
            status = -1
            result_log_i(i, x, (c'*x), status, stream)
            return x, s, p, status
        end

        k += 1
    end

    pwrite(stream, "Maximum number of iterations($maxit) exceeded!")
    result_log_i(maxit, x, (c'*x), 0, stream)

    return x, s, p, 0
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

function check_unbounded(A, b, c, x)
    # println("A=$A")
    # println("b=$b")
    # println("c=$c")
    # println("x = $x")
    x=round(x,3)

    xb = x[x.>0]
    B = A[:, x.>0]
    cb = c[x.>0]
    
    xn = x[x.==0]
    N = A[:, x.==0]
    cn = c[ x.==0]

    return all(cb'*B^-1*N .< cn')
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
    ratio_x = round(-x ./ dx, 3)
    ratio_x[ -ratio_x .>= 0] = Inf
    βp = minimum([1, α*minimum(ratio_x)]) 
    
    # compute step β for dual slack
    ratio_s = round(-s ./ ds, 3)
    ratio_s[-ratio_s .>= 0] = Inf
    βd = minimum([1, α*minimum(ratio_s)]) 
    
    # compute step β for dual
    ratio_p = round(-p ./ dp, 3)
    ratio_p[-ratio_p .>= 0] = Inf
    βd = minimum([1, α*minimum(ratio_p)]) 
    
    # update variables
    x = x + βp * dx
    s = s + βd * ds
    p = p + βd * dp
    unbounded = false

    return x, s, p, unbounded
end

function update_rho(ρ, x_new, x_old)
    if maximum(abs.(x_new - x_old)) < 1e-3 && ρ > 0.5
        ρ -= 0.05
    end
    return ρ
end 

function open_log_i(A::Array{Float64,2}, b::Array{Float64,1}, c::Array{Float64,1})
    fname = "Interior_Points.log"
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

function get_log_i()
    fname = "Interior_Points.log"
    stream = open(fname, "a")
    pwrite(stream, "Interior Points")
    pwrite(stream, "--------------")

    return stream
end

function result_log_i(it::Int, x::Array{Float64,1}, z::Float64, status::Int, stream::IOStream)
    pwrite(stream, "iter $it:")
    pwrite(stream, "x = $x")
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
        pwrite(stream, "| x = $x")
        pwrite(stream, "| z = $z")
        pwrite(stream, "| status = $status")
        pwrite(stream, "")
    elseif status == 0 
        pwrite(stream, "| Solucao subotima obtida:")
        pwrite(stream, "| -------------------------")
        pwrite(stream, "| x = $x")
        pwrite(stream, "| z = $z")
        pwrite(stream, "| status = $status")
        pwrite(stream, "")
    elseif status == -2 
        pwrite(stream, "| Problema e inviavel:")
        pwrite(stream, "| -------------------------")
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
    A = float([2 1 1 0; 1 2 0 1])
    b = float([4 ; 4])
    c = float([4 ; 3; 0; 0])
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