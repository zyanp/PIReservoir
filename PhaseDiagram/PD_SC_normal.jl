using Random, Distributions,LinearAlgebra,Statistics, HDF5, Distributed,Printf,Dates
#Random Distributions LinearAlgebra Statistics HDF5 Distributed Printf Dates


N_PROCESS = 68
addprocs(N_PROCESS)
@show nprocs(), workers()



@everywhere module myModule

    using Random, Distributions,LinearAlgebra,Statistics, HDF5, Distributed,Printf,Dates

    function Run_Reservoir(p_number::Int64, N_PROCESS::Int64, n_ens::Int64, nv_ite::Int64, nm_ite::Int64, Vs::Vector{Float64}, Ms::Vector{Float64}, N_res::Int64, Nt::Int32, alpha::Float64, PUT::Float64)

        r1::Vector{Float64} = zeros(N_res)
        r1s::Matrix{Float64} = zeros(n_ens,N_res)
        M::Matrix{Float64} = zeros(N_res, N_res)
        J::Float64 = 0.0
        J0::Float64 = 0.0
        mean_ind::Int64 = 0
        var_ind::Int64 = 0


        for ind in p_number:N_PROCESS:nv_ite*nm_ite-1

            mean_ind = div(ind,nv_ite) + 1
            var_ind = rem(ind,nv_ite) + 1
            J = 1.0 / Vs[var_ind]
            J0 = Ms[mean_ind] * J
            # d = Gamma(J0^2/J^2/N_res, J^2/J0)

            r1s = zeros(n_ens,N_res)


            for ens_ind in 1:n_ens

                # println(@sprintf("%dth ens of (%.2f,%.2f) started",ens_ind,Ms[mean_ind],Vs[var_ind]))

                # M = reshape(rand(d, N_res*N_res), N_res, N_res)
                M = randn(N_res, N_res) .* (J/sqrt(N_res)) .+ (J0/N_res)

                r1 = rand(N_res)

                for ti in 1:Nt
                    r1 = alpha * (M * tanh.(r1)) + (1.0-alpha) * r1
                end

                r1s[ens_ind,:] = r1

            end

            h5open(@sprintf("./Data/rs_normal_%f_(%.2f,%.2f).h5",PUT,Ms[mean_ind],Vs[var_ind]), "w") do file
                write(file, "rs", r1s)
            end

        end

    end

    function main(N_PROCESS::Int64)

        Nt::Int32 = 5000

        N_res::Int64 = 500
        alpha::Float64 = 0.2

        n_ens::Int64=25
        nm_ite::Int64=96
        nv_ite::Int64=96

        # M_max=1.00
        M_min::Float64=0.10
        # V_max=2.00
        V_min::Float64=0.10
        dV::Float64=0.02
        dM::Float64=0.02

        Vs::Vector{Float64} = range(V_min,step=dV,length=nv_ite)
        Ms::Vector{Float64} = range(M_min,step=dM,length=nm_ite)

        PUT::Float64 = datetime2unix(Dates.now())
        pmap(p_number -> Run_Reservoir(p_number, N_PROCESS, n_ens, nv_ite, nm_ite, Vs, Ms, N_res, Nt, alpha, PUT), 0:N_PROCESS-1)

    end
end

using .myModule
@time myModule.main(N_PROCESS)



# @time main()
