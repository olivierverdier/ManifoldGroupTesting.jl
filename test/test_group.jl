

grp_rep(::Any, x) = x
grp_rep(G, ::Identity) = grp_rep(G, identity_element(G))
alg_rep(::Any, x) = x

@testset "Group Testing.jl" for G in [SpecialOrthogonal(3)]
    χ1, χ2 = [rand(rng, G) for i in 1:2]
    ξ1, ξ2 = [rand_lie(rng, G) for i in 1:2]
    @test GT.check_exp_lie_point(G, ξ1)
    @test GT.check_adjoint_action_in_alg(G, χ1, ξ1)
    @test GT.check_grp_rep_Identity(G, grp_rep)
    @test GT.check_grp_rep_compose(G, grp_rep, χ1, χ2)
    @test GT.check_alg_rep(G, alg_rep, ξ1, ξ2)
    @test GT.check_zero_Identity(G) broken=true # fails for SpecialOrthogonal
    @test GT.check_exp_ad(G, ξ1, ξ2)
    @test GT.check_adjoint_action(G, grp_rep, alg_rep, χ1, ξ1)
    @test GT.check_inv_rep(G, grp_rep, χ1)
    @testset "$side" for side in [LeftSide(), RightSide()]
        @test GT.check_apply_diff_group_at_id(G, ξ1, side, Identity)
        @test GT.check_apply_diff_group_at_id(G, ξ1, side, identity_element)
        @test GT.check_inv_diff(G, χ1, ξ1, side)
    end
end
