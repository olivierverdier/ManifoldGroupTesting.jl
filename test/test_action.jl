
@testset "Action $name" for (name, A) in [
    "Rot plane" => RotationAction(Euclidean(3), SpecialOrthogonal(3)),
    "SO L" => GroupOperationAction(SpecialOrthogonal(3), Manifolds.LeftForwardAction()),
    "SO R*" => GroupOperationAction(SpecialOrthogonal(3), Manifolds.RightForwardAction()),
    "SO L*" => GroupOperationAction(SpecialOrthogonal(3), Manifolds.LeftBackwardAction()),
]
    G = base_group(A)
    χ1, χ2 = [rand(rng, G) for i in 1:2]
    ξ1, ξ2 = [rand_lie(rng, G) for i in 1:2]
    p = rand(rng, group_manifold(A))
    @test GT.check_action_morphism(A, χ1, χ2, p)
    @test GT.check_apply_morphism_Identity(A, p)
    @test GT.check_trivial_infinitesimal_action(A, p, identity_element)
    @test GT.check_switch_action_direction(A, χ1, p)
    @test GT.check_apply_diff_group(A, χ1, ξ1, p) broken = A isa RotationAction
end
