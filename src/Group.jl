using ManifoldGroupUtils
import ManifoldGroupUtils: matrix_from_lin_endomorphism
"""
A collection of general purpose methods
for testing Lie groups.
"""

"""
``Ad`` and ``exp`` commute:
```math
Ad_{exp(ξ_1)}ξ_2 = exp(ad_{ξ_1}) ξ_2
```
"""
function check_exp_ad(G, vel, tvel)
    χ = exp_lie(G, vel)
    Ad_exp = adjoint_action(G, χ, tvel)

    lie_bracket(G, vel, tvel)
    B = DefaultOrthogonalBasis()
    der = matrix_from_lin_endomorphism(G, ξ -> lie_bracket(G, vel, ξ), B)
    mor = exp(der)
    tvel_coord = get_coordinates_lie(G, tvel, B)
    exp_ad = get_vector_lie(G, mor * tvel_coord, B)
    return isapprox(algebra(G), Ad_exp, exp_ad)
end

"""
    check_adjoint_action(G, grp_rep, alg_rep, χ, ξ)

The group representation ``ρ`` and algebra representation ``ρ``
 commute with the adjoint action:
```math
ρ(χ) ρ(ξ) ρ(χ)^{-1} = ρ(Ad_χ (ξ))
```
"""
check_adjoint_action(G, grp_rep, alg_rep, χ, ξ) = begin
    mat = grp_rep(G, χ)
    matinv = inv(mat)
    expected = mat * alg_rep(G, ξ) * matinv
    computed = alg_rep(G, adjoint_action(G, χ, ξ))
    return isapprox(computed, expected)
end


"""
    check_inv_rep(G, grp_rep, χ)

The group representation ``ρ`` commutes with the inverse.
```math
ρ(χ)^{-1} = ρ(χ^{-1})
```
"""
check_inv_rep(G, grp_rep, χ) = begin
    computed = grp_rep(G, inv(G, χ))
    expected = inv(grp_rep(G, χ))
    return isapprox(computed, expected)
end



_switch_sign(ξ, ::LeftSide) = ξ
_switch_sign(ξ, ::RightSide) = -ξ
"""
    check_apply_diff_group_at_id(G, side::GroupActionSide)

The left group operation action on itself ``α(χ_1)χ_2``
is either (left side)
```math
α(χ_1)χ_2 = χ_1 χ_2
```
or (right side)
```math
α(χ_1)χ_2 = χ_2 χ_1^{-1}
```
Now fix ``χ_2 = 1`` (1 is the identity of ``G``) and define ``f : G → G`` by ``f(χ) := α(χ) 1``. Since ``f(1) = 1``,
its differential at identity is a map ``Alg(G) → Alg(G)``.
This map is either
- ``Id`` (left side)
- ``-Id`` (right side)
"""
check_apply_diff_group_at_id(G, ξ, side::Manifolds.GroupActionSide, id_func=Identity) = begin
    id = id_func(G)
    ξ_ = apply_diff_group(GroupOperationAction(G, (LeftAction(), side)), id, ξ, id)
    return isapprox(algebra(G), ξ, _switch_sign(ξ_, side))
end

"""
    check_exp_lie_point(G, ξ)

The Lie group exponential sends the vector ξ
to an element in the group.
"""
check_exp_lie_point(G, ξ) = is_point(G, exp_lie(G, ξ))

"""
    check_adjoint_action_in_alg(G, χ, ξ)

The adjoint action of χ on ξ is an element of Alg(G):
```math
Ad_{χ}ξ ∈ Alg(G)
```
"""
check_adjoint_action_in_alg(G, χ, ξ) = is_vector(G, identity_element(G), adjoint_action(G, χ, ξ))

"""
    check_grp_rep_Identity(G)

The representation works at the Identity(G) point.
"""
check_grp_rep_Identity(G, grp_rep) = begin
    expected = grp_rep(G, identity_element(G))
    computed = grp_rep(G, Identity(G))
    return isapprox(expected, computed)
end

"""
    check_grp_rep_compose(G, ρ, χ1, χ2)

The group representation ``ρ`` commutes with composition.
```math
 ρ(χ_1 χ_2) = ρ(χ_1) ρ(χ_2)
```
where the latter is a matrix multiplication.
"""
check_grp_rep_compose(G, grp_rep, χ1, χ2) = begin
    m1, m2 = [grp_rep(G, p) for p in [χ1, χ2]]
    expected = m1 * m2
    computed = grp_rep(G, compose(G, χ1, χ2))
    return isapprox(expected, computed)
end

"""
    check_alg_rep(G, alg_rep, ξ1, ξ2)

The algebra representation ``ρ`` is an algebra morphism.
```math
ρ([ξ_1, ξ_2]) = [ρ(ξ_1), ρ(ξ_2)]
```
where the latter is a matrix commutator.
"""
check_alg_rep(G, alg_rep, ξ1, ξ2) = begin
    m1, m2 = [alg_rep(G, v) for v in [ξ1,ξ2]]
    expected = m1*m2 - m2*m1
    computed = alg_rep(G, lie_bracket(G, ξ1, ξ2))
    return isapprox(expected, computed)
end

check_zero_Identity(G) = isapprox(algebra(G),
                                  zero_vector(G, Identity(G)),
                                  zero_vector(G, identity_element(G)))

"""
    check_exp_invariant(G, exp, χ, v, χ_)

The invariant exponential of  a Lie group fulfils
```math
χ' \exp_{χ}(v) = \exp_{χ'χ}(χ' v)
```

There is a right version which is not implemented. It would check that
```math
\exp_χ(v) χ' = \exp_{χχ'}(v χ')
```
"""
check_exp_invariant(G, exp, χ, v, χ_) = begin
    χ1 = translate(G, χ_, exp(G, χ, v), (LeftAction(), LeftSide()))
    v_ = translate_diff(G, χ_, χ, v, (LeftAction(), LeftSide()))
    χ_χ = translate(G, χ_, χ, (LeftAction(), LeftSide()))
    χ2 = exp(G, χ_χ, v_)
    return isapprox(G, χ1, χ2)
end

"""
    check_exp_exp_(G, exp, exp!, χ_, χ, v)

Compare `exp` and `exp!`.
"""
check_exp_exp_(G, exp, exp!, χ_, χ, v) = begin
    χ1 = exp!(G, χ_, χ, v)
    χ2 = exp(G, χ, v)
    return χ1 === χ_ && isapprox(G, χ1, χ2)
end
"""
    check_log_log_(G, log, log!, v_, χ1, χ2)

Compare `log` and `log!`.
"""
check_log_log_(G, log, log!, v_, χ1, χ2) = begin
    v__ = log!(G, v_, χ1, χ2)
    v = log(G, χ1, χ2)
    return v__ === v_ && isapprox(TangentSpace(G, χ1), v, v__)
end

"""
    check_exp_log(G, exp, log, χ1, χ2)

Check the identity
```math
exp_{χ_1}(log_{χ_1}(χ_2)) = χ_2
```
"""
check_exp_log(G, exp, log, χ1, χ2) = begin
    v = log(G, χ1, χ2)
    χ_ = exp(G, χ1, v)
    return isapprox(G, χ2, χ_)
end
"""
    check_log_exp(G, log, exp, χ, v)

Check the identity
```math
log_{χ}(exp_{χ}(v)) = v
```
"""
check_log_exp(G, log, exp, χ, v) = begin
    χ_ = exp(G, χ, v)
    v_ = log(G, χ, χ_)
    return isapprox(TangentSpace(G, χ), v, v_)
end
