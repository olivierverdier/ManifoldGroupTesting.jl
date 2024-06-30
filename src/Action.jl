
_compose_dir(G, χ1, χ2, ::LeftAction) = compose(G, χ1, χ2)
_compose_dir(G, χ1, χ2, ::RightAction) = compose(G, χ2, χ1)

"""
    check_action_morphism(α::GroupAction, χ1, χ2, p)

For a group action `α`,
```math
α(χ_1, α(χ_2, p)) = α(m(χ_1, χ_2), p)
```
where for a *left* action
```math
m(χ_1, χ_2) = χ_1 χ_2
```
and for a *right* action
```math
m(χ_1, χ_2) = χ_2 χ_1
```
"""
check_action_morphism(A::AbstractGroupAction{TAD}, χ1, χ2, p) where {TAD} = begin
    G = base_group(A)
    M = group_manifold(A)
    prod_first = apply(A, _compose_dir(G, χ1, χ2, TAD()), p)
    act_first = apply(A, χ1, apply(A, χ2, p))
    return isapprox(M, prod_first, act_first)
end



"""
    check_apply_morphism_Identity(α::AbstractGroupAction, p)

For any action ``α`` and ``p`` in the manifold, one has
```math
α(1, p) = p
```
where ``1`` denotes ``Identity(G)`` where ``G`` is the
action group.
"""
check_apply_morphism_Identity(A::AbstractGroupAction, p) = begin
    G = base_group(A)
    p_ = apply(A, Identity(G), p)
    return isapprox(group_manifold(A), p, p_)
end

"""
    check_trivial_infinitesimal_action(A::GroupAction, p)

For an action ``α`` and a point ``p∈ M``,
consider the function ``f : G → M`` defined as
``f(χ) := α(χ, p)``.
The tangent map at identity ``D := T f(1)`` is a linear map
and sends zero to zero:
```math
D 0 = 0
```
"""
check_trivial_infinitesimal_action(A::AbstractGroupAction, p, id=Identity) = begin
    G = base_group(A)
    ξ = zero_vector(G, id(G))
    computed = apply_diff_group(A, Identity(G), ξ, p)
    M = group_manifold(A)
    expected = zero_vector(M, p)
    TM = TangentSpace(group_manifold(A), p)
    return isapprox(TM, computed, expected)
end


"""
    check_switch_action_direction(α::GroupAction, χ, p)

For any group action ``α`` and its dual ``α'`` obtained
by `switch_direction` one has
```math
α(χ,p) = α'(χ^{-1}, p)
```
"""
check_switch_action_direction(A::AbstractGroupAction, χ, p) = isapprox(group_manifold(A), apply(switch_direction(A), χ, p), apply(A, inv(base_group(A), χ), p))
