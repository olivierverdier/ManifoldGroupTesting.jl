
#--------------------------------


_get_side_from_action_dir(::LeftAction) = RightSide()
_get_side_from_action_dir(::RightAction) = LeftSide()

_transporter(G, χ, ξ, dir) = translate_from_id(G, χ, ξ, _get_side_from_action_dir(dir))

"""
    check_apply_diff_group(
        A::AbstractGroupAction{TAD}, # group action G ⊂ Diff(M)
        χ, # element of G
        ξ, # element of Alg(G)
        p, # element of M
        id_func, # `identity_element` or `Identity`
    )

This should hold for *any* group action ``α`` on any manifold.
If you define ``π_p(χ) := α(χ, p)`` for ``χ ∈ G`` and ``p ∈ M``,
and define, for ``ξ ∈ Alg(G)``,
 ``T_R(χ, ξ) := ξχ`` (the right translation),
and ``T_L(χ, ξ) := χξ`` (the left translation), then we have the identity:
```math
⟨Dπ_{p}(χ), T(χ, ξ)⟩ = ⟨Dπ_{α(χ,p)}(1), ξ⟩
```
where, for a *left* action, ``T`` is the *right* translation,
and for a *right* action, ``T`` is the *left* translation.
"""
check_apply_diff_group(A::AbstractGroupAction{TAD}, χ, ξ, p, id_func=Identity) where {TAD} = begin
    G = base_group(A)
    p_ = apply(A, χ, p)
    v1 = apply_diff_group(A, χ, _transporter(G, χ, ξ, TAD()), p)
    v2 = apply_diff_group(A, id_func(G), ξ, p_)
    return isapprox(TangentSpace(G, p_), v1, v2)
end

#--------------------------------

"""
    check_inv_diff(
      G, # Group
      χ, # group element
      ξ, # Lie algebra element
      side::Manifolds.GroupActionSide,
      )

Test the differential of the inverse on a Lie group `G`.
Denote this inverse by ``I(χ) := χ^{-1}``.
If the left and right transports are ``T_L(χ,ξ) := χξ``
and ``T_R(χ,ξ) := ξχ`` respectively, then
```math
⟨DI(χ), T_L(χ,ξ)⟩ = -T_R(χ^{-1}, ξ)
```
and
``` math
⟨DI(χ), T_R(χ,ξ)⟩ = -T_L(χ^{-1}, ξ)
```
"""
check_inv_diff(G, χ, ξ, side) = begin
    χ_ = inv(G, χ)
    computed = inv_diff(G, χ, translate_from_id(G, χ, ξ, side))
    expected = -translate_from_id(G, χ_, ξ, switch_side(side))
    return isapprox(TangentSpace(G, χ_), computed, expected)
end
