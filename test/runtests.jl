import ManifoldGroupTesting as GT
using Test
using Manifolds
import ManifoldGroupUtils: rand_lie, translate_from_id
import Random
rng = Random.default_rng()


include("test_group.jl")
include("test_action.jl")
