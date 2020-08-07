# get hold of all the predefined math constants in Julia
import Base.MathConstants: π, pi, ℯ, e, γ, eulergamma, catalan, φ, golden
Base.@irrational twoπ  6.28318530717958647692528676655900576839433879  big(2.0)*pi
