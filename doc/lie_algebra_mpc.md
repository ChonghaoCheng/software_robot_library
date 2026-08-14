# SerialLinkLieAlgebraMPC

`SerialLinkLieAlgebraMPC` is a time-indexed Cartesian velocity controller. It
is intended as a controlled comparison between `SerialLinkMPC` and
`SerialLinkRMPCC`:

| controller | pose state | reference index | optimized progress |
|---|---|---|---|
| `SerialLinkMPC` | absolute `[p; rotation_vector]` in base | time | no |
| `SerialLinkLieAlgebraMPC` | `Log(T_ref^-1 T)` in reference/body coordinates | time | no |
| `SerialLinkRMPCC` | `Log(T_ref(s)^-1 T)` | path progress `s` | yes |

Thus a Cartesian-MPC versus Lie-algebra-MPC comparison isolates the pose-error
geometry without also changing the reference from time tracking to contouring
control.

## Error dynamics

Let

```text
T_dot     = T     hat(xi)
T_ref_dot = T_ref hat(xi_ref)
E         = inverse(T_ref) T   # left-invariant under a shared world transform
e         = Log(E)
```

where `xi` and `xi_ref` are body twists ordered `[v; omega]`. The exact group
error dynamics are

```text
E_dot = E hat(xi) - hat(xi_ref) E.
```

Writing `xi = xi_ref + delta_xi` and retaining first-order terms around
`E = I` gives

```text
e_dot = -ad(xi_ref) e + delta_xi.
```

Forward-Euler discretization produces the time-varying linear prediction

```text
e[k+1] = (I - dt ad(xi_ref[k])) e[k] + dt delta_xi[k].
```

For the library's `[v; omega]` ordering,

```text
ad([v; omega]) = [ skew(omega)  skew(v)     ]
                 [ 0            skew(omega) ].
```

The QP minimizes the stacked pose error and twist correction with the same
default horizon, sampling time, position/orientation weights, control weights,
and Cartesian velocity limits as `SerialLinkMPC`. Its box bounds apply to the
actual body twist `xi_ref + delta_xi`.

The selected body twist is rotated into base axes at the endpoint origin and
then passed to `SerialLinkVelocityBase::resolve_endpoint_twist()`. This rotation
does not use the full transform adjoint: adding `p x omega` would change the
linear velocity point and would be inconsistent with the manipulator geometric
Jacobian.

## Reference implementation boundary

This is the convex, first-order Lie-algebra error-state MPC comparator. It is
not a nonlinear manifold MPC, does not optimize joint dynamics, and does not
claim the configuration-independent stability result of a particular rigid
body dynamics paper. The library test checks first-order consistency against
direct SE(3) multiplication; closed-loop experiments determine whether the
approximation is useful for the UR5e endpoint-velocity task.

Primary methodological references:

- Sangli Teng and Maani Ghaffari, *A Lie Algebraic Model Predictive Control for
  Legged Robot Control: Implementation and Stability Analysis*.
- Jiawei Tang et al., *GMPC: Geometric Model Predictive Control for Wheeled
  Mobile Robot Trajectory Tracking*, 2024.
