# Numerical Calculus Library - Implementation Roadmap

## Executive Summary

Transform `algebraic_integrators` into a comprehensive **hybrid symbolic-numeric calculus library** that:
1. Preserves analytical relationships through generic programming
2. Provides numerical methods as fallback
3. Supports multivariate integration (cubature + Monte Carlo)
4. Maintains composable architecture

## What We Have (v1.0.0)

✅ **Univariate Integration** - Multiple methods, adaptive, parallel
✅ **Numerical Differentiation** - 8th-order finite difference, gradients
✅ **ODE Solvers** - Euler, RK4 for 1st/2nd order
✅ **Antiderivatives** - With analytical derivative preservation
✅ **Accumulators** - 5 precision strategies
✅ **118+ Unit Tests** - Comprehensive coverage
✅ **Documentation** - MkDocs with API reference

## What We're Missing

❌ **Multivariate Integration** - 2D/3D/ND cubature and Monte Carlo
❌ **Analytical Preservation** - Generalized beyond antiderivatives
❌ **Root Finding** - If "solvers" means this
❌ **Clear Library Identity** - Name and scope clarification
⚠️ **Legacy Code** - Line integrals need integration or removal

## Phase 1: Library Identity & Cleanup (Week 1)

### 1.1 Define Scope
- [ ] Update README with clear mission statement
- [ ] Decide on name: `numerical_calculus` vs keep `algebraic_integrators`
- [ ] Document what "algebraic operations" means (composability)
- [ ] Clarify "solvers" scope (ODEs only? or add root finding?)

### 1.2 Cleanup Legacy Code
- [ ] **Option A:** Integrate line integrals into composable architecture
  - Move from `numerical_integration` namespace to `algebraic_integrators`
  - Add tests
  - Document
- [ ] **Option B:** Remove line integrals (save for future multivariate work)
  - Delete `include/numerical_integrators/line_integral.hpp`
  - Remove from any references

### 1.3 Update Documentation
- [ ] Rewrite README with:
  - Clear purpose statement
  - Hybrid symbolic-numeric vision
  - Quick start examples
  - Architecture overview
- [ ] Update CLAUDE.md with multivariate plans
- [ ] Link to VISION.md and MULTIVARIATE_DESIGN.md

**Deliverable:** Clear, focused library with no ambiguous features

## Phase 2: Multivariate Integration - Cubature (Weeks 2-3)

### 2.1 Core Infrastructure
- [ ] Define multivariate concepts
  ```cpp
  template <typename F, typename T, int Dim>
  concept MultivariateFunction;

  template <typename R, int Dim>
  concept IntegrationRegion;
  ```
- [ ] Implement region types:
  - [ ] `hyperrectangle<Dim>`
  - [ ] `simplex<Dim>` (triangles, tetrahedra)
  - [ ] `ball<Dim>`
- [ ] Create `integration_result_nd<T, Dim>` with metadata

### 2.2 Tensor Product Integrator
- [ ] 2D tensor product using existing 1D Gauss-Kronrod
- [ ] 3D tensor product
- [ ] Transform support for non-rectangular domains
- [ ] Tests:
  - [ ] Polynomial integrals (exact)
  - [ ] Gaussian functions
  - [ ] Separable functions

### 2.3 Adaptive Cubature
- [ ] Implement Genz-Malik embedded rule
- [ ] Priority queue for adaptive subdivision
- [ ] Error estimation and stopping criteria
- [ ] Tests:
  - [ ] Peak functions
  - [ ] Corner singularities
  - [ ] Oscillatory integrands

**Deliverable:** Production-ready 2D/3D cubature integration

## Phase 3: Multivariate Integration - Monte Carlo (Week 4)

### 3.1 Random Number Generation
- [ ] Implement Sobol sequence generator
- [ ] Implement Halton sequence generator
- [ ] Wrapper for std::mt19937_64
- [ ] Thread-safe RNG for parallel execution

### 3.2 Basic Monte Carlo
- [ ] Simple Monte Carlo with variance estimation
- [ ] Stratified sampling
- [ ] Works with arbitrary regions (rejection sampling)
- [ ] Tests:
  - [ ] High-dimensional Gaussian
  - [ ] Sphere/ball volumes
  - [ ] Discontinuous functions

### 3.3 Quasi-Monte Carlo
- [ ] QMC with Sobol sequences
- [ ] Comparison with regular MC
- [ ] Adaptive sampling based on variance
- [ ] Tests against Genz benchmark suite

### 3.4 Unified Interface
- [ ] Automatic method selection based on dimension
- [ ] `integrate_nd<Dim>` template interface
- [ ] Convenience functions: `integrate_2d`, `integrate_3d`

**Deliverable:** Full-featured multivariate integration up to ~10 dimensions

## Phase 4: Analytical Preservation (Week 5)

### 4.1 Extend Current deriv(antideriv) Pattern
- [ ] Create `derivative_of<F>` wrapper (analogue of `antiderivative_of`)
- [ ] Implement `integrate(derivative_of<F>)` → returns F
- [ ] Tests for round-tripping: `deriv(integrate(f))` and `integrate(deriv(f))`

### 4.2 Concept-Based Dispatch
- [ ] Define `has_analytical_derivative<T>` concept
- [ ] Define `has_analytical_integral<T>` concept
- [ ] Generic `deriv()` with compile-time dispatch
- [ ] Generic `integrate()` with compile-time dispatch

### 4.3 Known Function Library (Optional)
- [ ] `sin_function`, `cos_function`, `exp_function`, etc.
- [ ] Each knows its analytical derivative/integral
- [ ] Automatic composition
- [ ] Example:
  ```cpp
  auto f = sin_function{};
  auto df = deriv(f);  // Returns cos_function{} analytically!
  ```

**Deliverable:** Framework for preserving analytical relationships

## Phase 5: Advanced Features (Weeks 6-7)

### 5.1 Sparse Grids (If time permits)
- [ ] Smolyak algorithm implementation
- [ ] Clenshaw-Curtis nested rules
- [ ] Good for 4-10 dimensions
- [ ] Tests on smooth high-dimensional functions

### 5.2 Root Finding (If in scope)
- [ ] Bisection method
- [ ] Newton-Raphson (uses `deriv()`)
- [ ] Secant method
- [ ] Brent's method
- [ ] Tests on polynomials and transcendental equations

### 5.3 Enhanced ODE Solvers (If in scope)
- [ ] Adaptive step size (RK45, Dormand-Prince)
- [ ] Stiff ODE solvers (implicit methods)
- [ ] Event detection
- [ ] Tests on standard ODE test suite

**Deliverable:** Complete numerical calculus toolkit

## Phase 6: Polish & Release (Week 8)

### 6.1 Performance Optimization
- [ ] Benchmark all integrators
- [ ] Profile hot paths
- [ ] SIMD optimization where applicable
- [ ] Parallel execution for multivariate

### 6.2 Documentation
- [ ] Update all docs with multivariate examples
- [ ] Add tutorial for hybrid symbolic-numeric usage
- [ ] Create Jupyter notebook examples
- [ ] API reference completeness

### 6.3 Testing
- [ ] Achieve >90% code coverage
- [ ] Add fuzzing tests
- [ ] Performance regression tests
- [ ] Integration tests with real-world problems

### 6.4 Release
- [ ] Version 2.0.0 (breaking changes to API)
- [ ] Migration guide from 1.0
- [ ] Announcement blog post
- [ ] Submit to Boost for review (optional)

**Deliverable:** Production-ready v2.0.0 release

## Success Metrics

| Metric | Target | Current | Status |
|--------|--------|---------|--------|
| Test Coverage | >90% | ~85% | 🟡 |
| Dimensions Supported | 1-10 | 1 | 🔴 |
| Integration Methods | 5+ | 3 (1D) | 🟡 |
| Analytical Preservation | Derivatives + Integrals | Derivatives only | 🟡 |
| Documentation Pages | 15+ | 9 | 🟡 |
| Examples | 20+ | ~10 | 🟡 |
| Performance vs GSL | Within 2x | N/A | ⚪ |

## Decision Points

### 1. Library Name
**Options:**
- Keep `algebraic_integrators` (current)
- Rename to `numerical_calculus` (clearer scope)
- Rename to `hybrid_calculus` (emphasizes symbolic-numeric)

**Recommendation:** Rename to `numerical_calculus` for clarity

### 2. Root Finding Scope
**Options:**
- Add root finding (makes "solvers" clear)
- Keep ODE-only (narrower scope)

**Recommendation:** Add basic root finding (bisection, Newton, Brent)

### 3. Symbolic Capability Level
**Options:**
- Minimal (just preserve operations) ← RECOMMENDED for now
- Medium (expression templates, lazy eval)
- Full (symbolic algebra engine)

**Recommendation:** Start minimal, extend based on user feedback

### 4. Legacy Line Integrals
**Options:**
- Integrate into composable architecture
- Remove completely
- Keep but deprecate

**Recommendation:** Remove now, reimplement properly in Phase 2 as part of multivariate

## Resource Requirements

### Developer Time
- **Phase 1:** 1 week (cleanup, docs)
- **Phase 2:** 2 weeks (cubature implementation)
- **Phase 3:** 1 week (Monte Carlo)
- **Phase 4:** 1 week (analytical preservation)
- **Phase 5:** 2 weeks (advanced features)
- **Phase 6:** 1 week (polish, release)

**Total:** ~8 weeks for full roadmap

### Testing Infrastructure
- Genz benchmark suite (standard multivariate test functions)
- Reference solutions from Mathematica/Maple
- Performance comparison with: GSL, Cuba, HIntLib

### Dependencies (Optional)
- Boost.Math (for special functions in tests)
- Eigen (for linear algebra in advanced ODE solvers)
- Intel TBB (for advanced parallelization)

## Risk Mitigation

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| Cubature too slow | High | Medium | Profile early, optimize hot paths |
| MC variance too high | Medium | Low | Use QMC, stratified sampling |
| API breaking changes | Medium | High | Semantic versioning, migration guide |
| Sparse grids too complex | Low | Medium | Make optional, good fallbacks |
| Tests too slow | Low | Medium | Parallel test execution, caching |

## Next Immediate Steps

1. **This week:** Decide on library name and scope
2. **This week:** Clean up or remove line integrals
3. **Next week:** Start Phase 2 (tensor product 2D integration)
4. **Ongoing:** Update documentation as we go

## Questions for User

1. **Name:** Keep `algebraic_integrators` or rename to `numerical_calculus`?
2. **Root finding:** In scope or defer to separate library?
3. **Priority:** Start with cubature (Phase 2) or analytical preservation (Phase 4)?
4. **Line integrals:** Remove completely or try to salvage?
5. **Release cadence:** One big v2.0 or incremental v1.x releases?
