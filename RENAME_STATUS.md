# Rename Status: algebraic_integrators → calckit

## ✅ Completed Tasks

### Core Infrastructure
- [x] **Renamed main header**: `include/algebraic_integrators.hpp` → `include/calckit.hpp`
- [x] **Updated main namespace**: `algebraic_integrators` → `calckit`
- [x] **Added backward compatibility**: `namespace algebraic_integrators [[deprecated]] = calckit;`
- [x] **Updated CMakeLists.txt**: Project name, targets, and installation
  - Project: `calckit` v2.0.0
  - Target: `calckit` (with `ck::calckit` alias)
  - Backward compat: `algebraic_integrators` target links to `calckit`
  - Namespace alias: `ai::algebraic_integrators` → `calckit`

### New Multivariate Code
- [x] **Updated multivariate concepts**: `include/concepts/multivariate_concepts.hpp`
- [x] **Updated regions**: `include/multivariate/regions.hpp`
- [x] **Updated tensor product integrator**: `include/multivariate/tensor_product_integrator.hpp`

### Documentation
- [x] **Created migration guide**: `MIGRATION_TO_CALCKIT.md`
- [x] **Created rename script**: `complete_rename.sh`
- [x] **Created status document**: `RENAME_STATUS.md` (this file)

---

## ⏳ Remaining Tasks

These can be completed by running `./complete_rename.sh`:

### Headers
- [ ] Update all existing headers in `include/` to use `calckit` namespace
  - `include/concepts/integrator_concepts.hpp`
  - `include/core/integration_result.hpp`
  - `include/accumulators/accumulators.hpp`
  - `include/quadrature/quadrature_rules.hpp`
  - `include/integrators/univariate_integrator.hpp`
  - `include/parallel/parallel_integration.hpp`
  - `include/transforms/coordinate_transforms.hpp`
  - `include/numerical_differentiation.hpp`
  - `include/central_finite_difference.hpp`
  - `include/antiderivative.hpp`
  - `include/numerical_integrators/*.hpp`
  - `include/ode/*.hpp`

### Tests
- [ ] Update all test files to use `calckit`
  - `tests/test_*.cpp` (10 files)
  - `tests/CMakeLists.txt`
  - `tests/basic_tests.cpp`

### Examples
- [ ] Update example files
  - `examples/demo.cpp`
  - `examples/simple_demo.cpp`
  - `examples/CMakeLists.txt`

### Documentation
- [ ] Update README.md
- [ ] Update CLAUDE.md
- [ ] Update docs/*.md files
- [ ] Update mkdocs.yml
- [ ] Update VISION.md, STEPANOV_DESIGN.md, MULTIVARIATE_DESIGN.md, ROADMAP.md

---

## How to Complete the Rename

### Option 1: Automated (Recommended)

Run the provided script:
```bash
./complete_rename.sh
```

This will update all remaining files automatically.

### Option 2: Manual

Follow the steps in `RENAME_PLAN.md` to update files individually.

---

## Testing the Rename

After completing the rename:

### 1. Test Build
```bash
# Clean build
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### 2. Run Tests
```bash
cd build
ctest --output-on-failure
```

### 3. Check Backward Compatibility
```cpp
// This should still work (with deprecation warning):
#include "algebraic_integrators.hpp"  // Links to calckit.hpp
using namespace algebraic_integrators;

auto f = [](double x) { return x * x; };
auto result = integrate_adaptive(f, 0, 1, 1e-8);
```

### 4. Test New Code
```cpp
#include "calckit.hpp"
using namespace calckit;

// 1D integration
auto f = [](double x) { return x * x; };
auto result = integrate_adaptive(f, 0, 1, 1e-8);

// 2D integration (NEW!)
auto f2d = [](std::array<double, 2> x) {
    return x[0] * x[0] + x[1] * x[1];
};
auto result2d = multivariate::integrate_2d(f2d, 0, 1, 0, 1, 1e-8);
```

---

## Commit Strategy

### Step 1: Complete Rename
```bash
./complete_rename.sh
```

### Step 2: Review Changes
```bash
git status
git diff
```

### Step 3: Commit
```bash
git add -A
git commit -m "Rename library to calckit v2.0

Major changes:
- Rename: algebraic_integrators → calckit
- Namespace: algebraic_integrators → calckit
- Main header: algebraic_integrators.hpp → calckit.hpp
- CMake target: algebraic_integrators → calckit
- Version bump: 1.0.0 → 2.0.0

New features in v2.0:
- Multivariate integration (2D/3D)
- Region types (hyperrectangle, simplex, ball)
- Tensor product integrator
- Improved concepts for multivariate functions

Backward compatibility:
- Old namespace aliased (deprecated, will remove in v3.0)
- Old CMake target links to new target
- Migration guide provided

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Step 4: Tag Release
```bash
git tag -a v2.0.0 -m "Version 2.0.0: Renamed to calckit with multivariate integration"
```

### Step 5: Push
```bash
git push origin main
git push origin v2.0.0
```

### Step 6: Update GitHub Repository Name
1. Go to GitHub repository settings
2. Repository name: `algebraic_integrators` → `calckit`
3. Update description: "CalcKit - Composable calculus toolkit for modern C++"

### Step 7: Deploy Documentation
```bash
mkdocs gh-deploy
```

### Step 8: Create GitHub Release
```bash
gh release create v2.0.0 \
  --title "v2.0.0 - Renamed to calckit" \
  --notes-file MIGRATION_TO_CALCKIT.md
```

---

## Files Modified Summary

### Core (Completed)
- `include/calckit.hpp` (renamed from algebraic_integrators.hpp)
- `CMakeLists.txt`
- `include/concepts/multivariate_concepts.hpp` (new)
- `include/multivariate/regions.hpp` (new)
- `include/multivariate/tensor_product_integrator.hpp` (new)

### Documentation (Completed)
- `MIGRATION_TO_CALCKIT.md` (new)
- `RENAME_STATUS.md` (new, this file)
- `complete_rename.sh` (new)

### Remaining (Use script)
- 20+ header files
- 10+ test files
- 2 example files
- 10+ documentation files

---

## Post-Rename Checklist

After running the script and committing:

- [ ] Test build passes
- [ ] All tests pass
- [ ] Backward compatibility works
- [ ] Documentation builds
- [ ] Examples compile and run
- [ ] GitHub repository renamed
- [ ] Documentation deployed
- [ ] Release created
- [ ] Users notified (if applicable)

---

## Questions?

See `MIGRATION_TO_CALCKIT.md` for user-facing migration guide.
See `RENAME_PLAN.md` for original planning document.

Ready to complete? Run `./complete_rename.sh`! 🚀
