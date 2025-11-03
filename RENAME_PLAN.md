# Library Rename Execution Plan

## Decision: [PENDING USER CHOICE]

Top candidates:
1. **algebraic_calculus** ⭐ (recommended)
2. **calckit** (modern alternative)
3. **calculus** (simple, conflicts)

## Files Requiring Updates

### Core Headers
- [ ] `include/algebraic_integrators.hpp` → rename and update namespace
- [ ] `include/concepts/integrator_concepts.hpp` → update namespace
- [ ] All `include/**/*.hpp` → update namespace declarations

### CMake Build System
- [ ] `CMakeLists.txt`:
  - `project(algebraic_integrators ...)` → `project(NEW_NAME ...)`
  - `add_library(algebraic_integrators ...)` → `add_library(NEW_NAME ...)`
  - `namespace ai::` → `namespace NEW_ALIAS::`
- [ ] `cmake/*.cmake` → update package names
- [ ] `tests/CMakeLists.txt` → update target links

### Documentation
- [ ] `README.md` → rewrite with new name and clearmission
- [ ] `CLAUDE.md` → update all references
- [ ] `docs/index.md` → update branding
- [ ] `docs/getting-started.md` → update includes, namespaces
- [ ] `docs/architecture.md` → update terminology
- [ ] `docs/api/*.md` → update code examples
- [ ] `docs/examples.md` → update all code snippets
- [ ] `mkdocs.yml` → update site_name

### Test Files
- [ ] `tests/test_*.cpp` → update includes and namespaces
- [ ] Test assertions that check namespace names

### Examples
- [ ] `examples/demo.cpp` → update includes
- [ ] `examples/simple_demo.cpp` → update includes

### Metadata
- [ ] `.gitignore` → no changes needed
- [ ] LICENSE → no changes (copyright remains)

### New Files Created During This Session
- [ ] `include/concepts/multivariate_concepts.hpp` → uses `algebraic_calculus` namespace
- [ ] `include/multivariate/regions.hpp` → uses `algebraic_calculus` namespace
- [ ] `include/multivariate/tensor_product_integrator.hpp` → uses `algebraic_calculus` namespace

**Note:** These already use `algebraic_calculus` - may need to change if different name chosen!

---

## Namespace Strategy

### Option A: Complete Rename
```cpp
// Old
namespace algebraic_integrators { }

// New
namespace algebraic_calculus { }  // or chosen name
namespace acalc { }  // short alias
```

### Option B: Alias for Compatibility (Recommended)
```cpp
// New primary namespace
namespace algebraic_calculus {
    // All new code here
}

// Backward compatibility (deprecate in v3.0)
namespace algebraic_integrators = algebraic_calculus;

// Short aliases
namespace acalc = algebraic_calculus;
namespace ai = algebraic_calculus;  // Keep familiar alias
```

**Recommendation:** Option B for smooth migration

---

## Migration Timeline

### Version 2.0.0 (This Release)
- ✅ Introduce new name `algebraic_calculus`
- ✅ Add namespace alias: `namespace algebraic_integrators = algebraic_calculus;`
- ✅ Update all documentation to use new name
- ✅ Examples use new name
- ⚠️ Old namespace still works (deprecated with warnings)

### Version 2.x (Transition Period)
- Both namespaces fully supported
- Documentation emphasizes new name
- Deprecation warnings in headers:
  ```cpp
  [[deprecated("Use algebraic_calculus instead")]]
  namespace algebraic_integrators = algebraic_calculus;
  ```

### Version 3.0.0 (Future)
- Remove `algebraic_integrators` alias
- Clean break, new name only

---

## Detailed Update Script

### Step 1: Rename Repository (GitHub)
```bash
# GitHub Settings → Repository name
# algebraic_integrators → algebraic_calculus
```

### Step 2: Update CMakeLists.txt
```cmake
# OLD:
project(algebraic_integrators VERSION 2.0.0 LANGUAGES CXX)
add_library(algebraic_integrators INTERFACE)
add_library(ai::algebraic_integrators ALIAS algebraic_integrators)

# NEW:
project(algebraic_calculus VERSION 2.0.0 LANGUAGES CXX)
add_library(algebraic_calculus INTERFACE)
add_library(acalc::algebraic_calculus ALIAS algebraic_calculus)
add_library(ai::algebraic_calculus ALIAS algebraic_calculus)  # Compat alias
```

### Step 3: Update Main Header
```cpp
// OLD: include/algebraic_integrators.hpp
// NEW: include/algebraic_calculus.hpp

#pragma once

// Primary namespace
namespace algebraic_calculus {
    // ... all code ...
}

// Aliases
namespace acalc = algebraic_calculus;

// Backward compatibility (deprecated)
#ifdef ALGEBRAIC_CALCULUS_ENABLE_LEGACY_NAMESPACE
namespace algebraic_integrators [[deprecated("Use algebraic_calculus")]]
    = algebraic_calculus;
namespace ai = algebraic_calculus;
#endif
```

### Step 4: Update All Headers
```bash
# Find and replace in all headers
find include -name "*.hpp" -exec sed -i \
  's/namespace algebraic_integrators/namespace algebraic_calculus/g' {} \;
```

### Step 5: Update Documentation
```bash
# Update all markdown files
find . -name "*.md" -exec sed -i \
  's/algebraic_integrators/algebraic_calculus/g' {} \;

# Update mkdocs.yml
sed -i 's/Algebraic Integrators/Algebraic Calculus/g' mkdocs.yml
```

### Step 6: Update Tests
```bash
# Update test files
find tests -name "*.cpp" -exec sed -i \
  's/algebraic_integrators/algebraic_calculus/g' {} \;
```

### Step 7: Git Commit
```bash
git mv include/algebraic_integrators.hpp include/algebraic_calculus.hpp
git add -A
git commit -m "Rename library to algebraic_calculus

- Reflects broader scope (integration, differentiation, ODEs, root finding)
- Honors Stepanov's philosophy of algebraic structures
- Maintains backward compatibility via namespace alias
- Update all documentation and examples

Breaking changes:
- Primary namespace: algebraic_integrators → algebraic_calculus
- Main header: algebraic_integrators.hpp → algebraic_calculus.hpp
- CMake target: algebraic_integrators → algebraic_calculus

Migration:
- Old namespace aliased for compatibility (will be deprecated in v3.0)
- Update includes: #include \"algebraic_calculus.hpp\"
- Update using: using namespace algebraic_calculus; or namespace acalc;
"
```

### Step 8: Update Package Metadata
```bash
# Update conan package
# Update vcpkg port
# Update any other package managers
```

### Step 9: Deploy Updated Documentation
```bash
mkdocs gh-deploy
```

### Step 10: Create Release
```bash
git tag -a v2.0.0 -m "Version 2.0.0: Renamed to algebraic_calculus"
git push origin v2.0.0
gh release create v2.0.0 --title "v2.0.0 - Algebraic Calculus" \
  --notes "See CHANGELOG.md for details"
```

---

## User Communication

### Announcement Template

```markdown
# 🎉 algebraic_integrators → algebraic_calculus

We're excited to announce the evolution of `algebraic_integrators` to **`algebraic_calculus`**!

## Why the Change?

The library has grown beyond just integration:
- ✅ Numerical integration (1D and multivariate)
- ✅ Numerical differentiation
- ✅ ODE solvers
- ✅ Root finding (coming soon)
- ✅ Antiderivatives with analytical preservation

The name "integrators" was too narrow. "Algebraic calculus" better reflects:
1. **Comprehensive calculus operations** (not just integration)
2. **Stepanov's philosophy** (algebraic structures, generic programming)
3. **Hybrid approach** (analytical + numerical)

## Migration Guide

### Backward Compatibility
Your existing code continues to work! We've aliased the old namespace:
```cpp
namespace algebraic_integrators = algebraic_calculus;  // Deprecated
```

### Recommended Updates
```cpp
// OLD
#include "algebraic_integrators.hpp"
using namespace algebraic_integrators;

// NEW
#include "algebraic_calculus.hpp"
using namespace algebraic_calculus;
// Or use short alias:
namespace acalc = algebraic_calculus;
```

### CMake
```cmake
# OLD
find_package(algebraic_integrators REQUIRED)
target_link_libraries(your_target ai::algebraic_integrators)

# NEW
find_package(algebraic_calculus REQUIRED)
target_link_libraries(your_target acalc::algebraic_calculus)
```

## What's New in v2.0?

- 🚀 **Multivariate Integration** - 2D/3D cubature and Monte Carlo
- 🎯 **Improved API** - Cleaner, more composable
- 📚 **Better Documentation** - Complete rewrite
- 🧪 **More Tests** - 150+ unit tests
- ⚡ **Performance** - Optimizations throughout

## Timeline

- **v2.x**: Both namespaces supported (old deprecated)
- **v3.0** (future): Old namespace removed

Questions? Open an issue or discussion!
```

---

## Checklist

Before executing rename:
- [ ] User confirms chosen name
- [ ] Backup current repository state
- [ ] Create migration branch
- [ ] Test builds with new name
- [ ] Update all documentation
- [ ] Test backward compatibility
- [ ] Prepare announcement
- [ ] Update package managers
- [ ] Coordinate with any dependent projects

After rename:
- [ ] GitHub repository renamed
- [ ] Documentation deployed
- [ ] Release created
- [ ] Announcement published
- [ ] Package managers updated
- [ ] Users notified

---

## Rollback Plan

If issues arise:
1. Revert git commits
2. Restore old GitHub repository name
3. Redeploy old documentation
4. Issue apology and explanation
5. Plan more carefully for next attempt

---

## Questions for User

Before proceeding, confirm:
1. **Final name choice:** `algebraic_calculus`, `calckit`, or other?
2. **Version number:** 2.0.0 or different?
3. **Deprecation timeline:** Remove old namespace in v3.0?
4. **Backward compat:** Support both namespaces in v2.x?
5. **GitHub org:** Keep under queelius or move to new org?

---

Ready to execute when user confirms decision!
