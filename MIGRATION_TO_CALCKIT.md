# Migration Guide: algebraic_integrators → calckit

## 🎉 Welcome to calckit v2.0!

The library has been renamed from `algebraic_integrators` to **`calckit`** to better reflect its comprehensive scope and modern, composable toolkit approach.

---

## What Changed?

### Name & Branding
- **Old:** algebraic_integrators
- **New:** calckit
- **Why:** More memorable, signals toolkit/composability, broader than just "integrators"

### Namespace
```cpp
// OLD
#include "algebraic_integrators.hpp"
using namespace algebraic_integrators;
namespace ai = algebraic_integrators;

// NEW
#include "calckit.hpp"
using namespace calckit;
namespace ck = calckit;  // New short alias
```

### CMake
```cmake
# OLD
find_package(algebraic_integrators REQUIRED)
target_link_libraries(your_target ai::algebraic_integrators)

# NEW
find_package(calckit REQUIRED)
target_link_libraries(your_target ck::calckit)
```

---

## Backward Compatibility

**Good news:** Your existing code continues to work! We've provided compatibility aliases:

```cpp
// In calckit.hpp:
namespace algebraic_integrators [[deprecated]] = calckit;
namespace ai = calckit;  // Keep familiar short alias
```

### Deprecation Timeline

- **v2.0 - v2.9** (Current): Both namespaces work, old one deprecated
- **v3.0** (Future): `algebraic_integrators` alias removed

---

## Migration Checklist

### For Users

- [ ] Update `#include` statements
  ```cpp
  // Change:
  #include "algebraic_integrators.hpp"
  // To:
  #include "calckit.hpp"
  ```

- [ ] Update namespace declarations
  ```cpp
  // Change:
  using namespace algebraic_integrators;
  // To:
  using namespace calckit;
  ```

- [ ] Update CMake
  ```cmake
  find_package(calckit REQUIRED)
  target_link_libraries(your_target ck::calckit)
  ```

- [ ] Optional: Update short alias
  ```cpp
  namespace ck = calckit;  // Instead of namespace ai = ...
  ```

### For Contributors

See `RENAME_COMPLETED.md` for detailed list of files updated.

---

## What's New in v2.0?

Beyond the rename, calckit v2.0 brings major new features:

### 🆕 Multivariate Integration

```cpp
// 2D integration
auto f2d = [](std::array<double, 2> x) {
    return x[0] * x[0] + x[1] * x[1];
};
auto result = calckit::multivariate::integrate_2d(f2d, 0, 1, 0, 1, 1e-8);

// 3D integration
auto f3d = [](std::array<double, 3> x) {
    return x[0] * x[1] * x[2];
};
auto result = calckit::multivariate::integrate_3d(f3d, 0, 1, 0, 1, 0, 1, 1e-8);
```

### 🎯 Region Types

```cpp
using namespace calckit::multivariate;

// Rectangles, boxes, circles, spheres
rectangle<double> rect({0, 0}, {1, 1});
sphere<double> ball({0, 0, 0}, 1.0);

// Arbitrary regions
auto custom = arbitrary_region(bbox, [](auto x) {
    return x[0]*x[0] + x[1]*x[1] <= 1.0;  // Circle
});
```

### 📐 Tensor Product Integrators

```cpp
// Extend 1D rules to N dimensions
using Rule = calckit::quadrature::gauss_legendre_rule<double, 7>;
using Acc = calckit::accumulators::kahan_accumulator<double>;

calckit::multivariate::tensor_product_integrator<double, 2, Rule, Acc> integrator;
auto result = integrator(f, region);
```

### 🎨 Improved Concepts

```cpp
// Multivariate function concept
template <typename F, typename T, int Dim>
concept MultivariateFunction;

// Integration region concept
template <typename R, typename T, int Dim>
concept IntegrationRegion;
```

---

## Example Migration

### Before (v1.x)
```cpp
#include "algebraic_integrators.hpp"

int main() {
    using namespace algebraic_integrators;

    auto f = [](double x) { return x * x; };
    auto result = integrate_adaptive(f, 0.0, 1.0, 1e-8);

    std::cout << "Result: " << result.value << std::endl;
}
```

### After (v2.0)
```cpp
#include "calckit.hpp"

int main() {
    using namespace calckit;

    // 1D integration (unchanged API)
    auto f = [](double x) { return x * x; };
    auto result = integrate_adaptive(f, 0.0, 1.0, 1e-8);

    // NEW: 2D integration!
    auto f2d = [](std::array<double, 2> x) {
        return x[0] * x[0] + x[1] * x[1];
    };
    auto result2d = multivariate::integrate_2d(f2d, 0, 1, 0, 1, 1e-8);

    std::cout << "1D: " << result.value << std::endl;
    std::cout << "2D: " << result2d.value << std::endl;
}
```

---

## Troubleshooting

### Deprecation Warnings

If you see warnings like:
```
warning: 'algebraic_integrators' is deprecated: Use calckit namespace instead
```

**Solution:** Update to `calckit` namespace as shown above.

### CMake Can't Find calckit

**Solution:** Update your `find_package` call:
```cmake
find_package(calckit 2.0 REQUIRED)
```

### Include File Not Found

```
fatal error: algebraic_integrators.hpp: No such file or directory
```

**Solution:** Update include to:
```cpp
#include "calckit.hpp"
```

---

## Questions?

- **GitHub Issues:** https://github.com/queelius/calckit/issues
- **Discussions:** https://github.com/queelius/calckit/discussions
- **Documentation:** https://queelius.github.io/calckit/

---

## Philosophy: Why "calckit"?

1. **Kit** = Composable toolkit of calculus operations
2. **Memorable** = Easy to remember, spell, pronounce
3. **Modern** = Fresh branding for modern C++20 library
4. **Scope** = Not just integrators—full calculus toolkit

The library continues to honor Alexander Stepanov's philosophy of generic programming and algebraic structures (monoids, groups, fields), now with a name that's more accessible and memorable!

---

## Thank You!

Thank you for using calckit. We hope the new name better reflects the comprehensive, composable nature of the library, and that the new multivariate features open up new possibilities for your numerical computing needs!

Happy computing! 🚀

—The calckit Team
