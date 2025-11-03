# Library Name Analysis & Recommendation

## Availability Research Results

| Name | GitHub Status | Package Conflict | Assessment |
|------|---------------|------------------|------------|
| **calculus** | ⚠️ `water-chika/calculus` exists (small, 10 stars) | Low risk | **Available** but minor conflict |
| **algebraic_calculus** | ✅ Not found | None | **Fully available** |
| **calckit** | ⚠️ Python research tool exists | Medium risk | **Available** for C++ |
| **xcalc** | ⚠️ X Window calculator | Low risk | **Available** for C++ library |
| **generic_calculus** | ✅ Not found | None | **Fully available** |
| **composable_calculus** | ✅ Not found (not searched) | None | Likely **available** |

## Detailed Analysis

### 1. **calculus**
**Availability:** ⚠️ Moderate conflict
- Small existing C++ library: `water-chika/calculus` (10 stars, minimal activity)
- General term, likely many small projects with this name
- Would need disambiguation (e.g., on package managers)

**Pros:**
- ⭐ Simple, one word, memorable
- ⭐ Perfectly descriptive
- ⭐ Easy to explain

**Cons:**
- ⚠️ Generic term, hard to Google
- ⚠️ Conflicts with existing (small) project
- ⚠️ May confuse with other "calculus" libraries

**Verdict:** Good if existing conflicts acceptable

---

### 2. **algebraic_calculus** ⭐ RECOMMENDED
**Availability:** ✅ Fully available
- No C++ library with this name found
- No package manager conflicts

**Pros:**
- ⭐ Honors Stepanov philosophy explicitly
- ⭐ "Algebraic" = structures/generic programming (clear to CS audience)
- ⭐ Distinguishes from "numerical_calculus"
- ⭐ Maintains continuity with current name
- ⭐ Namespace: `algebraic_calculus::` or `acalc::`
- ⭐ Educational: teaches about algebraic thinking

**Cons:**
- Longer name (two words)
- Could be confused with symbolic algebra (but documentation clarifies)

**Verdict:** Best balance of clarity, availability, and philosophy

---

### 3. **calckit**
**Availability:** ⚠️ Moderate conflict
- "Calkit" (Python research automation tool) exists
- "CalcKit" (mobile calculator app) exists
- Different domains, but some name collision

**Pros:**
- ⭐ Modern, memorable portmanteau
- ⭐ "Kit" implies composability/toolkit
- ⭐ Single word, easy to type
- ⭐ Fun, approachable name

**Cons:**
- ⚠️ Possible confusion with existing tools
- Doesn't convey Stepanov philosophy
- Less academic/serious sounding

**Verdict:** Good alternative if want modern branding

---

### 4. **xcalc**
**Availability:** ⚠️ Low conflict
- Primarily X Window System calculator
- Some small GitHub calculator projects
- No major C++ library

**Pros:**
- ⭐ Follows Boost convention (xpressive, xeus, etc.)
- ⭐ "x" = extensible/generic
- ⭐ Short, memorable

**Cons:**
- ⚠️ Conflicts with X Window calculator (known tool)
- Doesn't convey calculus focus
- "x" meaning unclear without context

**Verdict:** Possible but not ideal

---

### 5. **generic_calculus**
**Availability:** ✅ Fully available

**Pros:**
- ⭐ Explicitly conveys generic programming
- ⭐ Clear connection to STL philosophy
- ⭐ No conflicts found

**Cons:**
- "Generic" sounds vague to non-CS audience
- Longer name
- Less elegant than "algebraic"

**Verdict:** Solid alternative to algebraic_calculus

---

### 6. **composable_calculus**
**Availability:** ✅ Likely available (not searched)

**Pros:**
- ⭐ Emphasizes key architecture feature
- ⭐ Clear value proposition
- ⭐ Approachable term

**Cons:**
- Long name
- "Composable" less precise than "algebraic"
- Doesn't convey Stepanov connection

**Verdict:** Descriptive but verbose

---

## Comparison Matrix

| Criterion | calculus | algebraic_calculus | calckit | xcalc | generic_calculus |
|-----------|----------|-------------------|---------|-------|------------------|
| **Availability** | ⚠️ Moderate | ✅ Excellent | ⚠️ Moderate | ⚠️ Low | ✅ Excellent |
| **Clarity** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐ |
| **Stepanov Philosophy** | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐ | ⭐ | ⭐⭐⭐⭐ |
| **Memorability** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ |
| **SEO/Googleability** | ⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐ |
| **Composability Signal** | ⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐ |
| **Academic Appeal** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐ |
| **Brevity** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **TOTAL** | 23/35 | **31/35** | 23/35 | 15/35 | 25/35 |

---

## Top 3 Recommendations

### 🥇 **algebraic_calculus** (BEST OVERALL)

**Why this wins:**
1. ✅ **Fully available** - No conflicts
2. 🎓 **Honors Stepanov** - "Algebraic" = structures, generic programming
3. 🔍 **Highly Googleable** - Unique, specific
4. 📚 **Educational value** - Teaches algebraic thinking
5. 🔄 **Maintains identity** - Evolution of current name
6. 🎯 **Scope flexibility** - Works for numeric + analytical hybrid

**Namespace options:**
```cpp
namespace algebraic_calculus { }  // Full
namespace acalc { }               // Abbreviated alias
```

**Tagline:**
> "Algebraic Calculus - Generic programming for numerical analysis"

---

### 🥈 **calckit** (MODERN ALTERNATIVE)

**Why this is second:**
1. ⭐ Modern, approachable branding
2. ⭐ "Kit" signals composability
3. ⚠️ Minor conflicts (different domains)
4. ✨ Fresh start, no baggage

**Good if:** You want mass appeal over academic precision

---

### 🥉 **generic_calculus** (SAFE CHOICE)

**Why this is third:**
1. ✅ Fully available
2. 📖 Clear generic programming connection
3. 🎓 Academic but less elegant than "algebraic"

**Good if:** Want explicit "generic" in name

---

## Final Recommendation

**Choose: `algebraic_calculus`**

### Rationale

1. **Philosophy alignment**: "Algebraic" perfectly captures Stepanov's focus on algebraic structures (monoids, groups, fields) without implying symbolic manipulation

2. **Clear differentiation**:
   - NOT "numerical_calculus" (we do analytical preservation too)
   - NOT "symbolic_calculus" (we use numerics as fallback)
   - YES "algebraic_calculus" (generic programming with structures)

3. **Educational mission**: Users learn that "algebraic" means:
   - Composable operations
   - Structural thinking
   - Generic programming
   - Stepanov's vision

4. **Marketing pitch**:
   > "From the creators of algebraic_integrators comes algebraic_calculus -
   > a comprehensive C++20 library for calculus using generic programming and
   > algebraic structures. Compose integrators, differentiators, and solvers
   > with the elegance of Stepanov's STL philosophy."

5. **Package names**:
   - GitHub: `algebraic_calculus`
   - Conan: `algebraic-calculus`
   - vcpkg: `algebraic-calculus`
   - Header: `<algebraic_calculus.hpp>`
   - Namespace: `algebraic_calculus::` or `acalc::`

6. **Documentation clarity**:
   ```markdown
   # Algebraic Calculus

   **Why "algebraic"?**
   Following Alexander Stepanov's philosophy, we use "algebraic" to mean
   programming with algebraic structures like monoids, groups, and fields.

   This is NOT symbolic algebra (though we preserve analytical relationships).
   This IS generic programming applied to calculus operations.
   ```

---

## Alternative: If Two Words is Too Long

If `algebraic_calculus` feels too long, consider **portmanteau**:

- **`algcalc`** - Algebraic Calculus compressed
- **`alcalc`** - AL-gebra + CALC-ulus
- **`calculus`** - Just go for it, document the conflict

But honestly, `algebraic_calculus` with alias `acalc` is perfect.

---

## Migration Strategy

If choosing `algebraic_calculus`:

1. **GitHub rename** - Preserve stars/forks
2. **Namespace alias** - Keep `algebraic_integrators` as deprecated alias
3. **CMake** - Support both names for 1-2 versions
4. **Documentation** - Clear migration guide
5. **Announcement** - Explain philosophy behind name

```cpp
// Backward compatibility
namespace algebraic_integrators = algebraic_calculus;
namespace ai = algebraic_calculus;  // Keep short alias
```

---

## Decision Time

**Vote:**
- [ ] `algebraic_calculus` ⭐ RECOMMENDED
- [ ] `calckit` (modern alternative)
- [ ] `generic_calculus` (explicit alternative)
- [ ] `calculus` (simple but conflicted)
- [ ] Other: _______________

What's your choice?
