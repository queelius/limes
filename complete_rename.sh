#!/bin/bash
# Complete Rename Script: algebraic_integrators → calckit
# This script updates all remaining files to use the calckit namespace

set -e  # Exit on error

echo "🚀 Completing rename to calckit..."

# Update all header files (namespace declarations)
echo "📝 Updating header namespaces..."
find include -name "*.hpp" -type f -exec sed -i \
    's/namespace algebraic_integrators/namespace calckit/g' {} \;

find include -name "*.hpp" -type f -exec sed -i \
    's/algebraic_integrators::/calckit::/g' {} \;

# Update test files
echo "🧪 Updating test files..."
find tests -name "*.cpp" -type f -exec sed -i \
    's/#include "algebraic_integrators.hpp"/#include "calckit.hpp"/g' {} \;

find tests -name "*.cpp" -type f -exec sed -i \
    's/using namespace algebraic_integrators/using namespace calckit/g' {} \;

find tests -name "*.cpp" -type f -exec sed -i \
    's/algebraic_integrators::/calckit::/g' {} \;

# Update test CMakeLists.txt
echo "📦 Updating test CMakeLists.txt..."
sed -i 's/algebraic_integrators/calckit/g' tests/CMakeLists.txt

# Update examples
echo "💡 Updating examples..."
find examples -name "*.cpp" -type f -exec sed -i \
    's/#include "algebraic_integrators.hpp"/#include "calckit.hpp"/g' {} \;

find examples -name "*.cpp" -type f -exec sed -i \
    's/using namespace algebraic_integrators/using namespace calckit/g' {} \;

find examples -name "*.cpp" -type f -exec sed -i \
    's/algebraic_integrators::/calckit::/g' {} \;

# Update examples CMakeLists.txt
sed -i 's/algebraic_integrators/calckit/g' examples/CMakeLists.txt

# Update documentation
echo "📚 Updating documentation..."
find docs -name "*.md" -type f -exec sed -i \
    's/algebraic_integrators/calckit/g' {} \;

sed -i 's/Algebraic Integrators/CalcKit/g' docs/index.md
sed -i 's/algebraic_integrators/calckit/g' mkdocs.yml
sed -i 's/site_name: .*/site_name: CalcKit/' mkdocs.yml

# Update README
echo "📖 Updating README..."
sed -i 's/algebraic_integrators/calckit/g' README.md
sed -i 's/Algebraic integrators and differentiators/CalcKit - Composable Calculus Toolkit/g' README.md

# Update CLAUDE.md
echo "📋 Updating CLAUDE.md..."
sed -i 's/algebraic_integrators/calckit/g' CLAUDE.md
sed -i 's/algebraic_integrators/calckit/g' CLAUDE.md
sed -i 's/algebraic-integrators/calckit/g' CLAUDE.md

# Update .gitignore if needed
if grep -q "algebraic_integrators" .gitignore 2>/dev/null; then
    echo "🔧 Updating .gitignore..."
    sed -i 's/algebraic_integrators/calckit/g' .gitignore
fi

echo "✅ Rename complete!"
echo ""
echo "Summary of changes:"
echo "  - Namespace: algebraic_integrators → calckit"
echo "  - Main header: algebraic_integrators.hpp → calckit.hpp"
echo "  - CMake target: algebraic_integrators → calckit"
echo "  - Version: 1.0.0 → 2.0.0"
echo ""
echo "Next steps:"
echo "1. Test build: cmake -S . -B build && cmake --build build"
echo "2. Run tests: cd build && ctest"
echo "3. Review changes: git diff"
echo "4. Commit: git add -A && git commit"
echo ""
echo "🎉 Welcome to calckit v2.0!"
