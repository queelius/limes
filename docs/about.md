# About {#about}

## Author

**Alex Towell**

- Email: [queelius@gmail.com](mailto:queelius@gmail.com)
- Blog: [metafunctor.com](https://metafunctor.com)
- GitHub: [github.com/queelius](https://github.com/queelius)

## The Project

limes grew out of a fascination with the intersection of:

- **Generic programming** — Stepanov's vision of algorithms as mathematical objects
- **Numerical analysis** — The beautiful machinery of quadrature and approximation
- **Expression templates** — C++'s ability to represent computation as types

The question that started it all: *Can we build something for integration analogous to what autograd does for differentiation?* The honest answer is no—differentiation is local (chain rule), integration is global (no analogous decomposition). But the exploration led somewhere interesting: a library where **integrals are algebraic objects** that compose according to mathematical laws.

## Philosophy

limes prioritizes **clarity over feature count**. Each component should be:

1. **Teachable** — A readable example of API design
2. **Composable** — Small pieces that combine well
3. **Inspectable** — You can see what's happening

This makes limes well-suited for:

- Learning numerical integration concepts
- Prototyping mathematical computations
- Exploring expression-based numeric programming

## Influences

- **Alex Stepanov** — Generic programming, algebraic thinking
- **From Mathematics to Generic Programming** — Stepanov & Rose
- **Elements of Programming** — Stepanov & McJones
- **Automatic Differentiation** — The expression-graph paradigm
- **Boost.Math** — High-quality numerical C++

## License

limes is released under the MIT License.

```
MIT License

Copyright (c) 2024 Alex Towell

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Contributing

Contributions are welcome! The best ways to help:

1. **Report issues** — Bug reports and feature requests
2. **Improve documentation** — Tutorials, examples, corrections
3. **Add tests** — Coverage for edge cases
4. **Suggest APIs** — Ideas for better composability

Please open issues or pull requests on GitHub.

## Acknowledgments

Thanks to the C++ community for decades of work on expression templates, concepts, and generic programming. Special thanks to the Boost.Math authors whose high-quality implementations inspired many of the numerical algorithms here.
