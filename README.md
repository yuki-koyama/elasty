# elasty

A research-oriented elastic body simulator

## Algorithms

### Frameworks

- [x] Position-based dynamics (PBD) [Muller et al. 2007]
- [ ] Extended position-based dynamics (XPBD) [Macklin et al. 2016]
- [ ] Projective dynamics [Bouaziz et al. 2014]
- [ ] Quasi-Newton dynamics [Liu et al. 2017]

### Constraints for PBD/XPBD

- [ ] Area conservation constraint
- [ ] Bending constraint
- [x] Distance constraint
- [x] Environmental collision constraint
- [x] Fixed-point constraint
- [x] Isometric bending constraint
- [ ] Tetrahedron strain constraint
- [ ] Triangle strain constraint
- [ ] Volume conservation constraint

## Dependencies

- bigger <https://github.com/yuki-koyama/bigger> [MIT]
  - bigg <https://github.com/JoshuaBrookover/bigg> [Unlicense]
    - bgfx.cmake <https://github.com/JoshuaBrookover/bgfx.cmake> [CC0]
      - bgfx <https://github.com/bkaradzic/bgfx> [BSD 2-Clause]
      - bimg <https://github.com/bkaradzic/bimg> [BSD 2-Clause]
      - bx <https://github.com/bkaradzic/bx> [BSD 2-Clause]
    - Dear ImGui <https://github.com/ocornut/imgui> [MIT]
    - GLFW <https://github.com/glfw/glfw> [Zlib]
    - GLM <https://github.com/g-truc/glm> [MIT]
  - random-util <https://github.com/yuki-koyama/rand-util> [MIT]
  - tinyobjloader <https://github.com/syoyo/tinyobjloader> [MIT]
- Eigen <http://eigen.tuxfamily.org/> [MPL2]

## Prerequisites

```bash
brew install eigen
```

## Build

```bash
git clone https://github.com/yuki-koyama/elasty.git --recursive
mkdir build
cd build
cmake ../elasty
make
```

## License

MIT License

## References

(TODO)
