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
- [x] Bending constraint
- [x] Distance constraint
- [x] Environmental collision constraint
- [x] Fixed-point constraint
- [x] Isometric bending constraint
- [ ] Long-range attachment constraint
- [ ] Tetrahedron strain constraint
- [ ] Triangle strain constraint
- [ ] Volume conservation constraint

## Additional Features

- Export simulated cloth meshes as Alembic

## Dependencies

### Core

- Alembic <https://github.com/alembic/alembic> [BSD 3-Clause]
  - OpenEXR <https://github.com/openexr/openexr> [BSD 3-Clause]
- Eigen <http://eigen.tuxfamily.org/> [MPL2]

### Demos

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
  - string-util <https://github.com/yuki-koyama/string-util> [MIT]
  - tinyobjloader <https://github.com/syoyo/tinyobjloader> [MIT]
- timer <https://github.com/yuki-koyama/timer> [MIT]

## Prerequisites

```bash
brew install eigen openexr
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

- Jan Bender, Dan Koschier, Patrick Charrier, and Daniel Weber. 2014. Position-based simulation of continuous materials. Comput. Graph. 44 (2014), 1-10. DOI: http://dx.doi.org/10.1016/j.cag.2014.07.004
- Jan Bender, Matthias Müller, and Miles Macklin. 2017. A survey on position based dynamics, 2017. In Proc. Eurographics '17 Tutorials, Article 6, 31 pages. DOI: https://doi.org/10.2312/egt.20171034
- Miles Macklin, Matthias Müller, and Nuttapong Chentanez. 2016. XPBD: position-based simulation of compliant constrained dynamics. In Proc. MIG '16, 49-54. DOI: https://doi.org/10.1145/2994258.2994272
- Matthias Müller, Bruno Heidelberger, Marcus Hennix, and John Ratcliff. 2007. Position based dynamics. J. Vis. Comun. Image Represent. 18, 2 (2007), 109-118. DOI=http://dx.doi.org/10.1016/j.jvcir.2007.01.005
- (TODO)
