# elasty

[![Build Status](https://travis-ci.com/yuki-koyama/elasty.svg?branch=master)](https://travis-ci.com/yuki-koyama/elasty)

A research-oriented elastic body simulator

![](./docs/cloth.gif)

## Algorithms

### Frameworks

- [x] Position-based dynamics (PBD) [Müller+07]
- [ ] Extended position-based dynamics (XPBD) [Macklin+16]
- [ ] Projective dynamics [Bouaziz+14]
- [ ] Quasi-Newton dynamics [Liu+17]

### Constraints for PBD/XPBD

- [ ] Area conservation constraint
- [x] Bending constraint [Müller+07]
- [x] Distance constraint [Müller+07]
- [x] Environmental collision constraint
- [x] Fixed point constraint
- [x] Isometric bending constraint [Bergou+06; Bender+14]
- [ ] Long range attachments constraint [Kim+12]
- [ ] Shape matching constraint [Müller+05]
- [ ] Tetrahedron strain constraint
- [ ] Triangle strain constraint
- [ ] Volume conservation constraint

## Additional Features

- Export simulated cloth meshes as Alembic

## Dependencies

### Core Library

- Alembic <https://github.com/alembic/alembic> [BSD 3-Clause]
  - OpenEXR <https://github.com/openexr/openexr> [BSD 3-Clause]
- Eigen <http://eigen.tuxfamily.org/> [MPL2]
- tinyobjloader <https://github.com/syoyo/tinyobjloader> [MIT]

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

## Contributing

Issue reports and pull requests are highly welcomed.

## References

- __[Bender+14]__ Jan Bender, Dan Koschier, Patrick Charrier, and Daniel Weber. 2014. Position-based simulation of continuous materials. Comput. Graph. 44 (2014), 1-10. DOI: http://doi.org/10.1016/j.cag.2014.07.004
- __[Bender+17]__ Jan Bender, Matthias Müller, and Miles Macklin. 2017. A survey on position based dynamics, 2017. In Proc. Eurographics '17 Tutorials, Article 6, 31 pages. DOI: https://doi.org/10.2312/egt.20171034
- __[Bergou+06]__ (TODO)
- __[Bouaziz+14]__ (TODO)
- __[Kim+12]__ (TODO)
- __[Liu+17]__ (TODO)
- __[Macklin+16]__ Miles Macklin, Matthias Müller, and Nuttapong Chentanez. 2016. XPBD: position-based simulation of compliant constrained dynamics. In Proc. MIG '16, 49-54. DOI: https://doi.org/10.1145/2994258.2994272
- __[Müller+05]__ (TODO)
- __[Müller+07]__ Matthias Müller, Bruno Heidelberger, Marcus Hennix, and John Ratcliff. 2007. Position based dynamics. J. Vis. Comun. Image Represent. 18, 2 (2007), 109-118. DOI=http://doi.org/10.1016/j.jvcir.2007.01.005
- (TODO)
