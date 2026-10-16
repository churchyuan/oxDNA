# 圆柱限域与纳米孔

这组示例用一个核苷酸检查外力接入。初始位置有意靠近侧壁，`particle_radius=0.4`
使它产生可测的排斥；这不是实际 DNA 穿孔的生产模拟，也没有施加牵引或电场。

从本目录运行：

```sh
../../build/bin/oxDNA input
../../build/bin/oxDNA input external_forces_file=nanopore.force
```

两次运行使用相同输出文件名。若需要保留两份轨迹，请在不同工作目录复制示例后运行。
使用已编译 CUDA 的 oxDNA 时，可以追加 `backend=CUDA backend_precision=mixed`。
`CUDA_DOUBLE=ON` 构建应使用 `backend_precision=double`。实际 GPU 验证可运行
`cmake --build build --target test_cylindrical_cuda`。

| 选项 | nanopore | cylinder_confinement |
| --- | --- | --- |
| 几何 | 有圆形通孔的无限横向膜 | 侧壁与上下端盖组成的封闭圆柱 |
| 轴向尺寸 | `thickness`，膜的完整厚度 | `height`，容器的完整高度 |
| 粒子半径 | 按球体到膜实体的距离计算，包含孔口棱边 | 可达半径与半高分别减去粒子半径 |
| 两端 | 开放，连通储液区 | 有端盖 |
| 排斥 | 到膜实体的有符号距离的谐和罚势 | 侧壁与端盖违反量的平方和 |

`center` 是几何中心，`dir` 是圆柱轴，自动归一化；尺寸、刚度和向量会严格校验。
这两种外力都按绝对坐标计算，不对外场做最小镜像，受力分子不会被 `fix_diffusion`
移回周期盒。DNA 自身的相互作用仍遵循模拟盒边界条件，应避免不希望的周期镜像接触。
当前支持固定盒的 NVT MC 与 MD；NPT/`use_barostat` 会明确报错，因为现有后端的
体积移动接受率没有一致地计入外场势能。

`particle_radius` 默认是 0，表示对核苷酸质心限域；0.4 只是本示例的显式选择，
不是 oxDNA 自动提供的真实核苷酸半径。它用各向同性球体近似排斥，不包含核苷酸取向，
也不直接检查两个核苷酸之间的键段是否穿过膜。

两者都是有限刚度的软排斥。刚度越大，越需要缩小时间步长；应检查穿透深度及步长收敛。
nanopore 的膜穿透势垒有限，尤其不能把增加刚度当作绝对不可穿膜的保证。
旧的 nanopore 输入仍可使用；默认参数保留原势能，修正了对称面和最近边界等距位置
任意选择正轴方向的行为。膜实体内的距离函数仍存在不可微位置，不能在那里声称力全局连续。

详细数学定义、边界处理、参数范围和测试说明见
[外力文档](../../docs/source/forces.md#cylindrical-geometry-nanopore-and-closed-confinement)。

已在 cityu19 的 RTX 3080 / CUDA 12.5 上通过 float、mixed、double 验证，包含
120 核苷酸 DNA2 构型的实际 MD。环境问题、修复内容和复现命令见
[CUDA 验证记录](../../test/Forces/VALIDATION_cityu19_20261016.md)。
