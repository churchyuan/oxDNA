# cityu19 CUDA 验证记录（2026-10-16）

## 结论

修复下列编译和首步初始化问题后，当前工作区的源码在 cityu19 上通过 CUDA 编译与 GPU 实测。
两个圆柱几何外力在 float、mixed、double 后端均通过实际 MD 检查。
用户重启服务器后，系统驱动恢复正常；在不加载私有驱动库的环境下重跑的 **14 项检查全部通过**。
测试使用独立源码快照，没有覆盖 `/home/qiu/oxDNA` 或 `/home/qiu/oxDNA_qiu`。

## 环境与隔离目录

- 主机：cityu19；GPU：NVIDIA GeForce RTX 3080，10 GiB，计算能力 8.6。
- 重启后 NVIDIA 内核模块与系统用户态驱动库均为 **580.178.04**；默认环境 `nvidia-smi` 正常，`cuInit(0)` 返回 0，识别到 1 块 GPU。
- CUDA 编译器：`/usr/local/cuda-12.5/bin/nvcc`，12.5.82。
- 主机编译器：GCC/G++ 11.4.0；CMake 3.22.1；构建目标为 `sm_86`。
- 测试根目录：`/home/qiu/oxdna_cuda_validation/20261016-cylindrical`。
- 源码：`source/`；普通 CUDA 构建：`build-mixed/`；双精度构建：`build-double/`。
- 本地基础提交：`b7074fc3`；测试包含工作区尚未提交的修改与 NA CUDA 源文件，具体文件散列见最终清单。

## 重启前的临时处理（历史记录）

首次验证时，服务器 NVIDIA 内核模块是 **580.173.02**，用户态库是 **580.178.04**。
原环境 `nvidia-smi` 报 `Driver/library version mismatch`，`cuInit` 返回 804。
因此首次验证从 NVIDIA 官方软件仓库下载匹配的 580.173.02 用户态包，仅解压到测试目录，
通过进程级 `LD_LIBRARY_PATH` 使用。`cuInit` 随后返回 0，成功识别 GPU。
首次验证没有安装系统包、更换系统驱动、重启机器或停止既有任务。
随后用户自行重启，系统驱动已一致，当前不再需要这些私有库。
旧 `run-with-private-driver.sh` 包装器保留作历史记录，并加上内核版本检查，防止在当前驱动上加载旧库。

私有库目录：
`private-driver/usr/lib/x86_64-linux-gnu`。

包来源与 SHA-256：

| 包 | 来源 | SHA-256 |
| --- | --- | --- |
| libnvidia-compute-580 580.173.02 | [NVIDIA 官方包](https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/libnvidia-compute-580_580.173.02-1ubuntu1_amd64.deb) | `31714f52224f32ed1b33cccab5b023163d3504bb588b33decd341aa5a4b6052f` |
| libnvidia-gpucomp-580 580.173.02 | [NVIDIA 官方包](https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/libnvidia-gpucomp-580_580.173.02-1ubuntu1_amd64.deb) | `dd3d6ab8d7de42243e0e84c4d215f6b71a558a55b4d3ae54d7aab112d7f2da05` |

## 首次验证发现并修复的问题

1. `CUDADRHInteraction` 直接访问 DNA helper 的 protected `_grooving`，阻断 CUDA 编译。
   增加显式配置方法，保留原配置值。
2. `CUDARNAInteraction.cu` 中存在重复参数、缺失过滤参数、重复函数声明和不完整括号。
   按现有 CUDA kernel 签名修正；RNA 的序列平均标志使用 `_average`。
3. CPU MD 初始化没有将外力写入首步 Verlet 力，导致首次半步与 CUDA 不一致，也影响受力状态下的重启。
   初始化时先计算外力，再叠加粒子间相互作用。
4. CUDA 对照测试需要完整粒子类型，并且不能跨编译单元直接设置 `MD_N` 设备常量。
   测试现在通过真实后端的初始化方法设置常量，调用已编译库中的外力计算 kernel。

## 验证范围与结果

| 检查 | 结果 |
| --- | --- |
| `CUDA_DOUBLE=OFF` 完整 oxDNA 与测试程序编译 | 通过 |
| `CUDA_DOUBLE=ON` 完整 oxDNA 与测试程序编译 | 通过 |
| 几何、势能梯度、输入、工厂、绝对坐标与动态添加外力 | 两种构建均通过 70,328 项检查 |
| 生产 GPU 外力 kernel 与 CPU 力/势能比较 | 普通 CUDA 与 double 构建各通过 4,096 组，含空槽和多个外力叠加 |
| 侧壁、上下端盖、棱边、孔口与储液区的 MD | float、mixed、double 均通过 |
| 平移/倾斜几何、`fix_diffusion`、仍受力时的重启 | float、mixed、double 均通过 |
| CPU NVE 步长收敛与 MC 集成 | 通过；步长减半后最大能量漂移比约 0.279 |
| 原 120 核苷酸 DNA2 构型，两个外力分别运行 1,000 步 | mixed、double 均完成；轨迹与势能有限 |
| RNA2 16 核苷酸回归，序列平均与序列依赖模式 | mixed 的 matrix/edge1/edge2，以及 double 的 matrix 均通过 CPU/CUDA 状态比较 |

重启后的 14 项检查包括两种构建各自的重新配置、增量构建、CPU 单元测试、GPU 力/势能对照，
加上 CPU 集成测试、三个精度的 CUDA MD 和 mixed/double 的 DNA2/RNA2 回归。
所有进程均清除 `LD_LIBRARY_PATH` 和 `LD_PRELOAD`，直接运行测试程序，没有调用私有库包装器。
实际加载的 CUDA 驱动库为 `/usr/lib/x86_64-linux-gnu/libcuda.so.580.178.04`。
本次复测没有修改程序源码；复测前对快照清单中的 878 个文件逐一校验，全部一致。

1,000 步单粒子测试中，CPU/CUDA 最大平移状态差（位置及速度分量）为：

| 精度 | 最大绝对差 |
| --- | --- |
| mixed | 约 `1.57e-5` |
| float | 约 `1.79e-4` |
| double | 约 `9.45e-8` |

这里比较的是固定初始状态下的短轨迹，并不表示长时间混沌轨迹应逐点相等。
现有 CUDA 后端即使在 double 构建中也用 float 保存时间步长，因此轨迹差的容差与逐点力容差分别设置。

RNA2 的最大状态分量相对尺度误差约为 mixed `1.21e-5`、double `9.63e-7`。
既有后端明确拒绝 `CUDA_DOUBLE + use_edge`，该组合不在支持范围内。
本次没有对 NA 混合核酸力场做完整物理验证；NA 修复限于阻断编译的接口访问错误。

## 复现

在 cityu19 默认系统驱动环境下直接运行：

```sh
V=/home/qiu/oxdna_cuda_validation/20261016-cylindrical
"$V/build-mixed/bin/cylindrical_force_tests"
"$V/build-mixed/bin/cylindrical_cuda_tests"
"$V/build-double/bin/cylindrical_force_tests"
"$V/build-double/bin/cylindrical_cuda_tests"
python3 "$V/source/test/Forces/cylindrical_smoke.py" "$V/build-mixed/bin/oxDNA"
python3 "$V/source/test/Forces/cylindrical_cuda_smoke.py" "$V/build-mixed/bin/oxDNA" --precision mixed
python3 "$V/source/test/Forces/cylindrical_cuda_smoke.py" "$V/build-mixed/bin/oxDNA" --precision float
python3 "$V/source/test/Forces/cylindrical_cuda_smoke.py" "$V/build-double/bin/oxDNA" --precision double
python3 "$V/source/test/Forces/cuda_polymer_smoke.py" "$V/build-mixed/bin/oxDNA" --precision mixed
python3 "$V/source/test/Forces/cuda_polymer_smoke.py" "$V/build-double/bin/oxDNA" --precision double
```

构建配置使用 `CUDA=ON`、`CUDA_COMMON_ARCH=OFF`、`CMAKE_CUDA_ARCHITECTURES=86`、
`CMAKE_CUDA_COMPILER=/usr/local/cuda-12.5/bin/nvcc`、`CUDAToolkit_ROOT=/usr/local/cuda-12.5`，
以及 `NATIVE_COMPILATION=OFF`、`Python=OFF`、`ZSTD_ENABLED=OFF`；两种构建分别设置 `CUDA_DOUBLE=OFF/ON`。
重启后的配置和构建也直接使用系统驱动，无需私有库环境。

原始日志、`validation_results.json`、`polymer_validation_results.json`、
驱动包来源记录和最终源码散列保留在测试根目录；本地副本位于 `build/cityu19_cuda_validation/results/`。
重启后结果单独保存在 `post-reboot-results.json`、`post-reboot-environment.json`、`post-reboot-*.log`，
完整命令及退出状态见结果 JSON。重启前后的证据分别归档，当前汇总状态见 `final_verification.json`。
