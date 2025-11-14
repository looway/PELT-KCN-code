# ARxAWEN — Radar AWEN Processing Scripts

简要说明
- 本目录包含用于仿真 FMCW 雷达信号并比较多种去幅包络（WEN / AWEN / AR）处理和 CFAR 检测的 MATLAB 脚本。
- 主要入口脚本：[`ARxAWEN.m`](ARxAWEN.m)。

主要文件
- [`ARxAWEN.m`](ARxAWEN.m) — 主脚本：生成仿真信号、应用 WEN/AWEN/AR 处理、计算频谱并用 CFAR 检测目标，最后保存多个图像。
- [`myCFAR.m`](myCFAR.m) — 项目自定义 CFAR 检测器（在主脚本中以 [`myCFAR`](myCFAR.m) 调用）。
- [`WEN_AWEN.m`](WEN_AWEN.m) — 与 WEN/AWEN 处理相关的示例/辅助代码。
- [`myfilgap.m`](myfilgap.m) — 用于缺失段填充/插值的自定义函数（在 AR 处理中有类似用途）。
- [`awen1.m`](awen1.m), [`awen2.m`](awen2.m), [`real_test.m`](real_test.m) — 相关实验/对比脚本和不同参数的实现示例。
- 其它示例：[`compare_test.m`](compare_test.m), [`atest.m`](atest.m)。

依赖
- MATLAB（建议含 Signal Processing Toolbox / Phased Array System Toolbox，根据脚本中对 `physconst`, `findchangepts` 等的调用）
- 脚本中使用的自定义函数：[`myCFAR`](myCFAR.m)、[`myfilgap`](myfilgap.m)、以及项目内其它辅助脚本（请确保工作目录或 MATLAB 路径包含本目录）。

快速运行
1. 在 MATLAB 中将当前目录切换到本文件夹：
2. 运行主脚本：
   - 在命令行输入: run('ARxAWEN.m') 或直接在编辑器中运行 [`ARxAWEN.m`](ARxAWEN.m)。
3. 脚本会生成若干图窗并保存为 1.jpg … 11.jpg（见脚本末尾的 saveas 调用）。

可调整的关键参数（在 [`ARxAWEN.m`](ARxAWEN.m) 中）
- rangeMax — 最大检测距离（m）
- rangeRes — 距离分辨率（m）
- A / AI — 仿真信号的幅值参数
- guard, train, Pf — CFAR 的 guard cells、training cells 和期望虚警率
- beta — AWEN / AR 中用于计算阈值的比例因子

输出
- 多个图形窗口（WEN / AWEN / AR 时间域与距离/频谱图）
- 保存文件：1.jpg … 11.jpg（保存在当前目录）

注意事项
- 若发生未定义函数/变量错误，请确保所有自定义脚本（例如 [`myCFAR.m`](myCFAR.m)、[`myfilgap.m`](myfilgap.m) 等）已在 MATLAB 路径中。
- 部分函数（如 `findchangepts`）需要 Signal Processing Toolbox。
- 脚本包含大量绘图与参数，可据需求去注释/调整。

参考与调试
- 查看主脚本：[`ARxAWEN.m`](ARxAWEN.m)
- 检查并调试 CFAR 实现：[`myCFAR.m`](myCFAR.m)
- WEN/AWEN 示例与对比：[`WEN_AWEN.m`](WEN_AWEN.m), [`awen1.m`](awen1.m), [`awen2.m`](awen2.m), [`real_test.m`](real_test.m)
