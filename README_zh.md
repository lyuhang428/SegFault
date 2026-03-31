**[English](README.md)**

# SegFault

## About

低性能，难使用的分子量子化学C++代码，目前处于非常初期阶段，仅支持闭壳层 HF, LDA, GGA 计算. 

设计逻辑是不使用输入文件（或未来使用json传递任务），用户需要将代码作为第三方库集成到计算任务中并依据需求使用相应的接口. 

## Dependencies

1. [xtensor](https://xtensor.readthedocs.io/en/latest/) (header only)
2. [libint-2.11.2](https://github.com/evaleev/libint)
3. [libxc-7.0.0](https://gitlab.com/libxc/libxc)
4. (optional but strongly recommended) intel mkl and intel icpx compiler (version 2025.3)

## Installation
```bash
source /path/to/intel/oneapi/setvars.sh
./worker.sh
cd ./install/bin/
./segfault /path/to/mol.xyz /path/to/gbs [X_id C_id]
```

## Example
see main.cc

## Contributing

该项目欢迎（急需）任何形式的贡献，包括但不限于新功能添加，性能优化，bug修复，文档更新等。
