# FASTQ_Compressor

**FASTQ_Compressor** 是一个支持两种压缩策略的 FASTQ 文件压缩与解压工具，目标是高效地处理 DNA 测序数据。它可以将 FASTQ 文件压缩为一个包含序列、质量信息与位置信息的 `.tar.gz` 包，并能恢复原始数据。

---

## 🚀 特性

- ✅ 支持两种 k-mer 压缩策略（标准与优化）与一种类二进制优化压缩法
- ✅ 使用 Biopython 处理 FASTQ 格式
- ✅ 命令行工具 `fastq-compress`，使用简单
- ✅ 解压时自动判断或手动指定压缩策略
- ✅ 可作为 Python 包安装（支持 `pip install -e .`）

---

## 📦 安装方法

建议在虚拟环境中使用：

```bash
git clone https://github.com/vesperale1112/FASTQ_Compressor.git
cd FASTQ_Compressor
pip install -e .
```

安装后你将获得一个命令行工具：

```bash
fastq-compress --help
```

---

## 🧪 使用示例

### 🔒 压缩

```bash
fastq-compress --type 1 --method 1 input.fastq output_prefix
```

- `--type 1`：执行压缩
- `--method 1` 或 `2`：选择压缩策略 (1代表标准 k-mer 压缩法， 2代表优化 k-mer 压缩法)
- `input.fastq`：待压缩的原始文件
- `output_prefix`：压缩后生成 `output_prefix.tar.gz`

示例：

```bash
fastq-compress --type 1 --method 2 sample.fastq compressed_data
```

压缩后将生成文件：

```
compressed_data.tar.gz
```

---

### 🔓 解压

```bash
fastq-compress --type 2 compressed_data.tar.gz recovered.fastq
```

也可以显式指定策略：

```bash
fastq-compress --type 2 --method 2 compressed_data.tar.gz recovered.fastq
```

- `--type 2`：解压模式
- `--method`：指定解压策略（可选，默认自动尝试）

---

## 📂 项目结构

```
FASTQ_Compressor/
├── fastq_compressor/
│   ├── __init__.py
│   └── __main__.py         # 主程序入口，包含 CLI 参数解析和压缩/解压逻辑
├── tests/
│   ├── test.fastq          # 测试数据
│   ├── test_result_method1.tar.gz  # 第一种压缩方法压缩文件
│   ├── test_result_method2.tar.gz  # 第二种压缩方法压缩文件
│   ├── recovered_method1.fastq  # 第一种解压缩方法还原文件
│   └── recovered_method2.fastq  # 第二种解压缩方法还原文件
├── pyproject.toml          # 构建配置
├── LICENSE                 # 许可证文件
├── README.md               # 项目说明（本文件）
├── README_en.md            # 项目说明英文版
└── requirements.txt        # 项目依赖
```

---

## 🛠️ 开发依赖

```bash
pip install -r requirements.txt
```

依赖项：

- `biopython`
- `numpy`

---

## 📖 许可证

本项目遵循 MIT License。你可以自由使用、修改、发布和商用，但需保留原始作者信息。

---

## 🙋‍♂️ 作者与贡献者

作者：[@vesperale1112](https://github.com/vesperale1112) [@hao791](https://github.com/hao791) [@fqZzzw](https://github.com/fqZzzw)

欢迎提交 Issue、Pull Request 和 Star ⭐️！
