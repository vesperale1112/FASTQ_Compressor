# FASTQ_Compressor

**FASTQ_Compressor** is a command-line tool designed for efficient compression and decompression of FASTQ files. It supports two compression strategies and outputs a `.tar.gz` archive containing sequence, quality, and positional data, which can be fully restored to a valid FASTQ file.

---

## 🚀 Features

- ✅ Two k-mer compression strategies (standard and optimized) and a binary-like compression strategy supported
- ✅ Uses Biopython for FASTQ parsing
- ✅ Easy-to-use CLI tool: `fastq-compress`
- ✅ Automatically detects or manually selects decompression strategy
- ✅ Can be installed as a Python package (`pip install -e .`)

---

## 📦 Installation

We recommend installing in a virtual environment:

```bash
git clone https://github.com/vesperale1112/FASTQ_Compressor.git
cd FASTQ_Compressor
pip install -e .
```

After installation, the following CLI tool becomes available:

```bash
fastq-compress --help
```

---

## 🧪 Usage Examples

### 🔒 Compression

```bash
fastq-compress --type 1 --method 1 input.fastq output_prefix
```

- `--type 1`: perform compression
- `--method 1` or `2`: choose compression strategy (1 represents the standard k-mer compression strategy, 2 represents the optimized k-mer compression strategy)
- `input.fastq`: original FASTQ file
- `output_prefix`: prefix of compressed output `.tar.gz`

Example:

```bash
fastq-compress --type 1 --method 2 sample.fastq compressed_data
```

Output:

```
compressed_data.tar.gz
```

---

### 🔓 Decompression

```bash
fastq-compress --type 2 compressed_data.tar.gz recovered.fastq
```

You can also explicitly specify the decompression strategy:

```bash
fastq-compress --type 2 --method 2 compressed_data.tar.gz recovered.fastq
```

- `--type 2`: perform decompression
- `--method`: optional; if omitted, auto-detection will be attempted

---

## 📂 Project Structure

```
FASTQ_Compressor/
├── fastq_compressor/
│   ├── __init__.py
│   └── __main__.py         # Main script entry, CLI logic and core functionality
├── tests/
│   ├── test.fastq          # Test Data
│   ├── test_result_method1.tar.gz  # Method 1 Compressing Example
│   ├── test_result_method2.tar.gz  # Method 2 Compressing Example
│   ├── recovered_method1.fastq  # Method 1 Decompressing Example
│   └── recovered_method2.fastq  # Method 2 Decompressing Example
├── pyproject.toml          # Build configuration
├── LICENSE                 # License
├── README.md               # Project documentation (Chinese Version)
├── README_en.md            # Project documentation (this file)
└── requirements.txt        # Dependency list
```

---

## 🛠️ Development Dependencies

```bash
pip install -r requirements.txt
```

Dependencies:

- `biopython`
- `numpy`

---

## 📖 License

This project is licensed under the MIT License. You are free to use, modify, distribute, and commercialize it, as long as you retain the original author information.

---

## 🙋‍♂️ Author & Contributors

Author: [@vesperale1112](https://github.com/vesperale1112) [@hao791](https://github.com/hao791) [@fqZzzw](https://github.com/fqZzzw)

Contributions, Issues, and Stars ⭐️ are welcome!
