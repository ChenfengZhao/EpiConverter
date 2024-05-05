# EpiConverter

<!-- <p align="center">
    <img width="30%" src="./EpiConverter_Logo.png" alt="EpiConverter"><br>
</p> -->

## Overview

EpiConverter is built to properly convert mzXML files to a usable input structure for downstream analysis. The converted data can be further searched using EpiProfile (version 2.1), as the generated “ms1” and “ms2” files can be examined to resolve the comprehensive histone PTM information, including the isobaric histone peptides, low-abundant PTMs, different histone mutations, different derivatization strategies, and histone sequences from various organisms. This tool can be installed on all kinds of platforms (e.g., Linux, macOS, Windows, etc).

## Requirements

- Recommended OS: macOS (>= 10.13), Linux (e.g. Ubuntu >= 18.04), or Windows (>= 10)
- Python3 (3.7 or higher is supported)
- Pip3
- Python dependencies: numpy, scipy, matchms, pyyaml
- Conda (optional): Miniconda or Anaconda

## Installation Guide

In order to use EpiConverter, you can install the dependent libraries either via package manager (Pip3) or by creating a new virtual environment with Conda

### Install via package

```
pip3 install scipy numpy matchms[chemistry] pyyaml
```

### Install via Conda

```
conda create --name matchms python=3.8
conda activate matchms
conda install --channel bioconda --channel conda-forge matchms
```

## Usage Example
1. Download the source code of EpiConverter. You can manully download the zip file and unzip it, or you can use the following command to directly download it.
```
git clone https://github.com/ChenfengZhao/EpiConverter.git
```
The github repository is composed of the following parts:
- **EpiConverter_mp.py** contains the python code of EpiConverter.
- **config_EpiConverter.ini** is to configure the internal parameters of the tool.

Here are the explaination on the parameters in config_EpiConverter.ini:
- **input_path**: input file folder location
- **output_path**: output file folder location
- **ion_mobility**: (under development) empolying ion mobility
- **core_num**: number of cores used for the conversion
- **filter_dataset**: precheck the data file
- **MS1 ion_injection_time**: ion accumulation time for MS1 scan
- **MS1 instrument_type**: MS analyzer for MS1 scan
- **MS2 ion_injection_time**: ion accumulation time for MS2 scan
- **MS2 activation_type**: fragmentation type
- **MS2 instrument_type**: MS analyzer for MS2 scan
- **acquisition_method**: type of MS analysis
  
1. Prepare your MS files following the format of the example. Put these files under the path defined in the **input_path** of config_EpiConverter.ini. EpiConverter automatically searches and processes all the MS files in right format in batches on multi-core CPUs. Feel free to skip this step if you are just intented to process the example data.

2. Execute EpiConverter using the following code:

If the dependencies are installed via package manager
```
cd <the path of EpiConverter_mp.py>
python3 EpiConverter_mp.py
```

If the dependencies are installed by Conda
```
conda activate matchms
cd <the path of EpiConverter_mp.py>
python3 EpiConverter_mp.py
```

3. The results will be generated under the path defined in the **output_path** of config_EpiConverter.ini (e.g., ./result/ in this example).


## License
[Apache_2.0_license]: http://www.apache.org/licenses/LICENSE-2.0

The source code of this project is released under the [Apache 2.0 License][Apache_2.0_license].

## Citation
If you think EpiConverter is helpful for your research, please cite the following paper:

<!-- Xie, Y; De Luna Vitorino, F.N.; Chen, Y; Lempiäinen, J. K.; Zhao, C.; Steinbock, R. T.; Liu, X.; Lin, Z.; Zahn, E.; Garcia, A. L.; Weitzman, M. D.; Garcia, B. A., SWAMNA: a comprehensive platform for analysis of nucleic acid modifications. *Chemical Communications* **2023** -->

<!-- [Add the paper reference in the format of above or keep it ANONYMOUS if required.] -->
[Removed to preserve anonymity]

<!-- [ChemComm version](https://doi.org/10.1039/D3CC04402E)

[ChemRxiv version](https://chemrxiv.org/engage/chemrxiv/article-details/64f6a89079853bbd781e9eb7) -->
