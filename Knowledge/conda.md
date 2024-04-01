# Anaconda

#### What is Anaconda

请参考[Anaconda 个人版_Anaconda 中文网](https://anaconda.org.cn/anaconda/)

Anaconda是一个包管理器，环境管理器，一个Python/R数据科学发行版，以及超过7500个开源包的集合。

当你运行任何程序进行分析，都需要下载相应的包或者依赖，我的理解是统称为工具。比如，分析单细胞数据，你需要用到R语言，需要用到R package的Seurat, ggplot, dplyr等，而且有的时候会想要用到不同的版本。在多组学的免疫组学的分析中，还需要用到一些终端的工具，比如序列比对，流程管理工具snakemake等等。如果我们在自己的电脑上直接下载，那下载路径将会非常的多且杂乱。所以，Anaconda提供一个管理器，在这个管理器中，你可以下载所有你用到的工具（前提是Anaconda有），然后在Anaconda中管理他们的路径等等。

#### How to install Anaconda

官网下载：[Free Download | Anaconda](https://www.anaconda.com/download)

以`Linux`系统为例。

1. 运行安装用的shell脚本

```shell
>>> chmod +x Anaconda3-2024.02-1-Linux-x86_64.sh 
>>> ./Anaconda3-2024.02-1-Linux-x86_64.sh 
```

2. 同意安装许可

```shell
Welcome to Anaconda3 2024.02-1
In order to continue the installation process, please review the license
agreement.
Please, press ENTER to continue:
```

按回车同意，并长按回车读到最后

```shell
Do you accept the licence terms? [yes|no]
>>> yes
```

3. 确认安装路径，如未制定，则按照默认路径安装。默认路径在控制台可以看到

```shell
Anaconda3 will now be installed into this location:
/.../.../.../

- Press ENTER to confirm the location
- Press CTRL-C to abort the installation
- Or specify a different location below

>>> your/path/dir
PREFIX=your/path/dir
```

4. 下载安装

```shell
Unpacking payload ...                                                                     Installing base environment...
Downloading and Extracting Packages:
Downloading and Extracting Packages:
Preparing transaction: done
Executing transaction: / 
    Installed package of scikit-learn can be accelerated using scikit-learn-intelex.
    More details are available here: https://intel.github.io/scikit-learn-intelex
    For example:
        $ conda install scikit-learn-intelex
        $ python -m sklearnex my_application.py
done
installation finished.
```

出现installation finished即表示下载完成

5. 设置是否自启动

输入yes以确保开启shell后自启动

```shell
Do you wish to update your shell profile to automatically initialize conda?
This will activate conda on startup and change the command prompt when activated.
If you'd prefer that conda's base environment not be activated on startup,
   run the following command when conda is activated:

conda config --set auto_activate_base false

You can undo this by running `conda init --reverse $SHELL`? [yes|no]
[no] >>> yes
```

输入no以阻止自启动

更改以使其不自启动

```shell
conda config --set auto_activate_base false
conda config --show
```

#### How to use conda (conda使用指南)

conda常用命令

1. 创建虚拟环境

```shell
conda create -n env_name
conda create -n env_name python=3.11
```

创建好的虚拟环境可以在Anaconda安装目录下的`envs`文件夹中被找到。在不指定python版本的时候，将会下载最新版python。

2. 查看虚拟环境列表

```shell
conda env list
```

3. 激活和退出当前虚拟环境

```shell
conda activate env_name
conda deactivate env_name
```



#### 之后会更新R语言，python，以及相关包的安装！！！

#### 并且会更新Miniconda的使用教程



