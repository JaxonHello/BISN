# mac git 指南

## 下载安装

访问官方网站[Git - Downloading Package (git-scm.com)](https://git-scm.com/download/mac)了解更多细节

one can visit [Git - Downloading Package (git-scm.com)](https://git-scm.com/download/mac) for more details.

### 安装命令

```bash
brew install git

# 验证是否下载成功
git --version
> git version 2.39.3 (Apple Git-145)
```

出现以上结果则代表下载成功

注：在此之前请确保下载好`homebrew`组件

### 下载问题及相应对策

```bash
brew install git
> Error: Git is unavailable

git --version
> xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun
```

如果出现上述情况，说明你的Xcode开发工具的命令行工具可能未被正确安装

解决方法：安装或者更新Xcode开发工具

```bash
xcode-select --install
```

等待数分钟即可安装成功，安装成功后验证

## git使用

###  案例：上传文件夹dir到Gitee

1. 在github上新建一个仓库 （GitHub主页）

以下操作均在`terminal终端`中进行

2. Git全局设置

```bash
git config --global user.name "JaxonHe"
git config --global user.email "jaxonhe1021@outlook.com"

# 查看全局设置
git config --list
# 初始配置如下
> credential.helper=osxkeychain
> init.defaultbranch=main
> user.name=JaxonHe
> user.email=jaxonhe1021@outlook.com
```

3. 进入想要上传的目标文件夹

```bash
cd Users/your_name/destination_dir
```

4. Git设置

```bash
# git初始化
git init
# 连接远程仓库
# 连接名称是 https://github.com/username/目标地址仓库名称.git
git remote add origin https://github.com/username/your_repository.git
```

在这一步可能出现以下错误：

```bash
git remote add origin https://github.com/JaxonHello/JaxonHello.github.io.git
> error: remote origin already exists.
```

这里表示Git已经存在一个名为origin的远程仓库配置，如果想要更改远程仓库的URL，请使用以下命令：

```bash
git remote set-url origin https://github.com/JaxonHello/JaxonHello.github.io.git
```

之后

```bash
# 添加文件夹中的所有文件到Git追踪
git add .

# 提交更改
git commit -m 'commit_description'
# 这一步会出现要更改的文件名
> [main (root-commit) 884cf29] git_test_commit
> 1 file changed, 1 insertion(+)
> create mode 100644 test_git.txt

# 推动到远程仓库的main分支
git push -u origin main

# 有可能在这个过程中需要输入github用户名和密码
```

会出现以下问题

```bash
Username for 'https://github.com': JaxonHello
Password for 'https://JaxonHello@github.com': 
remote: Support for password authentication was removed on August 13, 2021.
remote: Please see https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-https-urls for information on currently recommended modes of authentication.
fatal: Authentication failed for 'https://github.com/JaxonHello/JaxonHello.github.io.git/'
```

原因：GitHub移除了对密码身份验证的支持，而推荐使用个人访问令牌进行验证。

5. 生成个人访问令牌

`github` > `Setting` > `Developer settings` > `Personal access tokens`

注：该tokens至少要具有repo权限

6. 重新配置git

```bash
# 先移除原有的远程仓库配置
git remote remove origin  
git remote add origin https://JaxonHello:<YOUR_PERSONAL_ACCESS_TOKEN>@github.com/JaxonHello/JaxonHello.github.io.git
git push -u origin main
```

然后就可以在github上看到上传成功的文件