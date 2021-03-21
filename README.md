## 基于虚拟麦克风的实时跑调修复  
## 使用方法  
1.安装sox`sudo apt install pavucontrol libsox-dev`  
2.安装soundtouch  
```
git clone https://gitlab.com/soundtouch/soundtouch
cd soundtouch && make && sudo make install
```
3.用qt打开vmic.pro，编译  
4.执行`pactl load-module module-null-sink`加载模块  
5.执行`pavucontrol`，在录音中将源改成空输入  
6.执行qt编译出的程序
### 操作  
执行`./vmic.sh`查看帮助
## 乐谱  
在vmic同目录创建notes.txt文件，启动时将自动载入。也可以在运行时通过`./vmic.sh load [url]`载入  
## 乐谱格式
```
行数
开始时间（秒） 持续时间（秒） 音高（midi id）
开始时间（秒） 持续时间（秒） 音高（midi id）
……
```
例如  
```
3
0 1 47
1 2.5 49
3.5 2.5 53
```
## 乐谱生成  
到[midilib](https://midi.sinriv.com)上传midi文件，进入播放页面，复制下载的链接，将后缀.mid改成.vmic，访问网址即可  
