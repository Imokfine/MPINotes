# MPINotes
### 笔记参考来源
Slides from course  
[MPI tutorial](https://mpitutorial.com/tutorials/)  
[两小时入门MPI与并行计算系列](https://zhuanlan.zhihu.com/p/355652501)  
[Microsoft MPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi)  [(中文版)](https://learn.microsoft.com/zh-cn/message-passing-interface/microsoft-mpi)  

### 常用函数用法列表

##  Point-to-Point Communication 点对点通信  
* 4 Communication modes 四种通信模式
  * MPI_Send / MPI_Recv (Standard 标准)  
  * MPI_Ssend (Synchronous 同步)  
  * MPI_Bsend (Buffered 缓冲) 
  * MPI_Rsend (Ready 就绪)
* Deadlock 死锁
* Buffering 缓冲
* MPI_Sendrecv
* Blocking Message Passing 阻塞式通信  

* Non Blocking Message Passing 非阻塞式通信  
  * MPI_Isend / MPI_Irecv  
  * MPI_Request / MPI_Request_free
  * MPI_Test / MPI_Testall / MPI_Testany / MPI_Testsome
  * MPI_Wait / MPI_Waitall / MPI_Waitany / MPI_Waitsome
  * MPI_Cancel / MPI_Test_cancelled
  
## Collective Communication 集体通信  
* One-to-All
  * MPI_Bcast
  * MPI_Scatter
* All-to-One
  * MPI_Gather
  * MPI_Reduce
* All-to-All
  * MPI_Allgather

## Virtual Topologies 虚拟拓扑  
* Cartesian (grid) topologies 笛卡尔拓扑  
  * MPI_Cart_create
  * MPI_Cart_shift
