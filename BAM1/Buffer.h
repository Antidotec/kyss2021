#ifndef BAM_BUFFER_H
#define BAM_BUFFER_H
#include <mutex>
#include <thread>

#include <iostream>
#include <fstream>

using namespace  std;
class BufferConfig{
public:
    BufferConfig();
    BufferConfig(int Buffer_number, int consumerpack_number, int Maxn);
    ~BufferConfig();
public:
    int Buffer_number;
    int write_number;
    int consumerpack_number;
    int Maxn;
};

class Buffer {
public:
    Buffer();
    Buffer(BufferConfig *config,ofstream *fout);
    pair<char *,int> getCap(); //提取一个空的输出内存块
    void initoutput(int id,int pos); //将待输出的内存块id和输出长度放入write和pos数组中
    void output(); //输出已经准备好数据的内存块
    bool is_complete();
    void complete_thread();

public:
    BufferConfig *config;
    ofstream *fout; //输出流
    mutex mtx;
    char **buffer; //输出所用的内存块
    int *write; //待输出的内存块ids数组
    int *pos; //待输出的内存块的数据大小
    int write_bg; //
    int write_ed;
    int *capacity; //空内存块的id数组
    int cap_bg;
    int cap_ed;
};


#endif //BAM_BUFFER_H
