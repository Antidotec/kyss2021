#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <htslib/khash.h>
#include <stdint.h>
#include <chrono>
#include "config.h"
#include "BamBlock.h"
#include "Buffer.h"
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

typedef std::chrono::high_resolution_clock Clock;
#ifdef CYCLING
#define TDEF(x_) static unsigned long long int x_##_t0, x_##_t1;
    #define TSTART(x_) x_##_t0 = __rdtsc();
    #define TEND(x_) x_##_t1 = __rdtsc();
    #define TPRINT(x_, str) printf("%-20s \t%.6f\t M cycles\n", str, (double)(x_##_t1 - x_##_t0)/1e6);
#elif defined TIMING
#define TDEF(x_) chrono::high_resolution_clock::time_point x_##_t0, x_##_t1;
#define TSTART(x_) x_##_t0 = Clock::now();
#define TEND(x_) x_##_t1 = Clock::now();
#define TPRINT(x_, str) printf("%-20s \t%.6f\t sec\n", str, chrono::duration_cast<chrono::microseconds>(x_##_t1 - x_##_t0).count()/1e6);
#else
#define TDEF(x_)
#define TSTART(x_)
#define TEND(x_)
#define TPRINT(x_, str)
#endif
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

int n_thread=1;
long long NUM_N[50]={0};
long long NUM_M[50]={0};

void read_pack(BGZF *fp,BamBlock *block){
    pair<bam_block *,int> b;
    b=block->getEmpty();
    int count=0;
    while(read_block(fp,b.first)==0){
        block->inputblock(b.second);
        b=block->getEmpty();
    }
    block->ReadComplete();
}

void write_pack(Buffer *buffer){
    while(!buffer->is_complete()){
        std::this_thread::sleep_for(chrono::milliseconds(10));
        buffer->output();
    }
}

void consumer_pack(BamBlock *block,Buffer *buffer,int id){

    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    pair<bam_block *,int> comp;
    bam_block *un_comp;
    un_comp = (bam_block *)malloc(sizeof(bam_block));
    int pos=0;
    uint8_t* quality;
    uint8_t *seq;
    uint8_t* qname;
    pair<char *, int> thread_buffer;
    thread_buffer=buffer->getCap();
    while (true){
        comp=block->getCompressdata();
        if (comp.second<0) {
            break;
        }
        block_decode_func(comp.first,un_comp);
        block->backempty(comp.second); //这块内存中的压缩块已经信息提取完毕了，可以归还内存
        while (read_bam(un_comp,b,0)>=0) {  //处理提取的信息，格式化输出
            NUM_N[id]++;
            if (b->core.flag&2048) continue;
            NUM_M[id]++;
            seq=bam_get_seq(b);
            long long len = strlen(bam_get_qname(b));
            if (pos>buffer->config->Maxn)
            {
                buffer->initoutput(thread_buffer.second,pos);
                thread_buffer=buffer->getCap();
                pos=0;
            }
            thread_buffer.first[pos++]='@';
            memcpy(thread_buffer.first+pos,bam_get_qname(b),len); //文件名
            pos+=len;
            thread_buffer.first[pos++]='\n';
            if (b->core.flag&16){
                for (int i=0;i<b->core.l_qseq;i++) thread_buffer.first[pos+b->core.l_qseq-1-i] = BaseRever[bam_seqi(seq,i)];
                pos+=b->core.l_qseq;
                thread_buffer.first[pos++]='\n';thread_buffer.first[pos++]='+';thread_buffer.first[pos++]='\n';
                quality=bam_get_qual(b);
                for (int i=0;i<b->core.l_qseq;i++){ thread_buffer.first[pos+i]=(char)quality[b->core.l_qseq-i-1]+33;}
                pos+=b->core.l_qseq;
            }else{
                for (int i=0;i<b->core.l_qseq;i++) thread_buffer.first[pos+i] = Base[bam_seqi(seq,i)];
                pos+=b->core.l_qseq;
                thread_buffer.first[pos++]='\n';thread_buffer.first[pos++]='+';thread_buffer.first[pos++]='\n';
                quality=bam_get_qual(b);
                for (int i=0;i<b->core.l_qseq;i++){ thread_buffer.first[pos+i]=(char)quality[i]+33;}
                pos+=b->core.l_qseq;
            }
            thread_buffer.first[pos++]='\n';
        }
    }
    buffer->initoutput(thread_buffer.second,pos);
    buffer->complete_thread();
}

int main(int argc, char **argv){
    if(argc<4){
        printf("please input using this format: -i inputFile -o outputFile");
        exit(-1);
    }
    char *input_bam = argv[1];
    char *output_fq = argv[3];
    TDEF(fq)
    TSTART(fq)
    printf("Starting Running\n");
    n_thread=1; 
    samFile *sin;
    sam_hdr_t *hdr;
    bam1_t *b;
    ofstream fout;
    fout.open(output_fq);
    if ((sin=sam_open(input_bam, "r"))==NULL){
        printf("Can`t open this file!\n");
        return 0;
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
        return  0;
    }
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    BufferConfig config(50,n_thread,10000000);
    Buffer buffer(&config,&fout);
    BamBlockConfig bamconfig(4000);
    BamBlock block(&bamconfig);
    thread **Bam = new thread *[n_thread+2];
    Bam[0]=new thread(&read_pack,sin->fp.bgzf,&block);
    for (int i=1;i<=n_thread;i++)
        Bam[i]=new thread(&consumer_pack,&block,&buffer,i);
    Bam[n_thread+1]=new thread(&write_pack,&buffer);
    for (int i=0;i<n_thread+2;i++)
        Bam[i]->join();
    long long N=0,M=0;
    for (int i=0;i<50;i++) N+=NUM_N[i];
    for (int i=0;i<50;i++) M+=NUM_M[i];
    printf("total read is %lld\n",N);
    printf("totol process is %lld\n",M);
    sam_close(sin);
    TEND(fq)
    TPRINT(fq,"change time is : ");
    fout.close();
}
