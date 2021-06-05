//
// Created by 赵展 on 2021/3/10.
//
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
#include "BamStatus.h"
#include "Duplicate.h"
#include "Overrepresent.h"
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

int n_thread=1;
long long NUM_N[500]={0};
long long NUM_M[500]={0};

void read_pack(BGZF *fp,BamBlock *block){
    pair<bam_block *,int> b;
    b=block->getEmpty();
    int count=0;
    while(read_block(fp,b.first)==0){
        block->inputblock(b.second);
        //printf("read block is %d\n",++count);
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
    pair<char *,int> thread_buffer;
    thread_buffer=buffer->getCap();
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp=block->getCompressdata();
        //printf("%d is get One compressed data\n",id);
        if (comp.second<0) {
            //printf("%d is Over\n",id);
            break;
        }
        block_decode_func(comp.first,un_comp);
        block->backempty(comp.second);
        while (read_bam(un_comp,b,0)>=0) {
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
            memcpy(thread_buffer.first+pos,bam_get_qname(b),len);
            pos+=len;
            thread_buffer.first[pos++]='\n';
            if (b->core.flag&16){
                for (int i=0;i<b->core.l_qseq;i++) thread_buffer.first[pos+b->core.l_qseq-1-i] = StatusBaseRever[bam_seqi(seq,i)];
                pos+=b->core.l_qseq;
                thread_buffer.first[pos++]='\n';thread_buffer.first[pos++]='+';thread_buffer.first[pos++]='\n';
                quality=bam_get_qual(b);
                for (int i=0;i<b->core.l_qseq;i++){ thread_buffer.first[pos+i]=(char)quality[b->core.l_qseq-i-1]+33;}
                pos+=b->core.l_qseq;
            }else{
                for (int i=0;i<b->core.l_qseq;i++) thread_buffer.first[pos+i] = StatusBase[bam_seqi(seq,i)];
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
void basic_status_pack(BamBlock *block,BamStatus *status){
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        printf("bam1_t is not ready\n");
    }
    pair<bam_block *,int> comp;
    bam_block *un_comp;
    un_comp = (bam_block *)malloc(sizeof(bam_block));
    int pos=0;
    uint8_t* quality;
    uint8_t *seq;
    uint8_t* qname;

    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp=block->getCompressdata();
        //printf("%d is get One compressed data\n",id);
        if (comp.second<0) {
            //printf("%d is Over\n",id);
            break;
        }
        block_decode_func(comp.first,un_comp);
        block->backempty(comp.second);
        while (read_bam(un_comp,b,0)>=0) {
            status->statusbam(b);
        }
    }
}
void basic_duplicate_status_pack(BamBlock *block,BamStatus *status,Duplicate *duplicate){
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        printf("bam1_t is not ready\n");
    }
    pair<bam_block *,int> comp;
    bam_block *un_comp;
    un_comp = (bam_block *)malloc(sizeof(bam_block));
    int pos=0;
    uint8_t* quality;
    uint8_t *seq;
    uint8_t* qname;

    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp=block->getCompressdata();
        //printf("%d is get One compressed data\n",id);
        if (comp.second<0) {
            //printf("%d is Over\n",id);
            break;
        }
        block_decode_func(comp.first,un_comp);
        block->backempty(comp.second);
        while (read_bam(un_comp,b,0)>=0) {
            status->statusbam(b);
            duplicate->statusSeq(b);
        }
    }
}
void basic_duplicate_overrepresent_status_pack(BamBlock *block,BamStatus *status,Duplicate *duplicate,Overrepresent *overrepresent){
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        printf("bam1_t is not ready\n");
    }
    pair<bam_block *,int> comp;
    bam_block *un_comp;
    un_comp = (bam_block *)malloc(sizeof(bam_block));
    int pos=0;
    uint8_t* quality;
    uint8_t *seq;
    uint8_t* qname;
    int over_fg=1;
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp=block->getCompressdata();
        //printf("%d is get One compressed data\n",id);
        if (comp.second<0) {
            //printf("%d is Over\n",id);
            break;
        }
        block_decode_func(comp.first,un_comp);
        block->backempty(comp.second);
        while (read_bam(un_comp,b,0)>=0) {
            status->statusbam(b);
            duplicate->statusSeq(b);
            if (over_fg) overrepresent->insert(b);
            if (over_fg && overrepresent->Pos>=overrepresent->Total) {
                overrepresent->status();
                over_fg=0;
            }
        }
    }
    if (over_fg) overrepresent->status();
}
int main(int argc,char* argv[]){
    TDEF(fq)
    TSTART(fq)
    printf("Starting Running\n");
    n_thread=atoi(argv[1]);
    samFile *sin;
    sam_hdr_t *hdr;
    ofstream fout;
    fout.open("./BamStatus.html");
    // /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam
    // ../data/NC_T_1.sorted.bam
    // ../data/EA_WES_T_1.sorted.bam
    // ../data/SRR_sort.bam
    string path="../data/NC_T_1.sorted.bam";
    if ((sin=sam_open(path.c_str(), "r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
        return  0;
    }
    /*
     *  读取和处理准备
     */
    BamBlockConfig bamconfig(20000);
    BamBlock block(&bamconfig);

    /*
     * 分析准备
     */
    BamStatus **status=new BamStatus*[n_thread];
    Duplicate *duplicate=new Duplicate();
    Overrepresent *overrepresent=new Overrepresent(20000);


    thread **Bam = new thread *[n_thread+1];


    Bam[n_thread]=new thread(&read_pack,sin->fp.bgzf,&block);
    for (int i=0;i<n_thread;i++){
        status[i]=new BamStatus(path);
        if (i==0) {
            Bam[i]=new thread(&basic_duplicate_overrepresent_status_pack,&block,status[i],duplicate,overrepresent);
        } else Bam[i]=new thread(&basic_duplicate_status_pack,&block,status[i],duplicate);
    }
    for (int i=0;i<n_thread+1;i++) Bam[i]->join();
    for (int i=1;i<n_thread;i++) {
        status[0]->add(status[i]);
    }

    status[0]->statusAll();


    status[0]->reportHTML(&fout,duplicate,overrepresent);
    printf("Over Represent is %lf\n",overrepresent->OverrepresentDate);
    sam_close(sin);
    TEND(fq)
    TPRINT(fq,"change time is : ");
}
