//
// Created by 赵展 on 2021/3/10.
//

#ifndef BAM_BAMBLOCK_H
#define BAM_BAMBLOCK_H
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <thread>
#include <htslib/khash.h>
#include "header.h"
#include <stdint.h>
#include "config.h"
#include <cstdlib>
#include <mutex>
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

typedef struct
{
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
KHASH_MAP_INIT_INT64(cache, cache_t)

typedef struct bam_block bam_block;
struct bam_block
{
    unsigned int errcode;
    unsigned char data[BGZF_MAX_BLOCK_SIZE]; //0x1000
    unsigned int length;
    unsigned int pos;
    int64_t block_address;
};
struct bgzf_cache_t
{
    khash_t(cache) * h;
    khint_t last_pos;
};
int sam_realloc_bam_data(bam1_t *b, size_t desired);
inline int realloc_bam_data(bam1_t *b, size_t desired);
inline int possibly_expand_bam_data(bam1_t *b, size_t bytes);
int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning); // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                      hts_pos_t *rlen, hts_pos_t *qlen);
void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host);
inline int unpackInt16(const uint8_t *buffer);
int load_block_from_cache(BGZF *fp, int64_t block_address);
int check_header(const uint8_t *header);
const char *bgzf_zerr(int errnum, z_stream *zs);
//int bgzf_uncompress(uint8_t *dst, size_t *dlen,const uint8_t *src, size_t slen,uint32_t expected_crc);
int bgzf_uncompress(uint8_t *dst, unsigned int *dlen, const uint8_t *src, unsigned int slen, uint32_t expected_crc);
int fixup_missing_qname_nul(bam1_t *b);
int block_decode_func(struct bam_block *comp, struct bam_block *un_comp);
int read_block(BGZF *fp, struct bam_block *j);
int Rabbit_bgzf_read(struct bam_block *fq, void *data, unsigned int length);
int read_bam(struct bam_block *fq, bam1_t *b, int is_be);
using namespace std;
class BamBlockConfig
{
public:
    BamBlockConfig();
    BamBlockConfig(int Buffer_number);

public:
    int Buffer_number;
    int write_number;
    int complete;
};

class BamBlock
{
public:
    BamBlock();
    BamBlock(BamBlockConfig *config);
    pair<bam_block *, int> getEmpty(); //提取一个空的内存块
    void inputblock(int id); // 导入未解压的数据
    pair<bam_block *, int> getCompressdata(); //提取一个填入了压缩数据的内存块
    void backempty(int id); //归还一个已经解析完的内存块
    bool isComplete(); //是否读取结束了
    void ReadComplete(); //读取结束

public:
    BamBlockConfig *config;
    mutex mtx_read;
    mutex mtx_compress;
    bam_block **buffer; //内存块
    int *compress; //压缩数据内存块id数组 (循环数组)
    int compress_bg; //头指针
    int compress_ed;//尾指针
    int *read; //空内存块的id数组,其它同上
    int read_bg;
    int read_ed;
};

#endif //BAM_BAMBLOCK_H
