

#include "BamStatus.h"

BamStatus::BamStatus(string filename){
    this->filename=filename;
    this->ContentLen=1001;
    this->QulityLen=50;
    this->LengthSequence=new int[MAXLEN];
    memset(this->LengthSequence,0x00,MAXLEN*sizeof(int));
    this->QualitySequence=new int[ContentLen];
    memset(this->QualitySequence,0x00,ContentLen*sizeof(int));
    this->NumberList=new int*[MAXLEN];
    this->Qualitylist=new int*[MAXLEN];
    this->QualityPositonList=new int*[MAXLEN];
    this->Content=new int*[ContentLen];
    this->DoubleContent=new double*[ContentLen];
    for (int i=0;i<MAXLEN;i++){
        this->NumberList[i]=new int[8];
        memset(this->NumberList[i],0x00,8*sizeof(int));
        this->Qualitylist[i]=new int[8];
        memset(this->Qualitylist[i],0x00,8*sizeof(int));
        this->QualityPositonList[i]=new int[QulityLen];
        memset(this->QualityPositonList[i],0x00,QulityLen*sizeof(int));
    }
    for (int i=0;i<ContentLen;i++){
        this->Content[i]=new int[8];
        memset(this->Content[i],0x00,8*sizeof(int));
        this->DoubleContent[i]=new double[8];
        memset(this->DoubleContent[i],0x00,8*sizeof(double));
    }
    this->KmerBit=1<<(2*this->KmerBase);
    this->Kmer=new int*[this->KmerBit];
    for (int i=0;i<this->KmerBit;i++){
        this->Kmer[i]=new int[MAXLEN];
        memset(this->Kmer[i],0x00,MAXLEN*sizeof(int));
    }
    this->ChooseKmerNum=6;
    this->ChooseKmerPos=new int[this->ChooseKmerNum];
    this->ChooseKmerKey=new int[this->ChooseKmerNum];
    memset(this->ChooseKmerKey,0xFFFF,this->ChooseKmerNum*sizeof(int));
    this->total_number=0;
    this->total_aligen_number=0;
    this->max_len=0;
    this->min_len=-1;
}

BamStatus::~BamStatus(){
    for (int i=0;i<MAXLEN;i++){
        delete this->NumberList[i];
        delete this->Qualitylist[i];
    }
    for (int i=0;i<ContentLen;i++){
        delete this->Content[i];
        delete this->DoubleContent[i];
    }
    delete this->NumberList;
    delete this->Qualitylist;
}

void BamStatus::statusbam(bam1_t *b) {
    total_number++;
    max_len=max_len>b->core.l_qseq?max_len:b->core.l_qseq;
    min_len=min_len==-1?b->core.l_qseq:min(min_len,b->core.l_qseq);
    LengthSequence[b->core.l_qseq]++;
    this->Chr.insert(b->core.tid);
    if (b->core.flag&2048) return;
    total_aligen_number++;


    /*
     *  序列内部信息
     */
    uint8_t* quality;
    uint8_t *seq;
    uint8_t* qname;
    seq=bam_get_seq(b);
    int number[8];memset(number,0x00,8*sizeof(int));
//    printf("this is number before\n");
//    for (int i=0;i<8;i++){
//        printf("%d %d\n",i,number[i]);
//    }
    int total_qual=0;
    long long len = strlen(bam_get_qname(b));
    if (b->core.flag&16){
        int kmer=0;
        int last_n=-1;
        for (int i=0;i<b->core.l_qseq;i++) {
            NumberList[b->core.l_qseq-i-1][StatusBaseRever[bam_seqi(seq,i)]&0x07]++;
            number[StatusBaseRever[bam_seqi(seq,i)]&0x07]++;
            kmer<<=2;
            kmer+=DupvalAGCT[StatusBaseRever[bam_seqi(seq,i)]&0x07];
            kmer&=0x03FF;
            if ((DupvalAGCT[StatusBaseRever[bam_seqi(seq,i)]&0x07])==-1) {
                kmer=0;
                last_n=i;
            }
            if (i-last_n>=kmer) Kmer[kmer][i]++;
        }
        quality=bam_get_qual(b);
        for (int i=0;i<b->core.l_qseq;i++)  {
            Qualitylist[b->core.l_qseq-i-1][StatusBaseRever[bam_seqi(seq,i)]&0x07]+=quality[b->core.l_qseq-i-1];
            QualityPositonList[b->core.l_qseq-i-1][quality[b->core.l_qseq-i-1]]++;
            total_qual+=quality[b->core.l_qseq-i-1];
        }
    }else{
        int kmer=0;
        int last_n=-1;
        for (int i=0;i<b->core.l_qseq;i++)  {
            NumberList[i][StatusBase[bam_seqi(seq,i)]&0x07]++;
            number[StatusBaseRever[bam_seqi(seq,i)]&0x07]++;
            kmer<<=2;kmer+=DupvalAGCT[StatusBaseRever[bam_seqi(seq,i)]&0x07];kmer&=0x03FF;
            if ((DupvalAGCT[StatusBaseRever[bam_seqi(seq,i)]&0x07])==-1) {
                kmer=0;
                last_n=i;
            }
            if (i-last_n>=kmer) Kmer[kmer][i]++;
        }
        quality=bam_get_qual(b);
        for (int i=0;i<b->core.l_qseq;i++)  {
            Qualitylist[i][StatusBaseRever[bam_seqi(seq,i)]&0x07]+=quality[i];
            QualityPositonList[i][quality[i]]++;
            total_qual+=quality[i];
        }
    }
    int total=0;
    for (int i=0;i<8;i++) total+=number[i];
    Content[((ContentLen-1)*number['A'&0x07])/total]['A'&0x07]++;
    Content[((ContentLen-1)*number['G'&0x07])/total]['G'&0x07]++;
    Content[((ContentLen-1)*number['C'&0x07])/total]['C'&0x07]++;
    Content[((ContentLen-1)*number['T'&0x07])/total]['T'&0x07]++;
    Content[((ContentLen-1)*number['N'&0x07])/total]['N'&0x07]++;
    Content[((ContentLen-1)*(number['G'&0x07]+number['C'&0x07]))/total][('C'+'G')&0x07]++;


    /*
     * 序列本身信息
     */
    QualitySequence[total_qual/b->core.l_qseq]++;

    //printf("%d is %d\n",b->core.l_qseq,LengthList[b->core.l_qseq]);
}

void BamStatus::add(BamStatus *b){
    for (set<int>::iterator i=b->Chr.begin();i!=b->Chr.end();i++){
        this->Chr.insert(*i);
    }
    for (int i=0;i<MAXLEN;i++){
        for (int j=0;j<8;j++){
            this->NumberList[i][j]+=b->NumberList[i][j];
            this->Qualitylist[i][j]+=b->Qualitylist[i][j];
        }
        for (int j=0;j<QulityLen;j++){
            this->QualityPositonList[i][j]+=b->QualityPositonList[i][j];
        }
    }
    for (int i=0;i<ContentLen;i++){e
        for (int j=0;j<8;j++){
            this->Content[i][j]+=b->Content[i][j];
        }
    }
    for (int i=0;i<MAXLEN;i++){
        this->LengthSequence[i]+=b->LengthSequence[i];
    }
    for (int i=0;i<KmerBit;i++){
        for (int j=0;j<MAXLEN;j++){
            this->Kmer[i][j]+=b->Kmer[i][j];
        }
    }
    this->total_number+=b->total_number;
    this->total_aligen_number+=b->total_aligen_number;
    this->max_len=this->max_len>b->max_len?this->max_len:b->max_len;
    this->min_len=this->min_len==1?b->min_len:min(this->min_len,b->min_len);
}
void BamStatus::statusAll() {
    /*
     * 这次先采用MAX-MIN的方法进行测试
     */
    for (int k=0;k<KmerBit;k++){
        int _min=-1,_max=-1;
        for (int i=0;i<MAXLEN;i++){
            _min=_min==-1?Kmer[k][i]:min(_min,Kmer[k][i]);
            _max=_max==-1?Kmer[k][i]:max(_max,Kmer[k][i]);
        }
        for (int i=0;i<ChooseKmerNum;i++){
            if (_max-_min>ChooseKmerKey[i]){
                for (int j=ChooseKmerNum-1;j>i;j--){
                    ChooseKmerPos[j]=ChooseKmerPos[j-1];
                    ChooseKmerKey[j]=ChooseKmerKey[j-1];
                }
                ChooseKmerKey[i]=_max-_min;
                ChooseKmerPos[i]=k;
                break;
            }
        }

    }
}
void BamStatus::contentstatus() {
    for (int i=0;i<101;i++){
        double total = 0.0;
        for (int j=0;j<8;j++) total+=Content[j][i];
    }
}
void BamStatus::print(){
    /*
     *  每个位置上的AGCTN的比例
     */
    {
        for (int i=0;i<max_len;i++){
            int num=0;
            for (int j=0;j<8;j++) num+=NumberList[i][j];
            printf("the Number postion %d\n",i);
            printf("A is %lf\n",NumberList[i]['A'&0x07]/(num+0.0));
            printf("G is %lf\n",NumberList[i]['G'&0x07]/(num+0.0));
            printf("C is %lf\n",NumberList[i]['C'&0x07]/(num+0.0));
            printf("T is %lf\n",NumberList[i]['T'&0x07]/(num+0.0));
            printf("N is %lf\n",NumberList[i]['N'&0x07]/(num+0.0));
        }
    }


    /*
     * AGCTN在整个序列上的比例
     */
    {
        int mod=1;
        for (int i=0;i<ContentLen;i++){
            printf("the Content postion %d\n",i);
            printf("A Content is %d\n",Content[i]['A'&0x07]*mod);
            printf("G Content is %d\n",Content[i]['G'&0x07]*mod);
            printf("C Content is %d\n",Content[i]['C'&0x07]*mod);
            printf("T Content is %d\n",Content[i]['T'&0x07]*mod);
            printf("N Content is %d\n",Content[i]['N'&0x07]*mod);
            printf("GC Content is %d\n",Content[i][('G'+'C')&0x07]*mod);
        }
    }


    /*
     * 序列长度的数量 或 比例
     */

    {
        int total_length_number=0;
        for (int i=0;i<this->max_len;i++) total_length_number+=LengthSequence[i];
        printf("total_length_number is %d\n",total_length_number);
        for (int i=0;i<this->max_len;i++){
            printf("%d length content is %d\n",i,LengthSequence[i]);
        }
    }

    {
        for (int i=0;i<41;i++){
            printf("%d  Mean Qulity number is  %d\n",i,QualitySequence[i]);
        }
    }

}
string HTMLHeader(){
    return string("<!DOCTYPE html>\n"
                       "<html lang=\"en\">\n"
                       "<head>\n"
                       "    <meta charset=\"UTF-8\">\n"
                       "    <title>Title</title>\n"
                       "    <script src=\"https://cdn.jsdelivr.net/npm/echarts@5.0.2/dist/echarts.min.js\"></script>\n"
                       "</head>\n"
                       "<body>\n");
}
string HTMLCss(){
    return "<style type=\"text/css\">\n"
           " @media screen {\n"
           "  div.summary {\n"
           "    width: 18em;\n"
           "    position:fixed;\n"
           "    top: 3em;\n"
           "    margin:1em 0 0 1em;\n"
           "  }\n"
           "  \n"
           "  div.main {\n"
           "    display:block;\n"
           "    position:absolute;\n"
           "    overflow:auto;\n"
           "    height:auto;\n"
           "    width:auto;\n"
           "    top:4.5em;\n"
           "    bottom:2.3em;\n"
           "    left:18em;\n"
           "    right:0;\n"
           "    border-left: 1px solid #CCC;\n"
           "    padding:0 0 0 1em;\n"
           "    background-color: white;\n"
           "    z-index:1;\n"
           "  }\n"
           "  \n"
           "  div.header {\n"
           "    background-color: #EEE;\n"
           "    border:0;\n"
           "    margin:0;\n"
           "    padding: 0.5em;\n"
           "    font-size: 200%;\n"
           "    font-weight: bold;\n"
           "    position:fixed;\n"
           "    width:100%;\n"
           "    top:0;\n"
           "    left:0;\n"
           "    z-index:2;\n"
           "  }\n"
           "\n"
           "  div.footer {\n"
           "    background-color: #EEE;\n"
           "    border:0;\n"
           "    margin:0;\n"
           "\tpadding:0.5em;\n"
           "    height: 1.3em;\n"
           "\toverflow:hidden;\n"
           "    font-size: 100%;\n"
           "    font-weight: bold;\n"
           "    position:fixed;\n"
           "    bottom:0;\n"
           "    width:100%;\n"
           "    z-index:2;\n"
           "  }\n"
           "  \n"
           "  img.indented {\n"
           "    margin-left: 3em;\n"
           "  }\n"
           " }\n"
           " \n"
           " @media print {\n"
           "\timg {\n"
           "\t\tmax-width:100% !important;\n"
           "\t\tpage-break-inside: avoid;\n"
           "\t}\n"
           "\th2, h3 {\n"
           "\t\tpage-break-after: avoid;\n"
           "\t}\n"
           "\tdiv.header {\n"
           "      background-color: #FFF;\n"
           "    }\n"
           "\t\n"
           " }\n"
           " \n"
           " body {    \n"
           "  font-family: sans-serif;   \n"
           "  color: #000;   \n"
           "  background-color: #FFF;\n"
           "  border: 0;\n"
           "  margin: 0;\n"
           "  padding: 0;\n"
           "  }\n"
           "  \n"
           "  div.header {\n"
           "  border:0;\n"
           "  margin:0;\n"
           "  padding: 0.5em;\n"
           "  font-size: 200%;\n"
           "  font-weight: bold;\n"
           "  width:100%;\n"
           "  }    \n"
           "  \n"
           "  #header_title {\n"
           "  display:inline-block;\n"
           "  float:left;\n"
           "  clear:left;\n"
           "  }\n"
           "  #header_filename {\n"
           "  display:inline-block;\n"
           "  float:right;\n"
           "  clear:right;\n"
           "  font-size: 50%;\n"
           "  margin-right:2em;\n"
           "  text-align: right;\n"
           "  }\n"
           "\n"
           "  div.header h3 {\n"
           "  font-size: 50%;\n"
           "  margin-bottom: 0;\n"
           "  }\n"
           "  \n"
           "  div.summary ul {\n"
           "  padding-left:0;\n"
           "  list-style-type:none;\n"
           "  }\n"
           "  \n"
           "  div.summary ul li img {\n"
           "  margin-bottom:-0.5em;\n"
           "  margin-top:0.5em;\n"
           "  }\n"
           "\t  \n"
           "  div.main {\n"
           "  background-color: white;\n"
           "  }\n"
           "      \n"
           "  div.module {\n"
           "  padding-bottom:1.5em;\n"
           "  padding-top:1.5em;\n"
           "  }\n"
           "\t  \n"
           "  div.footer {\n"
           "  background-color: #EEE;\n"
           "  border:0;\n"
           "  margin:0;\n"
           "  padding: 0.5em;\n"
           "  font-size: 100%;\n"
           "  font-weight: bold;\n"
           "  width:100%;\n"
           "  }\n"
           "\n"
           "\n"
           "  a {\n"
           "  color: #000080;\n"
           "  }\n"
           "\n"
           "  a:hover {\n"
           "  color: #800000;\n"
           "  }\n"
           "      \n"
           "  h2 {\n"
           "  color: #800000;\n"
           "  padding-bottom: 0;\n"
           "  margin-bottom: 0;\n"
           "  clear:left;\n"
           "  }\n"
           "\n"
           "  table { \n"
           "  margin-left: 3em;\n"
           "  text-align: center;\n"
           "  }\n"
           "  \n"
           "  th { \n"
           "  text-align: center;\n"
           "  background-color: #000080;\n"
           "  color: #FFF;\n"
           "  padding: 0.4em;\n"
           "  }      \n"
           "  \n"
           "  td { \n"
           "  font-family: monospace; \n"
           "  text-align: left;\n"
           "  background-color: #EEE;\n"
           "  color: #000;\n"
           "  padding: 0.4em;\n"
           "  }\n"
           "\n"
           "  img {\n"
           "  padding-top: 0;\n"
           "  margin-top: 0;\n"
           "  border-top: 0;\n"
           "  }\n"
           "\n"
           "  \n"
           "  p {\n"
           "  padding-top: 0;\n"
           "  margin-top: 0;\n"
           "  }\n"
           "</style>\n";
}
string insertDiv(string id){
    return "<div id=\""+id+"\" style=\"width: 800px;height:600px;\"></div>\n";
}
string insertTooltip(){
    return "tooltip: {\n"
           "        trigger: 'axis',\n"
           "        axisPointer: {\n"
           "            type: 'cross',\n"
           "            label: {\n"
           "                backgroundColor: '#6a7985'\n"
           "            }\n"
           "        }\n"
           "    },\n";
}
string insertChart(string id){
    return "\nvar "+id+"Chart = echarts.init(document.getElementById(\'"+id+"\'));\n";
}
string insertChartOption(string id){
    return id+"Chart.setOption("+id+"Option);\n";
}
string insertOptionBegin(string id){
    return "\nvar "+id+"Option = {\n";
}
string insertOptionEnd(){
    return "}\n";
}
string insertTitle(string text){
    return "title: {\n"
           "            text: \'"+text+"\',\n"
           "            left:\'center\',\n"
           "        },\n";
}
string insertxAxis(int len,int interval=1){
    string out("xAxis: {\n"
               "            type: 'category',\n"
               "            data: [");
    for (int i=0;i<len;i+=interval){
        out.append(to_string(i)+',');
    }
    out.append("]\n},\n");
    return out;
}
string insertxAxis(string name,int len,int interval=1){
    string out("xAxis: {\n"
               "            type: 'category',\n"
               "            name: \'"+name+"\',\n"
               "            nameLocation:'center',\n"
               "nameTextStyle: {\n"
               "                lineHeight: 50,\n"
               "                fontSize: 13,\n"
               "                fontFamily: \"monospace\",\n"
               "                fontWeight: \"bold\"\n"
               "            },\n"
               "            data: [");
    for (int i=0;i<len;i+=interval){
        out.append(to_string(i)+',');
    }
    out.append("]\n},\n");
    return out;
}
string insertyAxis(string type,string _min,string _max){
    return "yAxis: {\n"
           "            type: \'"+type+"\',\n"
           "            min:"+_min+",\n"
           "            max:"+_max+"\n"
           "        },\n";
}
string insertyAxis(string type){
    return "yAxis: {\n"
           "            type: \'"+type+"\',\n"
           "            axisLabel: {\n"
           "                formatter: function (value) {\n"
           "                    if (value < 100000) return value;\n"
           "                    var len=0;\n"
           "                    var last=value;\n"
           "                    while (value>9){\n"
           "                        value=value/10;\n"
           "                        len++;\n"
           "                    }\n"
           "                    return last/Math.pow(10,len)+'E+'+len;\n"
           "              },\n"
           "           },\n"
           "        },\n";
}
string insertSeriesBegin(){
    return "series: [\n";
}
string insertSeriesEnd(){
    return "],\n";
}
string insertSeriesData(string type,int *data,int len,int interval=1){
    string out("{\n"
               "            type: \'"+type+"\',\n"
               "            data: [\n");
    for(int i=0;i<len;i+=interval){
        out.append(to_string(data[i])+',');
    }
    out.append("\n]},\n");
    return out;
}
string insertSeriesData(string type,string name,int *data,int len,int interval=1){
    string out("{\n"
               "            type: \'"+type+"\',\n"
               "            name: \'"+name+"\',\n"
               "            data: [\n             ");
    for(int i=0;i<len;i+=interval){
        out.append(to_string(data[i])+',');
    }
    out.append("\n]},\n");
    return out;
}
string insertSeriesData(string type,double *data,int len,int interval=1){
    string out("{\n"
               "            type: \'"+type+"\',\n"
               "            data: [\n            ");
    for(int i=0;i<len;i+=interval){
        out.append(to_string(data[i])+',');
    }
    out.append("\n]},\n");
    return out;
}
string insertSeriesData(string type,string name,double *data,int len,int interval=1){
    string out("{\n"
               "            type: \'"+type+"\',\n"
               "            name: \'"+name+"\',\n"
               "            data: [\n");
    for(int i=0;i<len;i+=interval){
        out.append(to_string(data[i])+',');
    }
    out.append("\n]},\n");
    return out;
}
string insertSeriesMultiDataBegin(string type,string name){
    return "{\n"
           "            type: \'"+type+"\',\n"
           "            name:\'"+name+"\'"
           "            data: [\n";
}
string insertSeriesMultiDataBegin(string type){
    return "{\n"
           "            type: \'"+type+"\',\n"
           "            data: [\n";
}
string insertSeriesMultiDataEnd(){
    return "],\nitemStyle: {\n"
           "        color: \"#d7ab82\"\n"
           "    },\n},\n";
}
string insertSeriesOneData(int *data,int len){
    string out("[");
    for (int i=0;i<len;i++){
        out.append(to_string(data[i])+',');
    }
    out.append("],\n");
    return out;
}
string insertDataZoom(){
    return "dataZoom: [\n"
           "        {\n"
           "            id: 'dataZoomX',\n"
           "            type: 'slider',\n"
           "            xAxisIndex: [0],\n"
           "            filterMode: 'filter', // 设定为 'filter' 从而 X 的窗口变化会影响 Y 的范围。\n"
           "            start: 0,\n"
           "            end: 100\n"
           "        },\n"
           "        {\n"
           "            id: 'dataZoomY',\n"
           "            type: 'slider',\n"
           "            yAxisIndex: [0],\n"
           "            filterMode: 'empty',\n"
           "            start: 0,\n"
           "            end: 100\n"
           "        }\n"
           "    ],\n";
}
string insertLegend(string data){
    return "legend: {\n"
           "        data: ["+data+"],\n"
           "        left:\'85%\',\n"
           "    },\n";
}
string insertTableBegin(){
    return "<table>\n";
}
string insertTableEnd(){
    return "</table>\n";
}
string insertTableTitle(string str1,string str2){
    return "<thead>\n"
           "    <tr>\n"
           "        <th>"+str1+"</th>\n"
           "        <th>"+str2+"</th>\n"
           "    </tr>\n"
           "</thead>\n";
}
string insertTableTbobyBegin(){
    return "<tboby>\n";
}
string insertTableTbobyEnd(){
    return "</tboby>\n";
}
string insertTableTr(string str1,string str2){
    return "    <tr>\n"
           "        <td>"+str1+"</td>\n"
           "        <td>"+str2+"</td>\n"
           "    </tr>\n";
}
void BamStatus::reportHTML(ofstream *fout){
    string outhtml;

    int *tmp_int=new int[MAXLEN];
    double *tmp_double=new double[MAXLEN];

    string QualityScore("QualityScore");
    string LengthList("LengthList");
    string MeanQuality("MeanQuanlity");
    string AGCTContent("AGCTContent");
    string GCContent("GCContent");
    string NContent("NContent");
    string Overkmer("Overkmer");
    /*
     * 添加元素
     */
    outhtml.append(HTMLHeader());
    outhtml.append(HTMLCss());
    //Basic Status
    outhtml.append(insertTableBegin());
    outhtml.append(insertTableTitle("Measure","Value"));
    outhtml.append(insertTableTbobyBegin());
    outhtml.append(insertTableTr("FileName",filename));
    outhtml.append(insertTableTr("File Type","Conventional base calls"));
    outhtml.append(insertTableTr("Encoding","Sanger / Illumina 1.9"));
    outhtml.append(insertTableTr("Total Sequence",to_string(this->total_number)));
    outhtml.append(insertTableTr("Sequence flagged as poor quality","0"));
    outhtml.append(insertTableTr("Sequence Length",to_string(min_len)+"-"+to_string(max_len)));
    long long total_AGCT=0;long long total_GC=0;
    for (int i=0;i<max_len;i++){
        for (int j=0;j<8;j++) total_AGCT+=NumberList[i][j];
        total_GC+=NumberList[i]['G'&0x07];
        total_GC+=NumberList[i]['C'&0x07];
    }
    outhtml.append(insertTableTr("%GC",to_string((100*total_GC)/total_AGCT)));
    outhtml.append(insertTableTbobyEnd());
    outhtml.append(insertTableEnd());

    // Quality Scores cross all bases
    outhtml.append(insertDiv(QualityScore));

    // lengthlist
    outhtml.append(insertDiv(LengthList));

    //Mean Quanlity
    outhtml.append(insertDiv(MeanQuality));

    //AGCT Content
    outhtml.append(insertDiv(AGCTContent));

    //GC Content
    outhtml.append(insertDiv(GCContent));

    //N Content
    outhtml.append(insertDiv(NContent));

    //Over Kmer
    outhtml.append(insertDiv(Overkmer));

    outhtml.append("</body>\n");

    /*
     *  添加JS
     */
    outhtml.append("<script type=\"text/javascript\">\n");

    // Quality Scores cross all bases
    outhtml.append(insertChart(QualityScore));
    //option
    outhtml.append(insertOptionBegin(QualityScore));
    outhtml.append(insertTitle("Quality Scores cross all bases"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Position in read(bp)",max_len));
    outhtml.append(insertyAxis("value",to_string(0),to_string(41)));
    outhtml.append(insertSeriesBegin());
    for (int i=0;i<max_len;i++){
        long long total_q=0;
        long long total_n=0;
        for (int j=0;j<8;j++) total_n+=NumberList[i][j];
        for (int j=0;j<QulityLen;j++) total_q+=(long long)j*QualityPositonList[i][j];

        tmp_double[i]=(total_q+0.0)/total_n;
        //printf("%lld %lld\n",total_q,total_n);
    }
    outhtml.append(insertSeriesData("line",tmp_double,max_len));

    outhtml.append(insertSeriesMultiDataBegin("boxplot"));
    for (int i=0;i<max_len;i++){
        int max_q=0;
        int min_q=50;
        int q_25=0;
        int q_50=0;
        int q_75=0;
        long long total_n=0;
        long long total_now=0;
        for (int j=0;j<8;j++) total_n+=NumberList[i][j];
        for (int j=0;j<QulityLen;j++){
            if (QualityPositonList[i][j]){
                min_q=min(min_q,j);
                max_q=max(max_q,j);
                total_now+=QualityPositonList[i][j];
                if (!q_25 && total_now*4>=total_n) q_25=j;
                if (!q_50 && total_now*2>=total_n) q_50=j;
                if (!q_75 && total_now*4>=3*total_n) q_75=j;
            }
        }
        tmp_int[0]=min_q;
        tmp_int[1]=q_25;
        tmp_int[2]=q_50;
        tmp_int[3]=q_75;
        tmp_int[4]=max_q;
        outhtml.append(insertSeriesOneData(tmp_int,5));
    }
    outhtml.append(insertSeriesMultiDataEnd());
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(QualityScore));






    //lengthlist
    outhtml.append(insertChart(LengthList));
    //option
    outhtml.append(insertOptionBegin(LengthList));
    outhtml.append(insertTitle("Seqence Length List"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Sequence Length(bp)",max_len+1));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    outhtml.append(insertSeriesData("line",LengthSequence,max_len+1));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(LengthList));


    // Mean Quanlity
    outhtml.append(insertChart(MeanQuality));
    //option
    outhtml.append(insertOptionBegin(MeanQuality));
    outhtml.append(insertTitle("Mean Quanlity List"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Mean Sequence Quality(Phred Score)",42));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    outhtml.append(insertSeriesData("line",QualitySequence,42));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(MeanQuality));


    //AGCT Content
    outhtml.append(insertChart(AGCTContent));
    //option
    outhtml.append(insertOptionBegin(AGCTContent));
    outhtml.append(insertTitle("AGCT Content"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertLegend("\'A\',\'G\',\'C\',\'T\'"));
    outhtml.append(insertxAxis("Position in read(bp)",max_len));
    outhtml.append(insertyAxis("value",to_string(0),to_string(1)));
    outhtml.append(insertSeriesBegin());


    for (int i=0;i<max_len;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['A'&0x07]/total;
    }
    outhtml.append(insertSeriesData("line","A",tmp_double,max_len));
    for (int i=0;i<max_len;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['G'&0x07]/total;
    }
    outhtml.append(insertSeriesData("line","G",tmp_double,max_len));
    for (int i=0;i<max_len;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['C'&0x07]/total;
    }
    outhtml.append(insertSeriesData("line","C",tmp_double,max_len));
    for (int i=0;i<max_len;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['T'&0x07]/total;
    }

    outhtml.append(insertSeriesData("line","T",tmp_double,max_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(AGCTContent));

    //GC Content
    outhtml.append(insertChart(GCContent));
    //option
    int GCinterval=10;
    outhtml.append(insertOptionBegin(GCContent));
    outhtml.append(insertTitle("Per Sequence GC content"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Mean GC content(%)",101,1));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    memset(tmp_int,0x00,MAXLEN*sizeof(int));
    int last=0;
    int lastpos=-1;
    for (int i=0;i<ContentLen;i+=1){
        if (Content[i][('G'+'C')&0x07]){
            tmp_int[i]=Content[i][('G'+'C')&0x07];
            if (lastpos+1!=i){
                for (int k=lastpos+1;k<i;k++){
                    tmp_int[k]=last+(tmp_int[i]-last)*(k-lastpos)/(i-lastpos);
                }
            }
            lastpos=i;
            last=tmp_int[i];
        }else{
            if (i==ContentLen-1){
                for (int k=lastpos+1;k<i;k++){
                    tmp_int[k]=last+(tmp_int[i]-last)*(k-lastpos)/(k-lastpos);
                }
            }
        }
    }
    outhtml.append(insertSeriesData("line",tmp_int,ContentLen,GCinterval));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(GCContent));

    //N Content
    outhtml.append(insertChart(NContent));
    //option
    outhtml.append(insertOptionBegin(NContent));
    outhtml.append(insertTitle("Per Sequence N content"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("position in read(%)",max_len));
    outhtml.append(insertyAxis("value",to_string(0),to_string(1)));
    outhtml.append(insertSeriesBegin());
    for (int i=0;i<max_len+1;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['N'&0x07]/total;
    }
    outhtml.append(insertSeriesData("line",tmp_double,max_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(NContent));

    //Over Kmer
    outhtml.append(insertChart(Overkmer));
    //option
    outhtml.append(insertOptionBegin(Overkmer));
    outhtml.append(insertTitle("Kmer"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Position in read(bp)",max_len));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    outhtml.append(insertSeriesData("line",Kmer[0],max_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(Overkmer));



    outhtml.append("</script>");
    outhtml.append("</html>");
    fout->write(outhtml.c_str(),outhtml.length());
}
void BamStatus::reportHTML(ofstream *fout,Duplicate *duplicate,Overrepresent *overrepresent){
    string outhtml;

    int *tmp_int=new int[MAXLEN];
    double *tmp_double=new double[MAXLEN];
    char *tmp_char=new char[MAXLEN];

    string QualityScore("QualityScore");
    string LengthList("LengthList");
    string MeanQuality("MeanQuanlity");
    string AGCTContent("AGCTContent");
    string GCContent("GCContent");
    string NContent("NContent");
    string DuplicatePercent("DuplicatePercent");
    string Overkmer("Overkmer");
    /*
     * 添加元素
     */
    outhtml.append(HTMLHeader());
    outhtml.append(HTMLCss());
    //Basic Status
    outhtml.append(insertTableBegin());
    outhtml.append(insertTableTitle("Measure","Value"));
    outhtml.append(insertTableTbobyBegin());
    outhtml.append(insertTableTr("FileName",filename));
    outhtml.append(insertTableTr("File Type","Conventional base calls"));
    outhtml.append(insertTableTr("Encoding","Sanger / Illumina 1.9"));
    outhtml.append(insertTableTr("Total Sequence",to_string(this->total_number)));
    outhtml.append(insertTableTr("Sequence flagged as poor quality","0"));
    outhtml.append(insertTableTr("Sequence Length",to_string(min_len)+"-"+to_string(max_len)));
    long long total_AGCT=0;long long total_GC=0;
    for (int i=0;i<max_len;i++){
        for (int j=0;j<8;j++) total_AGCT+=NumberList[i][j];
        total_GC+=NumberList[i]['G'&0x07];
        total_GC+=NumberList[i]['C'&0x07];
    }
    outhtml.append(insertTableTr("%GC",to_string((100*total_GC)/total_AGCT)));
    outhtml.append(insertTableTr("Over Represent Date",to_string(overrepresent->OverrepresentDate*100)+"%"));
    outhtml.append(insertTableTbobyEnd());
    outhtml.append(insertTableEnd());

    // Quality Scores cross all bases
    outhtml.append(insertDiv(QualityScore));

    // lengthlist
    outhtml.append(insertDiv(LengthList));

    //Mean Quanlity
    outhtml.append(insertDiv(MeanQuality));

    //AGCT Content
    outhtml.append(insertDiv(AGCTContent));

    //GC Content
    outhtml.append(insertDiv(GCContent));

    //N Content
    outhtml.append(insertDiv(NContent));

    //Duplicate Percent
    outhtml.append(insertDiv(DuplicatePercent));

    //Over Kmer
    outhtml.append(insertDiv(Overkmer));

    outhtml.append("</body>\n");

    /*
     *  添加JS
     */
    outhtml.append("<script type=\"text/javascript\">\n");

    // Quality Scores cross all bases
    outhtml.append(insertChart(QualityScore));
    //option
    outhtml.append(insertOptionBegin(QualityScore));
    outhtml.append(insertTitle("Quality Scores cross all bases"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Position in read(bp)",max_len));
    outhtml.append(insertyAxis("value",to_string(0),to_string(41)));
    outhtml.append(insertSeriesBegin());
    for (int i=0;i<max_len;i++){
        long long total_q=0;
        long long total_n=0;
        for (int j=0;j<8;j++) total_n+=NumberList[i][j];
        for (int j=0;j<QulityLen;j++) total_q+=(long long)j*QualityPositonList[i][j];

        tmp_double[i]=(total_q+0.0)/total_n;
        //printf("%lld %lld\n",total_q,total_n);
    }
    outhtml.append(insertSeriesData("line",tmp_double,max_len));

    outhtml.append(insertSeriesMultiDataBegin("boxplot"));
    for (int i=0;i<max_len;i++){
        int max_q=0;
        int min_q=50;
        int q_25=0;
        int q_50=0;
        int q_75=0;
        long long total_n=0;
        long long total_now=0;
        for (int j=0;j<8;j++) total_n+=NumberList[i][j];
        for (int j=0;j<QulityLen;j++){
            if (QualityPositonList[i][j]){
                min_q=min(min_q,j);
                max_q=max(max_q,j);
                total_now+=QualityPositonList[i][j];
                if (!q_25 && total_now*4>=total_n) q_25=j;
                if (!q_50 && total_now*2>=total_n) q_50=j;
                if (!q_75 && total_now*4>=3*total_n) q_75=j;
            }
        }
        tmp_int[0]=min_q;
        tmp_int[1]=q_25;
        tmp_int[2]=q_50;
        tmp_int[3]=q_75;
        tmp_int[4]=max_q;
        outhtml.append(insertSeriesOneData(tmp_int,5));
    }
    outhtml.append(insertSeriesMultiDataEnd());
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(QualityScore));






    //lengthlist
    outhtml.append(insertChart(LengthList));
    //option
    outhtml.append(insertOptionBegin(LengthList));
    outhtml.append(insertTitle("Seqence Length List"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Sequence Length(bp)",max_len+1));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    outhtml.append(insertSeriesData("line",LengthSequence,max_len+1));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(LengthList));


    // Mean Quanlity
    outhtml.append(insertChart(MeanQuality));
    //option
    outhtml.append(insertOptionBegin(MeanQuality));
    outhtml.append(insertTitle("Mean Quanlity List"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Mean Sequence Quality(Phred Score)",42));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    outhtml.append(insertSeriesData("line",QualitySequence,42));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(MeanQuality));


    //AGCT Content
    outhtml.append(insertChart(AGCTContent));
    //option
    outhtml.append(insertOptionBegin(AGCTContent));
    outhtml.append(insertTitle("AGCT Content"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertLegend("\'A\',\'G\',\'C\',\'T\'"));
    outhtml.append(insertxAxis("Position in read(bp)",max_len));
    outhtml.append(insertyAxis("value",to_string(0),to_string(1)));
    outhtml.append(insertSeriesBegin());


    for (int i=0;i<max_len;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['A'&0x07]/total;
    }
    outhtml.append(insertSeriesData("line","A",tmp_double,max_len));
    for (int i=0;i<max_len;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['G'&0x07]/total;
    }
    outhtml.append(insertSeriesData("line","G",tmp_double,max_len));
    for (int i=0;i<max_len;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['C'&0x07]/total;
    }
    outhtml.append(insertSeriesData("line","C",tmp_double,max_len));
    for (int i=0;i<max_len;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['T'&0x07]/total;
    }

    outhtml.append(insertSeriesData("line","T",tmp_double,max_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(AGCTContent));

    //GC Content
    outhtml.append(insertChart(GCContent));
    //option
    int GCinterval=10;
    outhtml.append(insertOptionBegin(GCContent));
    outhtml.append(insertTitle("Per Sequence GC content"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Mean GC content(%)",101,1));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    memset(tmp_int,0x00,MAXLEN*sizeof(int));
    int last=0;
    int lastpos=-1;
    for (int i=0;i<ContentLen;i+=1){
        if (Content[i][('G'+'C')&0x07]){
            tmp_int[i]=Content[i][('G'+'C')&0x07];
            if (lastpos+1!=i){
                for (int k=lastpos+1;k<i;k++){
                    tmp_int[k]=last+(tmp_int[i]-last)*(k-lastpos)/(i-lastpos);
                }
            }
            lastpos=i;
            last=tmp_int[i];
        }else{
            if (i==ContentLen-1){
                for (int k=lastpos+1;k<i;k++){
                    tmp_int[k]=last+(tmp_int[i]-last)*(k-lastpos)/(k-lastpos);
                }
            }
        }
    }
    outhtml.append(insertSeriesData("line",tmp_int,ContentLen,GCinterval));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(GCContent));

    //N Content
    outhtml.append(insertChart(NContent));
    //option
    outhtml.append(insertOptionBegin(NContent));
    outhtml.append(insertTitle("Per Sequence N content"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("position in read(%)",max_len));
    outhtml.append(insertyAxis("value",to_string(0),to_string(1)));
    outhtml.append(insertSeriesBegin());
    for (int i=0;i<max_len+1;i++){
        double total=0.0;
        for (int j=0;j<8;j++) total+=NumberList[i][j];
        tmp_double[i]=NumberList[i]['N'&0x07]/total;
    }
    outhtml.append(insertSeriesData("line",tmp_double,max_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(NContent));

    //Duplicate Percent
    outhtml.append(insertChart(DuplicatePercent));
    //option
    outhtml.append(insertOptionBegin(DuplicatePercent));
    int* dupHist = NULL;
    double* dupMeanTlen = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    int histSize=32;
    dupHist = new int[histSize];
    memset(dupHist, 0, sizeof(int) * histSize);
    dupMeanGC = new double[histSize];
    memset(dupMeanGC, 0, sizeof(double) * histSize);
    dupRate = duplicate->statAll(dupHist, dupMeanGC, histSize);
    outhtml.append(insertTitle("duplicate rate "+to_string(dupRate*100.0)+ "%"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Sequence Duplication level",histSize));
    outhtml.append(insertyAxis("value",to_string(0),to_string(1)));
    outhtml.append(insertSeriesBegin());
    double dup_total=0.0;
    for (int i=0;i<histSize;i++) {
        dup_total+=dupHist[i];
    }
    for (int i=0;i<histSize;i++){
        tmp_double[i]=(dupHist[i]+0.0)/dup_total;
    }
    outhtml.append(insertSeriesData("bar",tmp_double,histSize));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(DuplicatePercent));

    //Over Kmer
    outhtml.append(insertChart(Overkmer));
    //option
    outhtml.append(insertOptionBegin(Overkmer));
    outhtml.append(insertTitle("Relative enrichment over read length"));
    outhtml.append(insertTooltip());
    int tmp_char_pos=0;
    for (int i=0;i<ChooseKmerNum;i++) {
        tmp_char[tmp_char_pos++]='\'';
        for (int k = 0, x = ChooseKmerPos[i]; k < KmerBase; k++) {
            tmp_char[tmp_char_pos++]=charAGCTN[x&0x03];
            x >>= 2;
        }
        tmp_char[tmp_char_pos++] = '\'';
        tmp_char[tmp_char_pos++] = ',';
    }
    tmp_char[tmp_char_pos++]='\0';
    outhtml.append(insertLegend(string(tmp_char)));
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Position in read(bp)",max_len));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    for (int i=0;i<ChooseKmerNum;i++){
        for(int k=0,x=ChooseKmerPos[i];k<KmerBase;k++){
            tmp_char[k]=charAGCTN[x&0x03];
            x>>=2;
        }
        tmp_char[KmerBase]='\0';
        outhtml.append(insertSeriesData("line",string(tmp_char),Kmer[ChooseKmerPos[i]],max_len));
    }
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(Overkmer));



    outhtml.append("</script>");
    outhtml.append("</html>");
    fout->write(outhtml.c_str(),outhtml.length());
}